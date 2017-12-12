'''
Module contains routines to perform grid optimization with two options:
  1. lsqr without constraint to solve Ax=b (scipy.sparse.linalg.lsqr) 
  2.  minimize function 1/2*||Ax-b||^2 using the L-BFGS-B algorithm (scipy.optimize.fmin_l_bfgs_b)
  
x is depth at nodes with respect to a nominal reference surface (which may be greater than sea level in upstream locations).

Regularization:

  1.  close to the original values for all nodes: damp 
  2.  minimize the 1st derivative of node elevation along land boundaries: damp_shoreline

Notes:

 * The function "cal_z_ref" only applies to the Bay Delta grid and need to be modified for other grids.
 * This script create an optimized gr3 file (_opt.gr3) only when one set of optimization parameters is specified.

First version, completed on 9/11/2013

'''

import numpy as np
from scipy.sparse import lil_matrix, vstack
from scipy.sparse.linalg import lsqr
from scipy.optimize import fmin_l_bfgs_b
import stacked_dem_fill as sdf
import os.path
import osgeo.ogr, osgeo.osr
from osgeo import ogr
import math
import matplotlib.pyplot as plt

def cal_length(a):
    return np.sqrt((a[0]-a[2])**2 + (a[1]-a[3])**2)
    
def cal_tri_area(a):
    return np.absolute((a[0]*(a[3]-a[5])+a[2]*(a[5]-a[1])+a[4]*(a[1]-a[3]))/2.0)

def cal_distance(x1,y1,x2,y2):
    return math.sqrt((x1-x2)**2+(y1-y2)**2)     
    
def cal_z_ref(xy_coor):
    '''
    Define reference water surface at locations specified in xy_coor
    based on the assumptions, which only appliable to Bay Delta grid
    (1) zref increases linerly from west (ocean) to east
    (2) east of (x_old_river, y_old_river), zref increases lineraly towards the south
    '''
    x_east = 653768.1
    x_west = 499546.3
    x_old_river = 647170.0  
    y_old_river = 4186000.0
    y_sj_end = 4171563.0
    z_ref = np.where(xy_coor[:,0]>x_old_river,3.0*(y_old_river-xy_coor[:,1])/(y_old_river-y_sj_end),0.0)
    z_ref = z_ref.clip(0.0)
    z_ref = z_ref + 0.5*(xy_coor[:,0]-x_west)/(x_east-x_west)
    return z_ref
    
def create_ref_surf(elements, nodes):    
    '''
    create reference surface at individual nodes:
    step 1: derive max levation within connected elements
    setp 2: compare with zref calculated by cal_z_ref() and use the higher value between the two
    '''
    nelm = len(elements)
    nnode = len(nodes)
    node_map = {i:[] for i in range(nnode)}
    for i in range(nelm):
        for j in range(3):
            node_map[elements[i][j]].append(i)
    elv_list = []
    ref_elv = []
    for i in range(nnode):
        for elm in node_map[i]:
            for node in elements[elm]:
                elv_list.append(nodes[node][2])
        ref_elv.append(max(elv_list))
        elv_list = []                    
    z1 = cal_z_ref(nodes)
    z2 = np.array(ref_elv, dtype='float64')
    ref_elv = np.maximum(z1,z2)
    return ref_elv
    
def construct_A_and_b(damp,damp_bnd,face_coeff,vol_coeff,A_tri,A_face,A_bnd,depth_q,face_hq,bnd_h,nnode,h0_node,solver='L-BFGS-B'):
    '''
    construct matrix A and vector b based on specified optimization parameters 
    damp: regularization to initial values
    A_bnd, damp_bnd: regularization to minimize gradient along the shoreline
    A_face, face_coeff: related to vertical faces of elements
    A_tri, vol_coeff: related to volumes of elements
    depth_q, face_hq: quadrature depths derived from elements and and faces
    bnd_h: target depth along the shoreline
    h0_node: original depth at nodes
    nnode: node number
    solver: either L-BFGS-B or lsqr
    '''  
    if  (face_coeff == 0 and vol_coeff == 0) or face_coeff < 0 or vol_coeff <0:
        print 'invalid vol_coeff and face_coeff, default values (1, 1) will be used'
        face_coeff = 1
        vol_coeff = 1   
    if solver != 'L-BFGS-B' and solver != 'lsqr':
        solver = 'L-BFGS-B'
        print 'invalid solver specified, L-BFGS-B will be used'
    if solver == 'lsqr': 
        if damp_bnd > 0.0:
            A = vstack([vol_coeff*A_tri,face_coeff*A_face,damp_bnd*A_bnd])              
            b = np.hstack((vol_coeff*(depth_q),face_coeff*(face_hq),damp_bnd*bnd_h))
        else:        
            A = vstack([vol_coeff*A_tri,face_coeff*A_face])              
            b = np.hstack((vol_coeff*(depth_q),face_coeff*(face_hq)))       
    else: 
        if damp > 0.0:
            A_node = lil_matrix((nnode,nnode))
            for i in range(nnode):
                A_node [i,i] = 1
        if damp_bnd > 0.0 and damp > 0.0 :
            A = vstack([vol_coeff*A_tri,face_coeff*A_face,damp_bnd*A_bnd,damp*A_node])              
            b = np.hstack((vol_coeff*(depth_q),face_coeff*(face_hq),damp_bnd*bnd_h,damp*h0_node))
        elif damp_bnd <= 0.0 and damp > 0.0:    
            A = vstack([vol_coeff*A_tri,face_coeff*A_face,damp*A_node])              
            b = np.hstack((vol_coeff*(depth_q),face_coeff*(face_hq),damp*h0_node))
        elif damp_bnd > 0.0 and damp <= 0.0 :
            A = vstack([vol_coeff*A_tri,face_coeff*A_face,damp_bnd*A_bnd])              
            b = np.hstack((vol_coeff*(depth_q),face_coeff*(face_hq),damp_bnd*bnd_h))
        else:        
            A = vstack([vol_coeff*A_tri,face_coeff*A_face])              
            b = np.hstack((vol_coeff*(depth_q),face_coeff*(face_hq)))
    return (A, b)    
    
def obj_func(x,A,b):  
    '''
    define the function to be minimize using L-BFGS-B algorithm
    1/2*||Ax-b||^2 
    '''
    y = A*x-b
    return 0.5*y.dot(y)
    
def fprime(x,A,b):
    '''
    The gradient of the objective function with respect to the heights (x). 
    This gradient is a vector the same size as x. 
    A^T (Ax-b)
    '''
    return A.transpose()*(A*x-b)    

def depth_opt(damp,damp_bnd,face_coeff,vol_coeff,A_tri,A_face,A_bnd,depth_q,face_hq,bnd_h,nnode,h0_node,solver='L-BFGS-B'):
    if solver != 'L-BFGS-B' and solver != 'lsqr':
        solver = 'L-BFGS-B'
        print 'invalid solver specified, L-BFGS-B will be used'
    print solver + ' optimization parameters used: '+str(damp)+', '+str(damp_bnd)+', '+str(face_coeff)+', '+str(vol_coeff)
    (A, b) = construct_A_and_b(damp,damp_bnd,face_coeff,vol_coeff,A_tri,A_face,A_bnd,depth_q,face_hq,bnd_h,nnode,h0_node,solver)    
    if solver == 'lsqr': 
        result = lsqr(A, b, damp)  
    else:
        result = fmin_l_bfgs_b(obj_func,h0_node,fprime=fprime,args=(A,b),approx_grad=0,bounds=[(0.,None)]*nnode, m=10, 
            factr=100.0, pgtol=1e-05, epsilon=1e-08, iprint=-1, maxfun=15000, maxiter=15000, disp=None, callback=None)
    print result  
    return result[0]
    
def grid_opt(nodes,elements,land_boundaries,dem_list_file,opt_parms,output_files=False,outfile_prefix='',remove_boundary_faces=False, solver='L-BFGS-B'):
    '''
    Perform grid optimization with two solvers:
    (1) L-BFGS-B (default): minimize function 1/2*||Ax-b||^2
    (2) linear least square solver (Ax=b). 
    nodes:  (x,y,z) for each node, numpy array, z is elevation (a place holder, valuse not used)
    elements: node numbers (3) for each element, numpy array
    land_boundaries: (land boubdary node pair index, node 1, node 2, segment no, sequence in segment), numpy array 
    dem_list_file:  a text file, containing a list of DEM used to derive elevation
    opt_parms: [[damp, damp_shoreline, face_coeff, volume_coeff]], list
    output_files: True or False, whether to create output files
    outfile_prefix: Prefix used for names of output files               
    '''
#===================================================================================
#Elements -- Define the evaluation point locations (in Barycentric coordinate) and 
#            weights for intergation on triangles
#===================================================================================
    nvert = 3 
    nquad = 7
    quadpts = np.zeros((nquad,3))
    quadpts = np.array([[0.333333333333333,0.333333333333333,0.333333333333333],
                        [0.797426985353087,0.101286507323456,0.101286507323456],
                        [0.101286507323456,0.797426985353087,0.101286507323456],
                        [0.101286507323456,0.101286507323456,0.797426985353087],
                        [0.059715871789770,0.470142064105115,0.470142064105115],
                        [0.470142064105115,0.059715871789770,0.470142064105115],
                        [0.470142064105115,0.470142064105115,0.059715871789770]])
                   
    quadwts = np.array([0.225,0.125939180544827,0.125939180544827,0.125939180544827,0.132394152788506,
                        0.132394152788506,0.132394152788506])
#==============================================================
#Vertical Faces -- 1D gaussian quadrature
#==============================================================
    nquad_1d = 5
    quadpts_1d = np.array([0.0,0.538469310105683,-0.538469310105683,0.906179845938664,-0.906179845938664])
    quadwts_1d = np.array([0.568888888888889,0.478628670499366,0.478628670499366,0.236926885056189,0.236926885056189])  

    nelement = len(elements)
    nnode = len(nodes)
    print "N nodes: %s\n" % nnode
    print "N element: %s\n" % nelement

#set z information to null and derive it from DEMs
    nodes[:,2] = np.nan   
    files = [x.strip() for x in open(dem_list_file,"r").readlines() if x and len(x) > 1 and not x.startswith("#")] 
    print 'Assign elevation to nodes......'
    nodes = sdf.stacked_dem_fill(files,nodes,require_all = False,na_fill = 2)
        
#create reference surface
    ref_surf = create_ref_surf(elements, nodes)  
    nodes_coor = np.hstack((nodes,ref_surf.reshape(nnode,1)))
    #np.savetxt('nodes_coor.txt',nodes_coor)
#calculate initial depth at each node (ref elev - node elev)
    h0_node = nodes_coor[:,3] - nodes_coor[:,2]
    
    if (np.min(h0_node) < 0.0):
        exit("Error: Some initial depth at nodes is less than 0, indicating somthing went wrong. Abort.....")

#populate (x,y) and elevation information for each triangular element
    tri = np.zeros((nelement,nvert,4))
    for itri in range(nelement):
        for ivert in range(nvert):
            tri[itri,ivert,:] = nodes_coor[elements[itri][ivert]]  

#calculate the initial averaged depth for individual trianglar element 
    depth_0 = tri[:,:,3]-tri[:,:,2]    
    depth_0 = np.average(depth_0, axis=1)    
    if (np.min(depth_0) < 0.0):
        exit("Error: Some initial averaged depth for elements is less than 0, indicating something went wrong.  Abort.....")
  
#derive (x,y) locations of the 7 quadrature points for volume integration 
    x_coor = np.dot(quadpts,tri[:,:,0].transpose()) #x-coor for the 7 pts for all tri, 7 by nelement
    y_coor = np.dot(quadpts,tri[:,:,1].transpose()) #y-coor for the 7 pts for all tri, 7 by nelement
    x_coor = x_coor.transpose().flatten()
    y_coor = y_coor.transpose().flatten()
    tri_pts = np.column_stack((x_coor,y_coor,np.zeros(x_coor.size)))
    tri_pts[:,2] = np.nan

#extract face information from elements
    tmp = [(elements[i][0],elements[i][1]) if elements[i][0]<=elements[i][1] else (elements[i][1],elements[i][0])
            for i in range (0,nelement)]
    tmp = tmp + [ (elements[i][1],elements[i][2]) if elements[i][1]<=elements[i][2] else (elements[i][2],elements[i][1])
            for i in range (0,nelement)]
    tmp = tmp + [ (elements[i][2],elements[i][0]) if elements[i][2]<=elements[i][0] else (elements[i][0],elements[i][2])
            for i in range (0,nelement)]
    tmp = set(tmp)
    face_nodes = np.array([ s for s in tmp ], dtype='int')
    nface = len(face_nodes)
    del tmp
#remove faces along the land boundaries 
    nbnd = len(land_boundaries)
    if nbnd > 0 and remove_boundary_faces == True:
        tmp1 = np.array([str(x[0])+str(x[1]) for x in face_nodes], dtype = 'int64')
        bnd_node_pairs = [(land_boundaries[i][1],land_boundaries[i][2]) if (land_boundaries[i][1]<=land_boundaries[i][2]) 
                                else (land_boundaries[i][2],land_boundaries[i][1]) for i in range(nbnd)]
        tmp2 = np.array([str(x[0])+str(x[1]) for x in bnd_node_pairs], dtype = 'int64')
        face_to_del = np.in1d(tmp1,tmp2)
        print ' Remove ' + str(len(face_to_del[face_to_del==True])) + ' faces along land boundaries'
        face_to_del=np.where(face_to_del==True)
        face_nodes = np.delete(face_nodes,face_to_del,0) 
        nface = len(face_nodes)
        del bnd_node_pairs
        del tmp1
        del tmp2
        del face_to_del

#populate (x,y) and elevation information for each vertical face
    face = np.zeros((nface,2,4))
    for iface in range(nface):
        for ivert in range(2):
            face[iface,ivert,:] = nodes_coor[face_nodes[iface][ivert]]  
        
#calculate the average depth based on two nodes for each face.
    face_h0 = face[:,:,3] - face[:,:,2]
    face_h0 = np.average(face_h0, axis=1)
    if (np.min(face_h0) < 0.0): 
        exit("Error: Some original depth at face is less than 0, indicating somthing went wrong. Abort......")
    
#derive (x,y) locations of the 5 points (1D Gaussian Quadrature)
    x_coor = 0.5*(face[:,1,0]+face[:,0,0]).reshape(nface,1) + 0.5*quadpts_1d*(face[:,1,0]-face[:,0,0]).reshape(nface,1) #x-coor for the 5 pts for all face, nface by 5
    y_coor = 0.5*(face[:,1,1]+face[:,0,1]).reshape(nface,1) + 0.5*quadpts_1d*(face[:,1,1]-face[:,0,1]).reshape(nface,1) #x-coor for the 5 pts for all face, nface by 5
    face_pts = np.column_stack((x_coor.flatten(),y_coor.flatten(),np.zeros(nface*nquad_1d)))
    face_pts[:,2] = np.nan
    del x_coor
    del y_coor

#calculate elevations for quadrature points (nquad*nelement) for volume and (nquad_1d*nface) for face
    points = np.vstack((tri_pts,face_pts))
    print 'Assign elevation to quadrature points for elements and faces......'
    points = sdf.stacked_dem_fill(files,points,require_all = False,na_fill = 2)
    tri_pts = points[0:nelement*nquad,:]
    face_pts = points[nelement*nquad:,:]
    #np.savetxt('tri_pts.txt',tri_pts)
    #np.savetxt('face_pts.txt',face_pts)
    del points

#calculate quadrature depth for each element by numerical integration    
    depth_dem = (np.dot(tri[:,:,3],quadpts.transpose()).reshape(nelement*7,) - tri_pts[:,2]).reshape(nelement,7)
    depth_less_zero = np.where(depth_dem<0)
    #np.savetxt('depth_dem.txt',np.array(depth_dem))
    depth_less_zero = set(depth_less_zero[0])
    if (len(depth_less_zero) > 1):
        print 'Warning: '+ str(len(depth_less_zero)) + ' elements containing quadrature pts with negative depth; set to zero' 
#set individual point depth to 0.0 if negative (i.e. above the reference surface)
    np.clip(depth_dem,0.0,99999,depth_dem) 
    depth_q = np.sum(depth_dem*quadwts, axis=1)
    del depth_dem 
#calculate the area for each element
    tri_areas = cal_tri_area(tri[:,:,0:2].reshape(nelement,6).transpose())

#calculate quadrature depth for each face 
    face_ref = 0.5*(face[:,1,3]+face[:,0,3]).reshape(nface,1) + 0.5*quadpts_1d*(face[:,1,3]-face[:,0,3]).reshape(nface,1)
    face_hq = face_ref.reshape(nface*nquad_1d) - face_pts[:,2]
    depth_less_zero = np.where(face_hq.reshape(nface,nquad_1d)<0)
    #np.savetxt('face_hq.txt',face_hq)
    depth_less_zero = set(depth_less_zero[0])
    if (len(depth_less_zero) > 1):
        print 'Warning: '+ str(len(depth_less_zero)) + ' faces containing quadrature pts with negative depth; set to zero' 
#set the depth at each node to 0.0 when negative before taking the average
    face_hq = face_hq.clip(0.0)
    face_hq = np.sum(face_hq.reshape(nface,nquad_1d)*quadwts_1d, axis=1)/2.0
    del face_ref
#calculate the length for each face
    face_len = cal_length(face[:,:,:2].reshape(nface,4).transpose())

#tempory testing script for excluding discontinuities, using cached data from process_bnd.py        
#    land_boundaries = np.loadtxt('D:/grid_chk/land_boundaries73.txt')
#    only include non-sms nodes boundary nodes in regularizations (i.e. land_boundaries[:][5] = 0) :
#    when land_boundaries[i-1][5] +land_boundaries[i][5]+ land_boundaries[i+1][5] == 0
#    results worse than original, drop the attempt
     
#assume two adjecent nodes along the land boundary have the same elevation, Zi=Zj
#the target value = (Zrefi-Zi)-(Zrefj-Zj) = Zrefi - Zrefj
#    bnd_h = np.array([nodes_coor[x[1]][3] - nodes_coor[x[2]][3] for x in land_boundaries])
#assume two adjecent nodes along the land boundary have the same 1st derivates at adjacent nodes:
#(Zi - Zi-1)/Li-1 = (Zi+1 - Zi)/Li+1
#if a = Li-1/Li+1, we have (1+a)*hi - a*hi+1 - hi-1 = (1+a)* Zrefi - a* Zrefi+1 - Zrefi-1
    bnd_seg_no = set(land_boundaries[:,3]) # find the number of land boundary segments
    bnd_eq = []   
    for i in bnd_seg_no:
        tmp = np.where(land_boundaries[:,3]==i)[0]
        if len(tmp) > 2:
            for i in tmp[1:-1]:
                node_center = land_boundaries[i][1]
                node_plus = land_boundaries[i+1][1]
                node_minus = land_boundaries[i-1][1]
                a = cal_distance(nodes_coor[node_center][0],nodes_coor[node_center][1],
                         nodes_coor[node_minus][0],nodes_coor[node_minus][1]) /      \
                    cal_distance(nodes_coor[node_center][0],nodes_coor[node_center][1],
                         nodes_coor[node_plus][0],nodes_coor[node_plus][1])  
                bnd_eq.append((node_center,node_plus,node_minus,a))                     
    bnd_h = np.zeros(len(bnd_eq))  
    for i in range(len(bnd_eq)):
        Zref_center = nodes_coor[bnd_eq[i][0]][3]
        Zref_plus = nodes_coor[bnd_eq[i][1]][3]  
        Zref_minus = nodes_coor[bnd_eq[i][2]][3]
        a = bnd_eq[i][3]
        if (Zref_center == Zref_plus) and (Zref_center == Zref_minus):
            bnd_h[i] = 0.
        elif (Zref_center == Zref_plus):
            bnd_h[i] = Zref_center - Zref_minus
        elif (Zref_center == Zref_minus):
            bnd_h[i] = a*(Zref_center - Zref_plus)
        else:   
            bnd_h[i] = Zref_center*(1.0+a) - a*Zref_plus - Zref_minus
      
#prepare components for matrix A
    A_tri = lil_matrix((nelement,nnode))
    for itri in range(nelement):
        for ivert in range(nvert):
            A_tri[ itri, elements[itri][ivert] ] = 1.0/3.0
    A_face = lil_matrix((nface,nnode))  
    for iface in range(nface):
        for ivert in range(2):
            A_face[ iface, face_nodes[iface][ivert] ] = 0.5     
    
    if nbnd > 0:     
        A_bnd = lil_matrix((len(bnd_eq),nnode))
        for ibnd in range(len(bnd_eq)):
            a = bnd_eq[ibnd][3]
            A_bnd[ ibnd, bnd_eq[ibnd][0]] = 1.0 + a
            A_bnd[ ibnd, bnd_eq[ibnd][1]] = -1.0 * a
            A_bnd[ ibnd, bnd_eq[ibnd][2]] = -1.0
    else:
        print 'No information on the boundary nodes, set optimization coeff. along the shoreline to 0.'
        A_bnd = []
        opt_parms = [[ parm[0],0,parm[2],parm[3]] for parm in opt_parms]    

#validate optimization parameters 
    flag = opt_parm_check(opt_parms)
    if flag <> 0:
        print 'No valid optimization paramaters specified, default values (0.2, 2, 1, 1) will be used'
        opt_parms = [[0.2, 2, 1, 1]]              
   
    if solver != 'L-BFGS-B' and solver != 'lsqr':
        solver = 'L-BFGS-B'
        print 'invalid solver specified, L-BFGS-B will be used'
#loop through optimization parameters     
    optimized_z = []
    for parm in opt_parms:
        print 'damp=' + str(parm[0])+', damp_shoreline=' + str(parm[1]) \
                +', face_coeff=' + str(parm[2]) + ', volume_coeff=' + str(parm[3])
#        node_h_new_lsqr = depth_opt(parm[0],parm[1],parm[2],parm[3],A_tri,A_face,A_bnd,depth_q,face_hq,bnd_h,nnode,h0_node,'lsqr')
        node_h_new = depth_opt(parm[0],parm[1],parm[2],parm[3],A_tri,A_face,A_bnd,depth_q,face_hq,bnd_h,nnode,h0_node,solver)
#        plt.plot(node_h_new_lsqr, 'go', label ='lsqr')
#        plt.plot(node_h_new, 'bo', label='L-BFGS-B')
#        plt.grid(b=True, which='both', color='0.65',linestyle='-')
#        plt.legend()
#        plt.show()
#        diff = node_h_new_lsqr-node_h_new
#        print 'difference = ' + str(diff.dot(diff))
#        print 'difference min and max = ' + str(np.min(diff)) + ', '+ str(np.max(diff))
#        print 'difference mean and std = ' + str(np.mean(diff)) +', '+ str(np.std(diff)) 
        optimized_z.append(nodes_coor[:,3] - node_h_new)
#new depth for individual triangular elements         
        h_new = A_tri * node_h_new       
#new depth for individual faces
        face_h_new = A_face * node_h_new 
        errors_list = [str(np.sqrt(np.sum((depth_0-depth_q)**2))), str(np.sqrt(np.sum((h_new-depth_q)**2))), \
          str(np.sqrt(np.sum((face_h0-face_hq)**2))), str(np.sqrt(np.sum((face_h_new-face_hq)**2))) ]
        print 'L2 volume old = ' + errors_list[0]
        print 'L2 volume new = ' + errors_list[1]
        print 'L2 face old   = ' + errors_list[2]
        print 'L2 face new   = ' + errors_list[3]   
   
        if output_files:
            print 'Writing detailed optimization results to files .......'        
        #output error summary 
            damp_str = str(parm[0])
            damp_bnd_str = str(parm[1])
            vol_coeff_str = str(parm[2])
            face_coeff_str = str(parm[3])
            output_file_name = outfile_prefix+'elm_sum_damp'+ damp_str + '_'+ damp_bnd_str +'vol'+vol_coeff_str + 'face'+ face_coeff_str + '.txt'
            f = open(output_file_name,'w')
            f.write('L2 volume old = ' + errors_list[0] +'\n')
            f.write('L2 volume new = ' + errors_list[1] +'\n')
            f.write('L2 face old   = ' + errors_list[2] +'\n')
            f.write('L2 face new   = ' + errors_list[3] +'\n')
            f.close()
        #output element summary 
            tmp = np.vstack( (np.arange(nelement), tri_areas, depth_0, depth_q, depth_0-depth_q, \
                h_new, h_new-depth_q, abs(h_new-depth_q)-abs(depth_0-depth_q) ) ).transpose() 
            output_file_name = outfile_prefix+'elm_sum_damp'+ damp_str + '_'+ damp_bnd_str + 'vol'+vol_coeff_str + 'face'+ face_coeff_str + '.csv'
            np.savetxt(output_file_name, tmp, delimiter=',')
            shapefile_name = outfile_prefix+'elm.shp'
            i = 2
            while os.path.exists(shapefile_name):
                shapefile_name = shapefile_name.replace('.shp','') + str(i) + '.shp'
                i = i + 1
            elm_attribute_def = [['h0','float'],['hq','float'],['hnew','float'],['h0_hq','float'], \
                                 ['hnew_hq','float'],['diff_new_0','float']]
            tmp = np.column_stack(( elements,depth_0,depth_q,h_new,depth_0-depth_q,h_new-depth_q, abs(h_new-depth_q)-abs(depth_0-depth_q)))                   
            create_element_shapefile(shapefile_name, tmp, nodes, n_attribute = 6, attribute_def = elm_attribute_def)
        #output .prop for element error
            tmp = np.vstack( (np.arange(nelement)+1, depth_0-depth_q )).transpose()
            output_file_name = outfile_prefix+'elm_err_before_damp'+ damp_str + '_'+ damp_bnd_str + 'vol'+vol_coeff_str + 'face'+ face_coeff_str + '.prop'
            np.savetxt(output_file_name, tmp, delimiter='   ', fmt='%10d  %10f')
            tmp = np.vstack( (np.arange(nelement)+1, h_new-depth_q )).transpose()
            output_file_name = outfile_prefix+'elm_err_after_damp'+ damp_str + '_'+ damp_bnd_str + 'vol'+vol_coeff_str + 'face'+ face_coeff_str + '.prop'
            np.savetxt(output_file_name, tmp, delimiter='   ', fmt='%10d  %10f' )
        #output xyz for element error
            output_file_name = outfile_prefix+'elm_err_after_damp'+ damp_str + '_'+ damp_bnd_str + 'vol'+vol_coeff_str + 'face'+ face_coeff_str + '.xyz'
            tmp = np.hstack( ( np.average(tri,axis=1)[:,0:2] , (h_new-depth_q).reshape(nelement,1) ) )
            np.savetxt(output_file_name, tmp, delimiter=',')
        #output face summary
            tmp = np.vstack( (np.arange(nface), nodes_coor[face_nodes[:,0]][:,0], nodes_coor[face_nodes[:,0]][:,1], \
                nodes_coor[face_nodes[:,1]][:,0], nodes_coor[face_nodes[:,1]][:,1], \
                face_len, face_h0, face_hq, face_h0-face_hq, face_h_new, face_h_new-face_hq,  \
                abs(face_h0-face_hq)-abs(face_h_new-face_hq) ) ).transpose()
            output_file_name = outfile_prefix+'face_sum_damp'+ damp_str +'_'+ damp_bnd_str +  'vol' + vol_coeff_str + 'face'+ face_coeff_str + '.csv'
            np.savetxt(output_file_name, tmp, delimiter=',')
        #output xyz for face error
            output_file_name = outfile_prefix+'face_err_after_damp'+ damp_str + '_'+ damp_bnd_str + 'vol'+vol_coeff_str + 'face'+ face_coeff_str + '.xyz'
            tmp = np.vstack( ( (nodes_coor[face_nodes[:,0]][:,0] + nodes_coor[face_nodes[:,1]][:,0])*0.5, \
                               (nodes_coor[face_nodes[:,0]][:,1] + nodes_coor[face_nodes[:,1]][:,1])*0.5, \
                                face_h_new-face_hq ) ).transpose()
            np.savetxt(output_file_name, tmp, delimiter=',')
        #output face shapefile    
            shapefile_name = outfile_prefix+'face.shp'
            i = 2
            while os.path.exists(shapefile_name):
                shapefile_name = shapefile_name.replace('.shp','') + str(i) + '.shp'
                i = i + 1
            face_attribute_def = [['h0','float'],['hq','float'],['hnew','float'],['h0_hq','float'],['hnew_hq','float'],['diff_new_0','float']]
            tmp = np.column_stack(( face_nodes, face_h0, face_hq, face_h_new, face_h0-face_hq,face_h_new-face_hq, \
                                    abs(face_h0-face_hq)-abs(face_h_new-face_hq)))                   
            create_face_shapefile(shapefile_name, tmp, nodes, n_attribute = 6, attribute_def = face_attribute_def)
        #output node summary
            tmp = np.vstack((np.arange(nnode), nodes_coor[:,0], nodes_coor[:,1], \
                nodes_coor[:,2], nodes_coor[:,3], h0_node, node_h_new, nodes_coor[:,3]-node_h_new ) ).transpose()
            output_file_name = outfile_prefix+'node_sum_damp'+ damp_str +'_'+ damp_bnd_str +  'vol' + vol_coeff_str + 'face'+ face_coeff_str + '.csv'
            np.savetxt(output_file_name, tmp, delimiter=',')
        #output xyz filr for nodes
            tmp = np.vstack( ( nodes_coor[:,0], nodes_coor[:,1], nodes_coor[:,3]-node_h_new ) ).transpose()
            output_file_name = outfile_prefix+'node_sum_damp'+ damp_str +'_'+ damp_bnd_str +  'vol' + vol_coeff_str + 'face'+ face_coeff_str + '.xyz'
            np.savetxt(output_file_name, tmp, delimiter=',')
        #output node shapefile    
            shapefile_name = outfile_prefix+'node.shp'
            i = 2
            while os.path.exists(shapefile_name):
                shapefile_name = shapefile_name.replace('.shp','') + str(i) + '.shp'
                i = i + 1
            node_attribute_def = [['z0','float'],['zref','float'],['h0','float'],['hnew','float'],['znew','float']]
            tmp = np.column_stack(( nodes_coor, h0_node, node_h_new, nodes_coor[:,3]-node_h_new ))                   
            create_point_shapefile(shapefile_name, tmp, n_attribute = 5, attribute_def = node_attribute_def)    
        #output node summary along the bank
            if nbnd > 0:
                tmp = np.array([(x[0]+1, x[3], x[4], nodes_coor[int(x[1])][2], nodes_coor[int(x[1])][3], \
                    nodes_coor[int(x[1])][3]-node_h_new[int(x[1])]) for x in land_boundaries])
                output_file_name = outfile_prefix+'node_bnd_damp'+ damp_str +'_'+ damp_bnd_str +  'vol' + vol_coeff_str + 'face'+ face_coeff_str + '.csv'
                np.savetxt(output_file_name, tmp, delimiter=',')
            del tmp 
    return optimized_z
        
def opt_parm_check(opt_parm_list):
    '''
    check whether specified optimization parameters are valid
    '''
    error = 0
    for parm in opt_parm_list:
        if len(parm) <> 4:
            print 'Error: incorrect number of optimization parameters.'
            error = 1
            print parm
        elif (min(parm) < 0) or (parm[2]==0 and parm[3] ==0):
            print 'Error: invalid optimization parameters.'
            print parm
            error = 1
    if (error == 1) :
        return -1
    else:   
        return 0
    
def read_opt_parm(opt_parm_file):  
    f = open(opt_parm_file, 'r')
    all_lines = f.readlines()
    f.close()      
    opt_parm_list = [ line.strip().split(',') for line in all_lines ]
    del all_lines
    opt_parm_list = [ map(float, x) for x in opt_parm_list ] 
    flag = opt_parm_check(opt_parm_list)    
    if flag == 0:
        return opt_parm_list
    else:
        return []      
    
def read_gr3(infile_gr3):
    f = open(infile_gr3,'r')
    all_lines = f.readlines()
    f.close()
    nelm = int(all_lines[1].strip().split()[0])
    nnode = int(all_lines[1].strip().split()[1])
    #read node information
    nodelines = [line.strip().split()[1:4] for line in all_lines[2:nnode+2]]
    nodes = np.array([(line[0:]) for line in nodelines],dtype='float64')
    del nodelines
    #read element information
    elementlines = [line.strip().split()[0:5] for line in all_lines[nnode+2:nnode+2+nelm]]
    elements = np.array([(line[2:]) for line in elementlines], dtype=np.int)    
    elements = elements - 1
    del elementlines 
    #read land boubdary information
    n_open_bnd_segments = int(all_lines[nnode+2+nelm].strip().split()[0])
    n_open_bnd_nodes = int(all_lines[nnode+2+nelm+1].strip().split()[0])
    n_lines_to_skip = nnode + 2 + nelm + n_open_bnd_nodes + n_open_bnd_segments + 2
    n_land_bnd_segments = int(all_lines[n_lines_to_skip].strip().split()[0])
    bnd_list = []
    if n_land_bnd_segments > 0 :
        #n_land_bnd_nodes = int(all_lines[n_lines_to_skip+1].strip().split()[0])
        n_bnd_node_pairs = 0
        idx_start = n_lines_to_skip + 2
        for i in range(n_land_bnd_segments):
            n_bnd_nodes = int(all_lines[idx_start].strip().split()[0])
            bnd_nodes = [line.strip().split()[0] for line in all_lines[idx_start+1:idx_start+1+n_bnd_nodes]]
            bnd_list = bnd_list + ([[n_bnd_node_pairs+j,int(bnd_nodes[j])-1,int(bnd_nodes[j+1])-1,i,j] for j in range(n_bnd_nodes-1)])
            idx_start = idx_start + n_bnd_nodes + 1
            n_bnd_node_pairs = n_bnd_node_pairs + n_bnd_nodes -1
        del bnd_nodes           
    land_boundaries = np.array(bnd_list, dtype=np.int)
    del bnd_list
    return (elements, nodes, land_boundaries) 
    
def read_2dm(infile_2dm):
    f = open(infile_2dm,'r')
    all_lines = f.readlines()
    f.close()
    #read element information
    elementlines = [line.strip().split()[1:5] for line in all_lines if line.startswith("E3T")]
    elements = np.array([line[1:] for line in elementlines], dtype=np.int)    
    elements = elements - 1
    del elementlines
    #read node information
    nodelines = [line.strip().split()[1:5] for line in all_lines if line.startswith("ND")]
    nodes = np.array([(line[1:]) for line in nodelines],dtype='float64')
    #to be consistent with the gr3 format 
    nodes = nodes*np.array([1,1,-1]) 
    del nodelines
    #read land boubdary information
    NSlines = [line.strip().split()[1:11] for line in all_lines if line.startswith("NS")]
    #create an empty list to hole the informaton for land boundaries    
    bnd_list = []
    if len(NSlines) > 0:
        group_id = 0
        group_idx = 0
        bnd_nodes = []      
        for line in NSlines:
            for i in range(len(line)):
                if int(line[i]) < 0:
                    bnd_nodes.append((group_id, group_idx, int(line[i].lstrip('-'))))
                    group_id = group_id +1
                    group_idx = 1
                    break
                else: 
                    bnd_nodes.append((group_id, group_idx, int(line[i])))
                group_idx = group_idx + 1
        bnd_idx1 = 0
        bnd_idx2 = 0
        bnd_list.append((bnd_idx1,int(bnd_nodes[0][2])-1,int(bnd_nodes[1][2])-1,bnd_nodes[0][0],bnd_idx2))  
        for i in range(1,len(bnd_nodes)-1):
            if bnd_nodes[i][0] == bnd_nodes[i-1][0]:
                bnd_idx2 = bnd_idx2 + 1
                bnd_idx1 = bnd_idx1 + 1
            else:
                bnd_list.pop()
                bnd_idx2 = 0
            bnd_list.append((bnd_idx1,int(bnd_nodes[i][2])-1,int(bnd_nodes[i+1][2])-1,bnd_nodes[i][0],bnd_idx2)) 
        del bnd_nodes    
    land_boundaries = np.array(bnd_list)
    del NSlines
    del bnd_list
    return (elements, nodes, land_boundaries) 
    
def read_boundary_node_csv(bnd_csv_file):
    '''Read boundary node list file 
       Reads as (no of adjacent node pairs, 1st node no, 2nd node no, boundary segment, sequence in segment) usually generated by other scripts, can be used with .2dm or .gr3 file without boundary node information. Mainly used during the developmnet phase, 
       kept for debugging purpose
    '''    
    f = open(bnd_csv_file,'r')
    all_lines = f.readlines()
    f.close()
    #nbnd = len(all_lines)
    bnd_list = [ line.split(',')[0:5] for line in all_lines ]
    bnd_list = [ map(int, x) for x in bnd_list ]
    bnd_list = [ (x[0],x[1],x[2],x[3],x[4]) for x in bnd_list ]
    return np.array(bnd_list, dtype=np.int)

def write_gr3(elements, nodes, outfile):
    '''
    create a .gr3 file with optimized elevations, but no boundary information 
    '''
    nelm = len(elements)
    nnode = len(nodes)
    f = open(outfile,"w")
    f.write("hgrid.gr3\n")
    f.write("%s %s    ! number of elements, nodes\n" % (nelm, nnode))
    padding = 2
    import math
    maxnum = int(math.log10(max(nelm,nnode))) + padding
    ifmt = "%" + ("%ii" % maxnum)
    ifmtj = "%-" + ("%ii" % maxnum)
    ffmt = "%18.8f"
    nfmt  = ifmtj + ffmt*3 + "\n" 
    for i in range(nnode):     
        f.write(nfmt % (i+1, nodes[i,0], nodes[i,1], -1*nodes[i,2]))
        
    shapecode = 3
    efmt  = ifmtj + "%2i"+ifmt*3 + "\n" 
    for i in range(nelm):
        f.write(efmt % (i+1, shapecode, elements[i,0]+1,elements[i,1]+1,elements[i,2]+1))
    
    f.write('0 !Number of open boundaries\n')
    f.write('0 !Number of open boundary nodes\n')
    f.write('0 !Number of land boundaries\n')
    f.write('0 !Number of land boundary nodes (including islands)\n')
    f.close()
    return
   
def create_driver():
    spatialReference = osgeo.osr.SpatialReference() 
    spatialReference.ImportFromProj4('+proj=utm +zone=10N +ellps=NAD83 +datum=NAD83 +units=m')
    driver = ogr.GetDriverByName('ESRI Shapefile')
    return (driver, spatialReference)
     
def create_point_shapefile(file_path, point_data, n_attribute = 0, attribute_def = None ):
    ''' Create point shapefile based on xy coordinates specifed in point_data,
    FID is zero-based, i.e. starting at 0

    Parameters
    ----------
    
    file_path  : str
        location and name of the shapefile to be created, string
    
    attribute_def  : str
        name and type of point attributes, list [[name, datatype]]
        valid datatype: float (=> OFTReal), int (=> OFTInteger), str (=>OFTString)
    
    point_data : array
        x,y and attribute values for each points, numpy array                

    '''
    if os.path.exists(file_path):
        print 'Error: cannot create point shapefile -- file already exists.'
        return -1
    if not file_path.endswith('.shp'):
        file_path = file_path + '.shp'   
    if point_data.shape[1] < 2:
        print 'Error: cannot create point shapefile -- insufficient point informaton.'
        return -1  
    if n_attribute > 0:
        if attribute_def == None or point_data.shape[1] - 2 < n_attribute or \
            ( attribute_def != None and len(attribute_def) < n_attribute ):
            print 'Error: cannot create point shapefile -- insufficient attribute informaton.'
            return -1  
        
    (driver, spatialReference) = create_driver()
    pointFile = driver.CreateDataSource(file_path)
    pointLayer = pointFile.CreateLayer('points',spatialReference, osgeo.ogr.wkbPoint)
    #create new field definitions for attributes
    if n_attribute > 0:
        for attr in attribute_def:
            if attr[1] == 'int':
                new_field = ogr.FieldDefn(attr[0], ogr.OFTInteger)
            elif attr[1] == 'float':
                new_field = ogr.FieldDefn(attr[0], ogr.OFTReal)
            elif attr[1] == 'str':
                new_field = ogr.FieldDefn(attr[0], ogr.OFTString)
            else:
                print 'Error: cannot create point shapefile -- incorrect attribute data type.'
                return -1 
            pointLayer.CreateField(new_field)
        
    layerDefn = pointLayer.GetLayerDefn()
    point = osgeo.ogr.Geometry(osgeo.ogr.wkbPoint)
    feature = osgeo.ogr.Feature(layerDefn) 
    for i in range(len(point_data)):
        point.AddPoint(float(point_data[i][0]),float(point_data[i][1]))  
        feature.SetGeometry(point) 
        #assign attribute values
        if n_attribute > 0:
            for j in range(len(attribute_def)):
                feature.SetField(j,point_data[i][2+j])
        a = pointLayer.CreateFeature(feature)
    pointFile.Destroy()
    return 0
  
def create_node_shapefile(file_path, nodes): 
    ''' Create shapefile for nodes
    
    Parameters
    ----------
    file_path  : str
        Location of the shapefile to be created, string
    
    nodes  : array
        x,y, depth for each points, numpy array  
        
    '''
    #convert depth to elevation
    nodes = nodes*np.array([1,1,-1]) 
    return create_point_shapefile(file_path, nodes, 1, [['z','float']])
  
def create_face_shapefile(file_path, face_data, node_data, n_attribute = 0, attribute_def = None):
    '''Create a shapefile for edges 
    
    Shapefile is based on based on node # specified in face_data and xy coordinates specifed in node_data, FID is zero-based, i.e. starting at 0

    Parameters
    ----------

    file_path  : str
        location and name of the shapefile to be created, string
    
    face_data  : array
        two node numbers for each face, numpy array
    
    node_data  : array
        x,y and attribute values for each points, numpy array 
    
    n_attribute : int
        number of attribute defined for the face feature
    
    attribute_def : List
        name and type of point attributes, list [[name, datatype]]
             valid datatype: float (=> OFTReal), int (=> OFTInteger), str (=>OFTString)               
    '''
    if os.path.exists(file_path):
        print 'Error: cannot create face shapefile -- file already exists.'
        return -1
    if not file_path.endswith('.shp'):
        file_path = file_path + '.shp'  
    if face_data.shape[1] < 2:
        print 'Error: cannot create face shapefile -- insufficient face node information.'
        return -1           
    if node_data.shape[1] < 2:
        print 'Error: cannot create face shapefile -- insufficient node xy-coordinate informaton.'
        return -1  
    if n_attribute > 0:
        if attribute_def == None or face_data.shape[1] - 2 < n_attribute or \
            ( attribute_def != None and len(attribute_def) < n_attribute ):
            print 'Error: cannot create face shapefile -- insufficient attribute informaton.'
            return -1 
    (driver, spatialReference) = create_driver()
    faceFile = driver.CreateDataSource(file_path)
    faceLayer = faceFile.CreateLayer('faces', spatialReference, ogr.wkbLineString)
    #create new field definitions for attributes
    for i in range(2):
        new_field = ogr.FieldDefn('node' + str(i), ogr.OFTInteger)
        faceLayer.CreateField(new_field)
    if n_attribute > 0:
        for attr in attribute_def:
            if attr[1] == 'int':
                new_field = ogr.FieldDefn(attr[0], ogr.OFTInteger)
            elif attr[1] == 'float':
                new_field = ogr.FieldDefn(attr[0], ogr.OFTReal)
            elif attr[1] == 'str':
                new_field = ogr.FieldDefn(attr[0], ogr.OFTString)
            else:
                print 'Error: cannot create element shapefile -- incorrect attribute data type.'
                return -1 
            faceLayer.CreateField(new_field)
        
    layerDefn = faceLayer.GetLayerDefn()
    feature = osgeo.ogr.Feature(layerDefn)
    myLine = ogr.Geometry(ogr.wkbLineString) 
    for i in range(len(face_data)):
        for j in range(2):
            node = int(face_data[i][j])
            myLine.AddPoint(node_data[node][0],node_data[node][1])
        a = feature.SetGeometry(myLine)
        #assign attribute values
        for j in range(2):
            feature.SetField(j,int(face_data[i][j]))
        if n_attribute > 0:
            for j in range(len(attribute_def)):
                feature.SetField(j+2,face_data[i][2+j])
        a = faceLayer.CreateFeature(feature)
        myLine.Empty()   
    faceFile.Destroy()
    return 0
  
def create_element_shapefile(file_path, element_data, node_data, n_node_per_elm = 3, n_attribute = 0, attribute_def = None):
    ''' Create a shapefile representing elements as polygons
    Each element one is defined  by a series of nodes (zero-based in element_data), 
    where xy coordinates are provided in node_data.

    Parameters
    ----------
    
    file_path  : str
        Location and name of the shapefile to be created, string
    
    element_data : array
        Node number for each element
    
    node_data : array
        x,y and attribute values for each node, numpy array 
    
    n_node_per_elm : int
        number of nodes (vertices) per element, deafult = 3
    
    n_attribute  : int
        number of attributes defined for individual elements, default = 0
    
    attribute_def : List
        name and type of point attributes, list [[name, datatype]]
    
    valid datatype: float (=> OFTReal), int (=> OFTInteger), str (=>OFTString)
        prototype datatype
    
    '''
    if os.path.exists(file_path):
        print 'Error: cannot create element shapefile -- file already exists.'
        return -1
    if not file_path.endswith('.shp'):
        file_path = file_path + '.shp'
    if element_data.shape[1] < n_node_per_elm:
        print 'Error: cannot create element shapefile -- insufficient element node information.'
        return -1    
    if node_data.shape[1] < 2:
        print 'Error: cannot create element shapefile -- insufficient node xy-coordinate informaton.'
        return -1
    if n_attribute > 0:
        if attribute_def == None or element_data.shape[1]-3 < n_attribute or \
                ( attribute_def <> None and len(attribute_def) < n_attribute ):
            print 'Error: cannot create element shapefile -- insufficient attribute information.'
            return -1
   
    (driver, spatialReference) = create_driver()
    elmFile = driver.CreateDataSource(file_path)
    elmLayer = elmFile.CreateLayer('elements', spatialReference, ogr.wkbPolygon)
    #create new field definitions for attributes
    for i in range(n_node_per_elm):
        new_field = ogr.FieldDefn('node' + str(i), ogr.OFTInteger)
        elmLayer.CreateField(new_field)
    if n_attribute > 0:
        for attr in attribute_def:
            if attr[1] == 'int':
                new_field = ogr.FieldDefn(attr[0], ogr.OFTInteger)
            elif attr[1] == 'float':
                new_field = ogr.FieldDefn(attr[0], ogr.OFTReal)
            elif attr[1] == 'str':
                new_field = ogr.FieldDefn(attr[0], ogr.OFTString)
            else:
                print 'Error: cannot create element shapefile -- incorrect attribute data type.'
                return -1 
            elmLayer.CreateField(new_field)
        
    layerDefn = elmLayer.GetLayerDefn()
    feature = osgeo.ogr.Feature(layerDefn)
    myRing = ogr.Geometry(type=ogr.wkbLinearRing)
    myPoly = ogr.Geometry(type=ogr.wkbPolygon)  
    for i in range(len(element_data)):
        for j in range(n_node_per_elm):
            node = int(element_data[i][j])
            myRing.AddPoint(node_data[node][0],node_data[node][1])
        #add the 1st node to close the polygon    
        myRing.AddPoint(node_data[int(element_data[i][0])][0],node_data[int(element_data[i][0])][1])
        a = myPoly.AddGeometry(myRing)
        a = feature.SetGeometry(myPoly)
        #assign attribute values
        for j in range(n_node_per_elm):
            feature.SetField(j,int(element_data[i][j]))
        if n_attribute > 0:
            for j in range(len(attribute_def)):
                feature.SetField(j+n_node_per_elm,element_data[i][3+j])
        a = elmLayer.CreateFeature(feature)
        myRing.Empty()
        myPoly.Empty()    
    elmFile.Destroy()
    return 0
  
def point_inside_polygon(x, y, poly):
    ''' Decide if a point is inside a polygon using the ray casting method

    Parameters
    ----------
    
    x,y : float
        Point to consider
    
    poly: List
        list of pairs (x,y) containing the coordinates of the polygon's vertices. 
    
    Returns
    -------
    inside : boolean
        True if inside
    
    '''
    n = len(poly)
    inside = False
    p1x, p1y = poly[0]
    for i in range(n):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y-p1y) * (p2x-p1x) / (p2y-p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y
    return inside
 
def create_partial_domain(elements, nodes, land_boundaries, polygon):
    ''' Select nodes, elements, and land_boundaries based on specified polygon
    '''
    nodes_mask = np.array([i if point_inside_polygon(nodes[i][0],nodes[i][1],polygon) else -1 for i in range(len(nodes))])
    sel_nodes = np.column_stack((nodes,nodes_mask))
    sel_nodes = sel_nodes[sel_nodes[:,3] > -1]
    sel_nodes = sel_nodes[:,0:3]
    sel_node_list = nodes_mask[nodes_mask > -1]
    #only keep the elements where all three nodes are in the selected node list
    sel_elms = np.array([x for x in elements if x[0] in sel_node_list and x[1] in sel_node_list and x[2] in sel_node_list])
    nelm = len(sel_elms)
    sel_elms = sel_elms.reshape(3*nelm)
    sel_elms = np.array([np.nonzero(sel_node_list == x)[0][0] for x in sel_elms])
    sel_elms = sel_elms.reshape(nelm,3)
    #only keep the land_boundaies where all two nodes are in the selected node list
    sel_land_bnds = np.array([x for x in land_boundaries if x[1] in sel_node_list and x[2] in sel_node_list])
    sel_land_bnds[:,1] = np.array([np.nonzero(sel_node_list==x)[0][0] for x in sel_land_bnds[:,1]])
    sel_land_bnds[:,2] = np.array([np.nonzero(sel_node_list==x)[0][0] for x in sel_land_bnds[:,2]])
    del nodes_mask
    del sel_node_list
    return (sel_elms, sel_nodes, sel_land_bnds)
    
def create_arg_parser():
    import argparse
    parser = argparse.ArgumentParser(description='Perform grid optimization with a *.2dm SMS mesh or gr3 file. An optimized gr3 file with extension _opt.gr3 will be created if only one set of optimization parameter specified.')
    parser.add_argument('filename',default = None, help = 'name of 2dm or gr3 file')
    parser.add_argument('demfile', default = None, help = 'file containing list of DEMs. These can be in any form that gdal accepts, which includes ESRI ascii format and GeoTiffs')
    parser.add_argument('optparm', default = None, help = 'file containing optimization parameters: damp, damp_shoreline, face_coeff, volume_coeff')
    parser.add_argument('--optfile', default = '', help = 'name for the gr3 file for the optimized results')
    parser.add_argument('--detailed_outputs',action = 'store_true', default=False, help='whether to write detailed optimization results to files, default=False')  
    parser.add_argument('--prefix', default = '', help = 'prefix used for output files when --detailed_outputs flag is set')
    parser.add_argument('--boundary_list', default = '', help = 'boundary node list information (in csv format) when missing from .2dm or .gr3 file, generated from Kijin''s script')
    parser.add_argument('--remove_boundary_faces',action = 'store_true', default=False, help='whether to remove faces along the land boundaries, default=False')
    parser.add_argument('--polygon_file', default = 'None', help = 'polygon definition file used to specify the area for partial grid optimization')
    parser.add_argument('--solver', default = 'L-BFGS-B', help = 'solver used for optimization, either L-BFGS-B (default) or lsqr')
    return parser
    
def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    infile = args.filename
    demlistfile = args.demfile
    optparmfile = args.optparm
    optfile = args.optfile
    output_details = args.detailed_outputs
    prefix_outfile = args.prefix
    bnd_csv_file = args.boundary_list
    remove_bnd_faces = args.remove_boundary_faces
    polygon_file = args.polygon_file
    solver = args.solver
    
    optparmlist = read_opt_parm(optparmfile)
    if len(optparmlist) == 0:
        print 'No valid optimization paramaters specified, default values (0.2, 2, 1, 1) will be used'
        optparmlist = [[0.2, 2, 1, 1]]
    if infile.endswith(".2dm"):
        (elements, nodes, land_boundaries) = read_2dm(infile)
    elif infile.endswith(".gr3"):
        (elements, nodes, land_boundaries) = read_gr3(infile)
    else:
        raise ValueError("Input file format not recognized (no gr3 or 2dm extension)")
    
    if len(land_boundaries) < 1 and len(bnd_csv_file) > 5:
        land_boundaries = read_boundary_node_csv(bnd_csv_file)
    
    if  polygon_file <> 'None':
        if len(polygon_file) == 0 or not os.path.exists(polygon_file): 
            exit('Error: missing polygon definition while parial_grid flag is set. Abort...............')
        else:   
            f = open(polygon_file,'r')
            all_lines = f.readlines()
            polygon = [line.strip().split(',') for line in all_lines]
            polygon = [(float(x[0]),float(x[1])) for x in polygon]
            del all_lines
            print 'Extracting partial grid information.....'
            (elements, nodes, land_boundaries) = create_partial_domain(elements, nodes, land_boundaries, polygon)
    
    if solver != 'L-BFGS-B' and solver != 'lsqr':
        solver = 'L-BFGS-B'
        print 'invalid solver specified, L-BFGS-B will be used'
    
    print 'Running optimization routine.......'    
    z_opt = grid_opt(nodes,elements,land_boundaries,demlistfile,optparmlist,output_details,prefix_outfile,remove_bnd_faces,solver)
    
    if len(z_opt) == 1:
        nodes[:,2] = z_opt[0]
        
        if len(optfile) < 1:
            optfile = infile[:-4] + '_opt.gr3'  
        elif not optfile.endswith('.gr3'):
            optfile = optfile + '.gr3'

        print 'Writing optimized gr3 file ' + optfile + ' .......'
        write_gr3(elements, nodes, optfile)
    print 'Grid optimization completed.'    
            
    
    
if __name__ == "__main__":
    main()
