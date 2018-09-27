#!/usr/bin/env python2.7

import os,sys,glob
sys.path.append(os.path.join(os.path.dirname(__file__),'..','Bay_delta_scripts')) 
import subprocess
import numpy as np
import netCDF4
from schism_setup import load_gr3
from schismIO import schismIO
import netCDF4
from collections import Counter

rho=1025 # kg/m3
def cal_tri_area(a):
    return np.absolute((a[0]*(a[3]-a[5])+a[2]*(a[5]-a[1])+a[4]*(a[1]-a[3]))/2.0)

def tri_area(x,y):
	# area from x, y which are (n,3)
	# tested in matlab version compared to polyarea.m same to to 10^-9m for  1500m^2 or larger triangles 
	# get lengths of 3 sides for herons formula
	# ross vennel May 2018
	a=np.sqrt((x[:, 1]-x[:,0])**2+(y[:, 1]-y[:,0])**2)
	b=np.sqrt((x[:, 2]-x[:,1])**2+(y[:, 2]-y[:,1])**2)
	c=np.sqrt((x[:, 0]-x[:,2])**2+(y[:, 0]-y[:,2])**2)

	s=(a+b+c)/2.0
	A=np.sqrt(s*(s-a)*(s-b)*(s-c))
	return  A

def get_areas(mesh,Elems):
	ref=np.zeros((mesh.nodes.shape[0],1))
	nodes_coor = np.hstack((mesh.nodes[:,0:2],ref))
#populate (x,y) and elevation information for each triangular element
	tri = np.zeros((len(Elems),3,3))
	for itri in range(len(Elems)):
		for ivert in range(3):
			tri[itri,ivert,:] = nodes_coor[mesh.elems[Elems[itri]][ivert]] 

	areas = cal_tri_area(tri[:,:,0:2].reshape(len(Elems),6).transpose())

	return areas


class power:
	def __init__(self,SC):
		self.sc=SC
		self.max_cycle=20
		self.dir= SC.dir
		self.bckg_value=0
		self.value=[]
		self.hgrid = load_gr3(os.path.join(SC.dir,'hgrid.gr3'))
		self.Xss=0.02
		areas=get_areas(self.hgrid.mesh,range(0,self.hgrid.mesh.n_elems()))
		self.farms={}
		self.farms['whole']={'nodes':range(0,len(self.hgrid.mesh.nodes[:,0])),\
				'limits':[[self.hgrid.mesh.nodes[:,0].min(),self.hgrid.mesh.nodes[:,0].max()],[self.hgrid.mesh.nodes[:,1].min(),self.hgrid.mesh.nodes[:,1].max()]],\
				'elements':  range(0,self.hgrid.mesh.n_elems()),\
				'areas': areas}

	def add_farm(self,name,vertices,nodes,elements,areas):
		#vertices=[float(x) for x in vertices.split(' ')]
		self.farms[name]={}
		self.farms[name]['limits']=[[float(x[0]) for x in vertices],[float(x[1]) for x in vertices]]
		self.farms[name]['nodes']=nodes
		self.farms[name]['elements']=elements
		self.farms[name]['areas']=areas
		
	

	def check_ss():
		pass


	def get_power(self,file_number):
		nc=netCDF4.Dataset(os.path.join(self.sc.dir,'outputs','schout_%i.nc' % file_number))
		ts=len(nc.variables['time'])
		dt=nc.variables['time'][2]-nc.variables['time'][1]
		U=nc.variables['dahv'][:,:,0] #U veloicity [time,nodes]
		V=nc.variables['dahv'][:,:,1] #V velocity [time,nodes]
		Cd=np.zeros((U.shape[0],U.shape[1]))+0.01#nc.variables['bottom_drag_coef'][:,:] # bottom drag [time,node]
		
		
		for farm in self.farms.keys():
			A=self.farms[farm]['areas'] # area in m2 for each element of the farm
			e=self.farms[farm]['elements'] # element inside the farm
			n=self.farms[farm]['nodes'] # nodes inside the farm
			tri=self.hgrid.mesh.elems[e,:]
			if farm == 'whole':
				bckg_value=0
			else:
				bckg_value=self.bckg_value

			Cde=Cd[0,tri]-bckg_value ## rmove backgroud value
			Ue=U[:,tri]
			Ve=V[:,tri]

			Cde = np.average(Cde, axis=1) # bottom drag at each element shape (86995,)
			Ue=np.average(Ue,axis=2)
			Ve=np.average(Ve,axis=2)
			spd=np.sqrt(Ue**2+Ve**2) # Speed time series at each element (92, 86995)
			
			P_ts=1025*np.sum(A[np.newaxis,:] *Cde*spd**3,axis=1) # Power time series at each element (92, )

			mean_P_ts=np.trapz(P_ts,  dx=1.0, axis=0)/(P_ts.shape[0]-1)

			self.farms[farm]['mean power']=mean_P_ts # average power over one tidal cycle over the whole farm (1,)
			self.farms[farm]['mean power ts']=P_ts# power timesries for the whole farm (92,)

	def export_nc(self,file_number,nTC,typ='power',outdir=None,params=None):
		if outdir is None:
			outdir=self.sc.dir

		old_filename=os.path.join(self.sc.dir,'outputs','schout_%i.nc' % file_number)
		new_filename=os.path.join(outdir,'schout_%i.nc' % file_number)
		grid_filename=os.path.join(outdir,'grid.nc')
		os.system('cp '+old_filename+' '+new_filename)

		#copy the hotstart file if any and add the laster Powr calculation for steady state
		hot_files=glob.glob(os.path.join(self.sc.dir,'outputs','hotstart_it=*.nc'))
		if len(hot_files)>0:
			os.system('cp '+hot_files[0]+' '+os.path.join(outdir,'hotstart.nc'))
			nc=netCDF4.Dataset(os.path.join(outdir,'hotstart.nc'),'r+')
			new_var = nc.createVariable('P2', 'f8', ('one'))
			new_var[0]=self.farms['whole']['mean power']
			nc.close()



		# crate th grid file
		os.system('ncks -O -v SCHISM_hgrid,SCHISM_hgrid_face_nodes,'\
			'SCHISM_hgrid_edge_nodes,SCHISM_hgrid_node_x,SCHISM_hgrid_node_y,'\
			'SCHISM_hgrid_face_x,SCHISM_hgrid_face_y,SCHISM_hgrid_edge_x,SCHISM_hgrid_edge_y,depth'\
			' %s %s' % (new_filename,grid_filename))
                
		print new_filename
		print grid_filename
		nc=netCDF4.Dataset(new_filename,'r+')

		if params is not None:
			for att in params.keys():
				nc.setncattr(att, str(params[att]))



		new_var = nc.createVariable('number_of_TC', 'i4', ('one'))
		new_var[0]=nTC

		self.value=[self.bckg_value]+self.value

		
		for i,farm in enumerate(self.farms.keys()):
			new_var = nc.createVariable(farm+'_mean_'+typ, 'f8', ('one'))
			new_var[0]=self.farms[farm]['mean '+typ]
			new_var.drag=self.value[i]

			new_var = nc.createVariable(farm+'_mean_timeSeries'+typ, 'f8', ('time'),fill_value=0)
			new_var[:]=self.farms[farm]['mean '+typ+' ts']



		nc.close()

if __name__ == "__main__":
	sc=schismIO('/home/remy/Buisness/0336_SSC_tides/make_hotstart2')
	#sc.create_nectdf_file('dahv.62',to_ts=20*3600*24)
	X=0.001
	pw=power(sc)
	dt=12.42*3600
	N=40 # 5 tidal cycle
	ts=15*60.
	P1=0.
	P2=1000.
	n=1
	while np.abs(P2-P1)/P2 > X :
		P1=P2
		pw.get_power(n)
		P2=pw.farms['whole']['mean power']
		print 'power at %i is %.f' % (n,P2)
		n+=1


	print 'Steady state reach at timestep: %i , P1=%f, P2=%f' % (n,P1,P2)
