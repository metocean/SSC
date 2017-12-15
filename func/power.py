#!/usr/bin/env python2.7

import os,sys
sys.path.append(os.path.join(os.path.dirname(__file__),'..','Bay_delta_scripts')) 
import subprocess
import numpy as np
import netCDF4
from schism_setup import load_gr3
from schismIO import schismIO
import netCDF4


rho=1025 # kg/m3
def cal_tri_area(a):
    return np.absolute((a[0]*(a[3]-a[5])+a[2]*(a[5]-a[1])+a[4]*(a[1]-a[3]))/2.0)


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
		self.dir= SC.dir
		self.hgrid = load_gr3(os.path.join(SC.dir,'hgrid.gr3'))
		self.Xss=0.02
		areas=get_areas(self.hgrid.mesh,range(0,self.hgrid.mesh.n_elems()))
		self.farms={}
		self.farms['whole']={'nodes':range(0,len(self.hgrid.mesh.nodes[:,0])),\
				'limits':[[self.hgrid.mesh.nodes[:,0].min(),self.hgrid.mesh.nodes[:,0].max()],[self.hgrid.mesh.nodes[:,1].min(),self.hgrid.mesh.nodes[:,1].max()]],\
				'elements':  self.hgrid.mesh.elems,\
				'areas': areas}

	def add_farm(self,name,vertices,nodes,elements,areas):
		vertices=[float(x) for x in vertices.split(' ')]
		self.farms[name]={}
		self.farms[name]['limits']=[[vertices[0::2]],[vertices[1::2]]]
		self.farms[name]['nodes']=nodes
		self.farms[name]['elements']=elements
		self.farms[name]['areas']=areas
		
	

	def check_ss():
		pass

	def get_power(self,file_number):
		nc=netCDF4.Dataset(os.path.join(self.sc.dir,'outputs','schout_%i.nc' % file_number))
		ts=len(nc.variables['time'])
		U=nc.variables['dahv'][:,:,0] #U veloicity [time,nodes]
		V=nc.variables['dahv'][:,:,1] #V velocity [time,nodes]
		Cd=nc.variabes['bottom_drag_coef'][:,:] # bottom drag [time,node]
		
		for farm in self.farms.keys():
			A=self.farms[farm]['areas'] # area in m2 for each element of the farm
			e=self.farms[name]['elements'] # element inside the farm
			n=self.farms[name]['nodes'] # nodes inside the farm
			tri=self.hgrid.mesh.elems[e,:]
			import pdb;pdb.set_trace()
			Cde = np.average(depth_0, axis=1) # bottom drag at each element


		spd=np.sqrt(U**2+V**2)
		P_ts=0.5*rho*spd**3
		P=sum(P_ts)/spd.shape[0]

		for farm in self.farms.keys():
			self.farms[farm]['mean power']=np.mean(P[self.farms[farm]['nodes']])
			self.farms[farm]['power']=P[self.farms[farm]['nodes']]
			self.farms[farm]['power ts']=P_ts[:,self.farms[farm]['nodes']]

	def export_nc(self,file_number,typ='power',outdir=None):
		if outdir is None:
			outdir=self.sc.dir

		old_filename=os.path.join(self.sc.dir,'outputs','schout_%i.nc' % file_number)
		new_filename=os.path.join(outdir,'schout_%i.nc' % file_number)
		os.system('cp '+old_filename+' '+new_filename)
		nc=netCDF4.Dataset(new_filename,'r+')
		for farm in self.farms.keys():
			new_var = nc.createVariable(farm+'_tidal_cycle_'+typ, 'f8', ('nSCHISM_hgrid_node'))
			new_var[self.farms[farm]['nodes']]=self.farms[farm][typ]
			new_var = nc.createVariable(farm+'_timeSeries'+typ, 'f8', ('time','nSCHISM_hgrid_node'))
			new_var[:,self.farms[farm]['nodes']]=self.farms[farm][typ+' ts']

		nc.close()

if __name__ == "__main__":
	sc=schismIO('/home/remy/Buisness/0336_SSC_tides/make_hotstart')
	#sc.create_nectdf_file('dahv.62',to_ts=20*3600*24)
	X=0.02
	pw=power(sc)
	dt=12.42*3600
	N=40 # 5 tidal cycle
	ts=15*60.
	P1=0.
	P2=1000.
	n=0
	while np.abs(P2-P1)/P2 > X :
		P1=P2
		pw.get_power(ts+dt*n,ts+dt*(n+1))
		P2=pw.farms['whole grid']['power']
		n+=1


	print 'Steady state reach at timestep: %i , P1=%f, P2=%f' % (ts+dt*n,P1,P2)
