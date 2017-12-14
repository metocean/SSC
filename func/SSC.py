#!/usr/bin/env python2.7
import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),'..','Bay_delta_scripts'))
sys.path.append(os.path.join(os.path.dirname(__file__),'..','func')) 
import yaml
import argparse
from param import param
import subprocess
import signal
import time
from schismIO import schismIO
from power import power
import numpy as np
from export_nc import export_nc
from polygons import *
import timeit


NPROC=3

def add_farm_parameter(filename,hgrid,default,value,nodes):
		attr_array = np.empty(hgrid.mesh.n_nodes())   
		attr_array.fill(default)    
		attr_array[nodes]=value
		hgrid.write_hgrid(filename, attr_array, False)

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


def get_nodes_elements(mesh,vertices):
	vertices = map(float, vertices.split())
	poly = Polygon(vertices, '', 0,'none')
	box=poly.box()
	nodes = mesh.find_nodes_in_box(box)
	Nodes=[]
	Elems=[]
	for node_i in nodes:
		node =mesh.nodes[node_i]
		flag = poly.check_point_inside_polygon(node[:2])
		if flag:
			Nodes.append(node_i)
			e=mesh.get_elems_i_from_node(node_i)
			Elems+=map( lambda x: int(x), e )

	Elems=list(set(Elems))

	return Nodes,Elems
def run_schism(mode,schism=None,proc=None,dirout=None):
#run the model in the background
	if mode is 'run':
		try:
			subprocess.Popen('/opt/mpich/bin/mpd')
		except:
			pass
		proc=subprocess.Popen('mpirun -np %i bash -c "ulimit -s unlimited && %s" &' % (NPROC,schism),\
						cwd="%s" % dirout,\
						shell=True,\
						preexec_fn=os.setsid)
		return proc
	elif mode is 'kill':
		os.killpg(proc.pid, signal.SIGTERM)


def set_params(options,param_file):
	## default ##
	options['params']['unit']='second'
	options['params']['dt']=162
	options['params']['hotstart']=2
	options['params']['rnday']=59
	options['params']['output dt']= 486 
	options['params']['file length']= 44712 #12.42 hours in seconds
	options['params']['output hotstart']= 0
	options['params']['hotstart dt']= 44712 # in hours 


	p=param()
	p.pass_values(options['params'])
	p.write(fileout=param_file)


def search_steady_state(dirout,pw,sc,X):
	P1=0.
	P2=1000.
	n=12
	## first filename
	tidal_cycle=os.path.join(dirout,'outputs','schout_0000_13.nc')
	
	# main loop while steady state not reach
	while np.abs(P2-P1)/P2 > X:
		print 'waiting for %s to be created' % tidal_cycle
		# wait that the next files get created
		start = timeit.default_timer()
		while not os.path.exists(tidal_cycle): 
			time.sleep(1)
		stop = timeit.default_timer()
		sec=stop-start
		m, s = divmod(sec, 60)
		h, m = divmod(m, 60)

		print 'one tidal cycle finished in %dhrs %02dmin %02dsec' % (h,m,s)
		# compile it
		sc.create_nectdf_file(n)
		# get the power for all the farms in the file
		P1=P2
		pw.get_power(n)

		# check power from the whole grid
		P2=pw.farms['whole']['mean power']
		n+=1

		# new file to wait
		tidal_cycle=os.path.join(dirout,'outputs','schout_0000_%i.nc' % (n+1))



	print 'Steady state reach at tidal cycle number : %i , P1=%f, P2=%f' % (n-1,P1,P2)
	return pw,n

	
def include_farm(run_parameters,pw):
	for farm in run_parameters['farms']:
		nodes,elements=get_nodes_elements(pw.hgrid.mesh,run_parameters['farms'][farm]['vertices'])
		areas=get_area(pw.hgrid.mesh,elements)
		pw.add_farm(farm,run_parameters['farms'][farm]['vertices'],nodes,elements,areas)

		for filename in run_parameters['farms'][farm]['params']:
			default=run_parameters['farms'][farm]['params'][filename]['default']
			value=run_parameters['farms'][farm]['params'][filename]['value']
			add_farm_parameter(os.path.join(run_parameters['run directory'],filename),\
				pw.hgrid,default,value,nodes)



def Wrapper(run_parameters):

	## check path and create it
	if not os.path.exists(run_parameters['run directory']):
		os.system('mkdir %s' %run_parameters['run directory'])


	## copy the inputs
	os.system('cp %s %s' %('/home/user/SSC/initial_files/*',run_parameters['run directory']))


	sc=schismIO(run_parameters['run directory']) # this will combine he file as it run
	pw=power(sc) # this wil get the power after 1 tidal cycle


	if not os.path.exists(run_parameters['saving directory']):
		os.system('mkdir %s' %run_parameters['saving directory'])

	## make an outputs directory for schism
	if not os.path.exists(os.path.join(run_parameters['run directory'],'outputs')):
		os.system('mkdir %s' % os.path.join(run_parameters['run directory'],'outputs'))
	else: # delete all the output
		os.system('rm %s' % os.path.join(run_parameters['run directory'],'outputs/*'))


	## create the parameter file
	set_params(run_parameters,os.path.join(run_parameters['run directory'],'param.in'))
	## create the GR3 with polygons and add the farms inside PW
	include_farm(run_parameters,pw)

	## run SCHISM
	proc=run_schism('run',schism='schism',proc=None,dirout=run_parameters['run directory'])
	print 'schism running in the background'


	## Main loop in search of a steady state
	pw,n=search_steady_state(run_parameters['run directory'],\
		pw,\
		sc,\
		run_parameters['params']['X'])

	## save to saving directory
	pw.export_nc(n-1,outdir=run_parameters['saving directory'])
	print 'schism data exported to %s' % run_parameters['saving directory']

	## kill schism
	run_schism('kill',proc=proc)
	print 'schism is killed'

	## delete file from previous run
	os.system('rm %s' % (os.path.join(run_parameters['run directory'],'param.in')))
	for farm in run_parameters['farms']:
		for filename in run_parameters['farms'][farm]['params']:
			os.system('rm %s' % (os.path.join(run_parameters['run directory'],filename)))


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Do a run')
	parser.add_argument('--yaml', metavar='yaml', type=str,
	                   help='Yaml file')

	args = parser.parse_args()

	## get all the option from the yaml file
	with open(args.yaml ,'r') as f:
		run_parameters = yaml.load(f)

	Wrapper(run_parameters)