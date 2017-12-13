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

MAXRUN=10
NPROC=3

def add_farm_parameter(filename,hgrid,default,value,nodes):
		attr_array = np.empty(hgrid.mesh.n_nodes())   
		attr_array.fill(default)    
		attr_array[nodes]=value
		hgrid.write_hgrid(filename, attr_array, False)

def get_nodes(mesh,vertices):
	vertices = map(float, vertices.split())
	poly = Polygon(vertices, '', 0,'none')
	box=poly.box()
	nodes = mesh.find_nodes_in_box(box)
	element_i=[]
	for node_i in nodes:
	    node =mesh.nodes[node_i]
	    flag = poly.check_point_inside_polygon(node[:2])
	    if flag:
	          element_i.append(node_i)

	return element_i
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
		print 'wait for %s to be created' % tidal_cycle
		# wait that the next files get created
		while not os.path.exists(tidal_cycle): 
			time.sleep(1)

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
		nodes=get_nodes(pw.hgrid.mesh,run_parameters['farms'][farm]['vertices'])
		pw.add_farm(farm,run_parameters['farms'][farm]['vertices'],nodes)

		for filename in run_parameters['farms'][farm]['params']:
			default=run_parameters['farms'][farm]['params'][filename]['default']
			value=run_parameters['farms'][farm]['params'][filename]['value']
			add_farm_parameter(os.path.join(run_parameters['run directory'],filename),\
				pw.hgrid,default,value,nodes)



def run(run_parameters):

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


	# NRUN=0
	# while NRUN<MAXRUN:

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
	os.system('rm %s' % (os.path.join(run_parameters['run directory'],'*.yaml')))
	for farm in run_parameters['farms']:
		for filename in run_parameters['farms'][farm]['params']:
			os.system('rm %s' % (os.path.join(run_parameters['run directory'],filename)))


		# ## wait that new yaml file has arived
		# while not os.path.exists(os.path.join(run_parameters['run directory'],'*.yaml')): 
		# 	time.sleep(1)


		# NRUN+=1



if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Do a run')
	parser.add_argument('--yaml', metavar='yaml', type=str,
	                   help='Yaml file')

	args = parser.parse_args()

	## get all the option from the yaml file
	with open(args.yaml ,'r') as f:
		run_parameters = yaml.load(f)

	run(run_parameters)