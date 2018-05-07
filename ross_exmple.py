#!/usr/bin/env python2.7
import sys,os
import numpy as np
import yaml
import subprocess

root='/home/remy/Buisness/0336_SSC_tides/test/' # My root directory


# This is where all functions are located
#root_script=os.path.dirname(os.path.realpath(__file__))
root_script=os.path.join(os.path.dirname(os.path.realpath(__file__)),'SSC')
print root_script+'/func/SSC.py'
# this is the basic parameter can add more to it
P=dict()
P['nproc']=3
params=dict()
params['itur']=0
params['X']=0.01
params['bfric']= 0
params['itr_met']= 1
params['bdrc61']= 1
P['params']=params




for farmid,farmdrag in enumerate(np.arange(.1,.3,.05)):
	farmparams=dict()
	farmparams['farms']=dict()

	farmname='farm_%i' % farmid

	# The farm specific paraemter
	farmparams['farms'][farmname]=dict()
	farmparams['farms'][farmname]['vertices']=np.array([[1730477,5400792],[1738253,5411596],[1749124,5407414],[1744441,5392080],[1734909,5391906]]).tolist()
	farmparams['farms'][farmname]['drag.gr3']=dict()
	farmparams['farms'][farmname]['drag.gr3']['default']=0.0025 # drag coefficient everywhere in the grid
	farmparams['farms'][farmname]['drag.gr3']['value']=float(farmdrag) # drag coefficient in the farm

	# Take care of the input an  output folder
	if not os.path.exists(os.path.join(root,farmname)):
		os.makedirs(os.path.join(root,farmname))
		outputFolder=os.path.join(root,farmname,'out')
		inputFolder=os.path.join(root,farmname,'in')
		os.makedirs(outputFolder)
		os.makedirs(inputFolder)
		os.system('chmod 777 '+inputFolder)
		os.system('chmod 777 '+outputFolder)
	else:
		outputFolder=os.path.join(root,farmname,'out')
		inputFolder=os.path.join(root,farmname,'in')


	# Save the paramter to a file 
	with open(os.path.join(inputFolder,farmname+'.yml'), 'w') as yaml_file: # this would write the yaml file that my function read probably best so we can track
		yaml.dump(P, yaml_file, default_flow_style=False)
		yaml.dump(farmparams, yaml_file, default_flow_style=False)

	# This would run the docker with folder mounted at the right place
	subprocess.call(["docker","run","-v",root_script+":/home/user/SSC/","-v",inputFolder+":/home/user/run/input","-v",outputFolder+":/home/user/run/output/"\
		,"metocean/ssc_tide","python","/home/user/SSC/func/SSC.py","--yaml",os.path.join('/home/user/run/input',farmname+'.yml')])


	# output for you to read should be in outputFolder


