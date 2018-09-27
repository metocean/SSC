#!/usr/bin/env python2.7

import os
import subprocess
import numpy as np
import netCDF4

soft='combine_output11'
hot='combine_hotstart7'

class schismIO:
	def __init__(self,dirIN):
		self.dir= dirIN

	# def get_file(self,ts):
	# 	F=np.ceil(ts/self.file_length)
	# 	return F
	# def get_timestep(self,ts):
	# 	file_end=np.ceil(ts/self.file_length)
	# 	timestep=np.floor((ts-((file_end-1)*self.file_length))/self.output_dt)
	# 	return timestep

	# def get_first_timestep(self):
	# 	f=os.path.join(self.dir,'mirror.out') 			
	# 	for line in open(f).readlines():
	# 		if 'TIME=' in line:
	# 			elapse=float(line.split(';')[1].split('\n')[0].replace('TIME=',''))
	# 			break
	# 	return elapse

	# def get_last_timestep(self):
	# 	f=os.path.join(self.dir,'mirror.out') 			
	# 	for line in reversed(open(f).readlines()):
	# 		if 'TIME=' in line:
	# 			elapse=float(line.split(';')[1].split('\n')[0].replace('TIME=',''))
	# 			break
	# 	return elapse

	def get_startfile(self):
		f=os.path.join(self.dir,'hotstart.nc') 
		nc=netCDF4.Dataset(f)
		ifile=nc.variables['ifile'][0]
		return ifile



	def create_nectdf_file(self,file_num):
		subprocess.call("%s -b %i -e %i" % (soft,file_num,file_num), cwd="%s" % os.path.join(self.dir,'outputs'),shell = True)

	def create_hotstart_file(self,file_num):
		subprocess.call("%s -i %i " % (hot,file_num), cwd="%s" % os.path.join(self.dir,'outputs'),shell = True)

if __name__ == "__main__":
	sc=schismIO('/home/remy/Buisness/0336_SSC_tides/make_hotstart')
	te=sc.get_last_timestep()
	sc.create_nectdf_file('dahv.62',to_ts=te)


