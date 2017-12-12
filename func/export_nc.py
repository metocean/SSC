#!/usr/bin/env python2.7
import sys,os
import netCDF4 as netCDF4

def export_nc(filename,pw,typ='power'):
	nc=netCDF4.Dataset(filename,'r+')
	for farm in pw.farms.keys():
		new_var = nc.createVariable(farm+'_'+typ, 'f8', ('time','nSCHISM_hgrid_node'))
		import pdb;pdb.set_trace()
		new_var[pw.farms[farm]['nodes']]=pw.farms[farm][typ]

	nc.close()