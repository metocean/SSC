#!/usr/bin/env python2.7
from mako.template import Template
from mako.runtime import Context
from StringIO import StringIO

tmpl='/home/user/SSC/param.tmpl'

def default_value():
	dflt={}
	dflt['year']=2017
	dflt['month']=1
	dflt['day']=1
	dflt['hour']=0

	dflt['dt']=60
	dflt['rnday']=0
	dflt['dramp']=1
	dflt['nspool']=60
	dflt['ihfskip']=2880

	dflt['iout_sta']=0
	dflt['nspool_sta']=1
	dflt['hotout']=1
	dflt['hotout_write']=86400
	dflt['nstep_wwm']=180
	dflt['icou_elfe_wwm']=0
	dflt['imm']=0
	dflt['itransport']=0

	dflt['h0']=0.01

	dflt['itr_met']=3

	dflt['mdc']=36
	dflt['msc']=15

	dflt['age_tracer']=0
	dflt['gen_tracer']=0
	dflt['sed_tracer']=0

	dflt['if_source']=0
	dflt['nramp_ss']=0
	dflt['dramp_ss']=1

	dflt['age_tracer']=0
	dflt['gen_tracer']=0
	dflt['sed_tracer']=0

	dflt['ic_GEN']=0
	dflt['ic_AGE']=0
	dflt['ic_SED']=0
	dflt['ic_SAL']=1
	dflt['ic_TEM']=1
	dflt['inu_GEN']=0
	dflt['hotstart']=0

	## Turbulence
	dflt['itur']=3
	dflt['dfv0']=1.e-6
	dflt['dfh0']=1.e-6
	dflt['turb_met']='"KE"'
	dflt['turb_stab']='"KC"'

	dflt['xlsc']=0.7

	dflt['inv_atm_bnd']=0
	dflt['inunfl']=0
	## Wind option
	dflt['nws']=0
	dflt['wtiminc']=3600
	dflt['nrampwind']=0
	dflt['drampwind']=1
	dflt['iwindoff']=0
	dflt['iwind_form']=-1

	dflt['nonhydro']=0
	dflt['ihdif']=0
	dflt['ibcc']=1
	dflt['bfric']=-1
	dflt['ihconsv']=0
	dflt['isconsv']=0
	
	dflt['ibcc_mean']=0
	dflt['inu_st']=0
	dflt['icst']=1
	
	dflt['ishapiro']=1
	dflt['ihorcon']=0
	dflt['indvel']=0
	dflt['inter_mom']=0

	## Gloabl output
	dflt['elev61']=1
	dflt['pres61']=0
	dflt['airt61']=0
	dflt['shum61']=0
	dflt['srad61']=0
	dflt['flsu61']=0
	dflt['fllu61']=0
	dflt['radu61']=0
	dflt['radd61']=0
	dflt['flux61']=0
	dflt['evap61']=0
	dflt['prcp61']=0
	dflt['bdrc61']=0
	dflt['wind62']=0
	dflt['wist62']=0
	dflt['dahv62']=1
	dflt['vert63']=0
	dflt['temp63']=0
	dflt['salt63']=0
	dflt['conc63']=0
	dflt['tdff63']=0
	dflt['vdff63']=0
	dflt['kine63']=0
	dflt['mixl63']=0
	dflt['qnon63']=0
	dflt['hvel64']=0

#	## Non standard ouput
	dflt['hvel67']=0
	
	modules=['GEN','AGE']
	for mo in modules:
		for n in range(0,10):
			dflt[mo+'_'+str(n+1)]=0



	return dflt


class param(object):
	def __init__(self):
		self.tmpl=tmpl
		self.opt=default_value()

	def pass_values(self,new_dict):
		for key in self.opt:
			if key in new_dict:
				self.opt[key]=new_dict[key]

		# Re-arange the timing
		unit=new_dict.get('unit','second')

		if unit=='hour':
			fac=3600
		elif unit=='second':
			fac=1
		elif unit=='minute':
			fac=60
		elif unit=='day':
			fac=3600*24.

		dt=self.opt['dt']
		if 'output dt' in new_dict:
			self.opt['nspool'] = int((new_dict.pop('output dt')*fac)/dt)

		if 'file length' in new_dict:
			self.opt['ihfskip']=int((new_dict.pop('file length')*fac)/dt)

		if 'station dt' in new_dict:
			self.opt['nspool_sta']=int((new_dict.pop('station dt')*fac)/dt)

		if 'hotstart dt' in new_dict:
			self.opt['hotout_write']=int((new_dict.pop('hotstart dt')*fac)/dt)

		if 'ramp' in new_dict:
			self.opt['dramp']=(new_dict.pop('ramp')*fac)/(3600.*24.)

		

		# get tracers
		modules=['gen_tracer','age_tracer']
		for mo in modules:
			N=new_dict.get(mo,0)
			if N>0:
				for n in range(0,10):
					if n<=N-1:
						self.opt[mo[0:3].upper()+'_'+str(n+1)]=1
					else:
						self.opt[mo[0:3].upper()+'_'+str(n+1)]=0


	def write(self,fileout='param.in'):



		cmdblank = Template(filename=self.tmpl)
		buf = StringIO()
		stri='ctx = Context(buf,'
		for key in self.opt:
			stri=stri+key+'='+str(self.opt[key])+','

		stri=stri[0:-1]+')'

		exec stri


		cmdblank.render_context(ctx)

		
		a=open(fileout,'w')
		
		a.write(buf.getvalue())
		a.close()
