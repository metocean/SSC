#!/usr/bin/env python2.7
import sys,os
sys.path.append(os.path.join(os.path.dirname(__file__),'Bay_delta_scripts'))
sys.path.append(os.path.join(os.path.dirname(__file__),'func')) 
import yaml
import argparse
from param import param
import subprocess
import signal
import time
from schismIO import schismIO
from power import power,get_areas
import numpy as np
from export_nc import export_nc
from polygons import *
import timeit

root='/home/remy/Buisness/0336_SSC_tides/test_run/input/' # My root directory

sc=schismIO(root) # this will combine he file as it run
pw=power(sc) # this wil get the power after 1 tidal cycle

pw.get_power(13)
pw.export_nc(13,outdir='/home/remy/Buisness/0336_SSC_tides/test_run/output/')
