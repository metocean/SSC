
from separate_species import *
from vtools.functions.api import *
from vtools.data.timeseries import *
import schism_setup
import numpy as np
import datetime
import struct
import matplotlib.pyplot as plt

################# command line application #####################
def create_arg_parser():
    import textwrap
    parser = argparse.ArgumentParser(
      formatter_class=argparse.RawDescriptionHelpFormatter,
      epilog = textwrap.dedent(
      """ 
         ============== Example ==================
      > gen_elev2D.py --stime=2009-03-12 --etime=2010-01-01 9415020_gageheight.csv 9413450_gageheight.csv
      """),
      description= """Script to create elev2D.th boundary condition from Point Reyes and Monterey NOAA file"""
    )
    parser.add_argument('--stime',default = None, required=False, help = "Start time in ISO-like format 2009-03-12T00:00:00. Time part and 'T' are optional.")
    parser.add_argument('--etime', default = None, required=False, help = 'End time.')
    parser.add_argument('--hgrid',default = 'hgrid.gr3',required=False,help = 'Name of hgrid file if not hgrid.gr3')    
    parser.add_argument('--outfile', default = 'elev2D.th', required = False, help = 'Name of output file if not elev2D.th')   
    parser.add_argument('pt_reyes',default = None,  help = 'Pt Reyes data file, must have a buffer of 16 days at either end of series')
    parser.add_argument('monterey',default = None,  help = 'Monterey data file, must have a buffer of 16 days at either end of series')
    return parser



def main():    
    parser = create_arg_parser()
    args = parser.parse_args()
    monterey_fpath = args.monterey
    pt_reyes_fpath = args.pt_reyes
    hgrid_fpath = args.hgrid
    fpath_out = args.outfile

    
    # UTM positions of Point Reyes, Monterey, SF
    pos_pr = np.array([502195.03, 4205445.47])
    pos_mt = np.array([599422.84, 4051630.37])
    pos_sf = np.array([547094.79, 4184499.42])

    var_subtidal = np.array([0.938, 0.905, 0.969])  # pr, mt, sf
    var_semi = np.array([0.554, 0.493, 0.580])

    tangent = np.array([1, -1])  # Assume 45 degree from north-west to south-east
    tangent = tangent / np.linalg.norm(tangent)  # Normalize
    normal = np.array([tangent[1], -tangent[0]]) # Rotate 90 cw to get normal vec
    print "tangent:", tangent
    print "normal:", normal

    mt_rel = pos_mt - pos_pr
    x_mt = np.dot(tangent, mt_rel) # In pr-mt direction
    y_mt = np.dot(normal, mt_rel)  # Normal to x-direction to the ocean

    # Grid
    s = schism_setup.load_gr3(hgrid_fpath)
    ocean_boundary = s.mesh.boundaries[0]  # First one is ocean
    buf = days(16)
    sdate=dtm.datetime(*map(int, re.split('[^\d]', args.stime))) # convert start time string input to datetime
    if args.etime:
        edate=dtm.datetime(*map(int, re.split('[^\d]', args.etime))) # convert start time string input to datetime
        bufend = edate + buf
    else:
        edate = None
        bufend = None

    # Data
    #start = datetime.datetime(2009, 3, 12)
    #end = datetime.datetime(2011, 1, 1)     
    #interval = datetime.timedelta(minutes=6)
    print "Reading Point Reyes"
    pt_reyes = read_ts(pt_reyes_fpath,sdate - buf, bufend,force_regular=True)
    ts_pr_subtidal, ts_pr_diurnal, ts_pr_semi, noise = separate_species(pt_reyes)
    del noise

    
    print "Reading Monterey"
    monterey = read_ts(monterey_fpath,sdate - buf, bufend,force_regular=True)
    if pt_reyes.interval != monterey.interval:
        raise ValueError("Point Reyes and Monterey time step must be the same in gen_elev2D.py")

    ts_mt_subtidal, ts_mt_diurnal, ts_mt_semi, noise = separate_species(monterey)
    del noise

    dt = monterey.interval.seconds
    
    print "Done Reading"

    print "Interpolating Point Reyes"
    ts_pr_subtidal = ts_pr_subtidal.window(sdate,edate) # interpolate_ts(ts_pr_subtidal.window(sdate,edate),step)
    ts_pr_diurnal = ts_pr_diurnal.window(sdate,edate)  #interpolate_ts(,step)
    ts_pr_semi = ts_pr_semi.window(sdate,edate)   #interpolate_ts(ts_pr_semi.window(sdate,edate),step)

    print "Intepolating Monterey"
    ts_mt_subtidal = ts_mt_subtidal.window(sdate,edate) #interpolate_ts(ts_mt_subtidal.window(sdate,edate),step)
    ts_mt_diurnal = ts_mt_diurnal.window(sdate,edate)   #interpolate_ts(ts_mt_diurnal.window(sdate,edate),step)
    ts_mt_semi =    ts_mt_semi.window(sdate,edate)      #interpolate_ts(ts_mt_semi.window(sdate,edate),step)


    # Grid
    boundaries = s.mesh.nodes[ocean_boundary.nodes]
    pos_rel = boundaries[:, :2] - pos_pr
    # x, y in a new principal axes
    x = np.dot(pos_rel, tangent.reshape((2, -1)))
    y = np.dot(pos_rel, normal.reshape((2, -1)))
    theta_x = x / x_mt
    theta_x_comp = 1. - theta_x
    theta_y = y / y_mt
    theta_y_comp = 1. - theta_y
    print "y", y.reshape((1,-1))
    print "y_mt", y_mt
    print "theta_x:", theta_x.reshape((1,-1))
    print "theta_y:", theta_y.reshape((1,-1))

    var_y = (theta_y_comp * var_semi[0] + theta_y * var_semi[1])

    adj_subtidal_mt = 0.08  # Adjustment in Monterey subtidal signal
    scaling_diurnal_mt = 0.95 # Scaling of Monterey diurnal signal (for K1/Q1)
    adj_subtidal_mt = 0.  # Adjustment in Monterey subtidal signal
    scaling_diurnal_mt = 1. # Scaling of Monterey diurnal signal (for K1/Q1)


    with open(fpath_out, 'wb') as f:
        for i in xrange(len(ts_pr_semi)):
            t = float(dt*i)
            # semi-diurnal
            # Scaling
            pr = ts_pr_semi[i].value
            mt = ts_mt_semi[i].value

            if pr is np.nan or mt is np.nan:
                raise ValueError("One of values is numpy.nan.")

            eta_pr_side = var_y / var_semi[0] * pr
            eta_mt_side = var_y / var_semi[1] * mt
            eta = eta_pr_side * theta_x_comp + eta_mt_side * theta_x

            # diurnal
            # Interpolate in x-direction only to get a better phase
            pr = ts_pr_diurnal[i].value
            mt = ts_mt_diurnal[i].value * scaling_diurnal_mt
            eta += pr * theta_x_comp + mt * theta_x

            # Subtidal
            # No phase change in x-direction. Simply interpolate in y-direction.
            pr = ts_pr_subtidal[i].value
            mt = ts_mt_subtidal[i].value + adj_subtidal_mt
            eta += pr * theta_y_comp + mt * theta_y

            # Write
            fmt = "f"
            buf = struct.pack(fmt, t)
            f.write(buf)
            fmt = "f" * len(eta)
            buf = struct.pack(fmt, *eta)
            f.write(buf)


if __name__ == "__main__":
    main()
