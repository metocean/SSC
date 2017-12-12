import matplotlib.pyplot as plt
from read_ts import read_cdec
from vtools.data.api import rts
import numpy as np
import sys
from scipy.signal import medfilt
from scipy.stats.mstats import mquantiles

def med_outliers(ts,level=4,scale = None,\
                 filt_len=7,range=(None,None),
                 quantiles = (0.25,0.75),
                 copy = True):
    """
    Detect outliers by running a median filter, subtracting it
    from the original series and comparing the resulting residuals
    to a global robust range of scale (the interquartile range).
    Individual time points are rejected if the residual at that time point is more than level times the range of scale. 

    The original concept comes from Basu & Meckesheimer (2007)
    although they didn't use the interquartile range but rather
    expert judgment. To use this function effectively, you need to
    be thoughtful about what the interquartile range will be. For instance,
    for a strongly tidal flow station it is likely to 
    
    level: Number of times the scale or interquantile range the data has to be
           to be rejected.

    scale: Expert judgment of the scale of maximum variation over a time step.
           If None, the interquartile range will be used. Note that for a 
           strongly tidal station the interquartile range may substantially overestimate the reasonable variation over a single time step, in which case the filter will work fine, but level should be set to 
           a number (less than one) accordingly.

    filt_len: length of median filter, default is 5
    
    quantiles : tuple of quantiles defining the measure of scale. Ignored
          if scale is given directly. Default is interquartile range, and
          this is almost always a reasonable choice.

    copy: if True, a copy is made leaving original series intact

    You can also specify rejection of values based on a simple range

    Returns: copy of series with outliers replaced by nan
    """
    import warnings
    ts_out = ts.copy() if copy else ts
    warnings.filterwarnings("ignore")
    if range[0]:
        ts_out.data[ts_out.data < range[0]] = np.nan
        
    if range[1]:
        ts_out.data[ts_out.data> range[1]] = np.nan

        
    filt = medfilt(ts_out.data, kernel_size=filt_len)
    res = ts_out.data - filt
    if not scale:
        low,high = mquantiles(filt[~ np.isnan(filt)],quantiles)
        scale = high - low 
    print np.sum(np.isnan(ts_out.data))
        
    ts_out.data[(res > level*scale) | (res < -level*scale)]= np.nan
    warnings.resetwarnings()
    print np.sum(np.isnan(ts_out.data))
    filt = rts(filt,ts.start,ts.interval)
    return ts_out, filt




if __name__ == '__main__':
    # Just an example
    station = sys.argv[1]
    ts = read_cdec("cdec_download/%s.csv"%station,start=None,end=None)

    filt = medfilt(ts.data, kernel_size=5)
    ts,filt = med_outliers(ts,quantiles=[0.2,0.8],range=[120.,None],copy = False)

    plt.plot(ts.times,ts.data,label="data")
    plt.plot(ts.times,filt.data,label="filt")
    plt.plot(ts.times,ts.data-filt.data,label="res")
    plt.legend()
    plt.show()