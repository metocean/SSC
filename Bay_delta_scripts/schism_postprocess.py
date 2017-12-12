""" Routines to post-process output requests
"""
# Contact:
# 2014, California Department of Water Resources
# Kijin Nam, Bay-Delta Office
# knam@water.ca.gov

import unit_conversions
import plot_default_formats as pd
from vtools.data.timeseries import *
from vtools.data.vtime import *
from vtools.functions.api import *
from vtools.functions.skill_metrics  import *
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl
import numpy as np
import datetime
import collections
import os

__all__ = ["plot_multiple_timeseries", "plot_instantaneous", "metrics_plot"]

def window_list_timeseries(list_timeseries, window=None, pad=None):
    """ Window multiple time series with some safety to check for None
    """
    if window is None:
        window = range_intersect(*list_timeseries[:2])
    window_to_cut = list(copy.deepcopy(window))
    if not pad is None:
        window_to_cut[0] -= pad
        window_to_cut[1] += pad
    for timeseries in list_timeseries:
        if not timeseries is None:
            if window_to_cut[0] < timeseries.start:
                window_to_cut[0] = timeseries.start
            if window_to_cut[1] > timeseries.end:
                window_to_cut[1] = timeseries.end
    list_timeseries_new = list()
    for timeseries in list_timeseries:
        if timeseries is None:
            list_timeseries_new.append(None)
        else:
            list_timeseries_new.append(timeseries.window(*window_to_cut))
    return list_timeseries_new


def plot_instantaneous(axes, list_timeseries, window=None, labels=None):
    """ Plot a instantaneous plot
    """
    if window is not None:
        window_to_cut = window
        list_timeseries_to_plot = window_list_timeseries(list_timeseries,
                                                         window_to_cut)
    else:
        list_timeseries_to_plot = list_timeseries
    plot_multiple_timeseries(axes, *list_timeseries_to_plot, labels=labels)
    axes.set_xlim(window)
    if labels:
        axes.legend()
        axes.legend(prop={'size': 12})
    set_xylabels(axes, list_timeseries_to_plot)

def filter_list_timeseries(list_timeseries, window=None):
    """
    """
    if window is None:
        window = range_intersect(list_timeseries[0], list_timeseries[1])
    list_timeseries_windowed = \
        window_list_timeseries(list_timeseries, window)
    # padding = datetime.timedelta(days=3)
    # print [(ts.start, ts.end) for ts in list_timeseries]
    list_timeseries_filtered = list()
    for timeseries in list_timeseries_windowed:
        timeseries_filtered = lanczos_filter(timeseries)
        timeseries_filtered = timeseries_filtered.window(*window)
        list_timeseries_filtered.append(timeseries_filtered)
    return list_timeseries_filtered

# def plot_filtered(axes, list_timeseries, window=None, labels=None):
#     if not window:
#         window = range_intersect(list_timeseries[0], list_timeseries[1])
#     padding = datetime.timedelta(days=3)
#     window_to_filter = (window[0] - padding, window[1] + padding)
#     list_timeseries_windowed = \
#         window_list_timeseries(list_timeseries, window_to_filter)
#     list_timeseries_filtered = list()
#     for timeseries in list_timeseries_windowed:
#         timeseries_filtered = lanczos_filter(timeseries)
#         list_timeseries_filtered.append(timeseries_filtered)

#     plot_multiple_timeseries(axes, *list_timeseries_filtered, labels=labels)
#     axes.set_xlim(window)
#     # set_xylabels(axes, list_timeseries_filtered)


def align_intersect_timeseries(base, target):
    """ Align target time series to the base one.
        It assumes a regular time series for the base.
    """
    from vtools.functions.api import interpolate_ts

    if base.interval < target.interval:
        raise ValueError("Untested with coarser target and windowing may be wrong")
    window_common = range_intersect(base, target)
    start = align(window_common[0], base.interval, 1)
    end = align(window_common[1], base.interval, -1)
    base_windowed = base.window(start, end)
    target_aligned = interpolate_ts(target,
                                    base_windowed.times,
                                    method=LINEAR)
    props = copy.deepcopy(target.props)
    target_aligned._props = props
    return base_windowed, target_aligned


def plot_scatter(axes, ts1, ts2, labels=None):
    """ Plot a scatter plot with a regression line
    """
    # if len(ts1_common) != len(ts2_common):
    #     raise ValueError("Length of the two time series must be the same.")
    d1 = ts1.data
    d2 = ts2.data

    # Remove nan pairs
    nonnan_flag = np.logical_not(np.logical_or(np.isnan(d1), np.isnan(d2)))
    d1 = d1[nonnan_flag]
    d2 = d2[nonnan_flag]

    # Draw
    h = axes.scatter(d1, d2)
    plt.setp(h, alpha=0.15, edgecolor='grey',
             facecolor=mpl.rcParams['axes.color_cycle'][0])

    if len(d1) < 2: return
    # Regression
    model, resid = np.linalg.lstsq(
        np.vstack([d1, np.ones(len(d1))]).T,
        d2)[:2]
    # r2 = 1. - resid / (len(d2) * d2.var())
    x = np.array([min(d1), max(d1)])
    y = model[0] * x + model[1]
    l, = axes.plot(x, y, color=mpl.rcParams['axes.color_cycle'][1])
    # Text info
    if model[1] >= 0.:
        eqn = "Y=%.3f*X+%.3f" % (model[0], model[1])
    # Calculate the linear regression
    else:
        eqn = "Y=%.3f*X-%.3f" % (model[0], -model[1])
    axes.legend([l,], [eqn,], loc='upper left', prop={'size': 10})

    xlim = axes.get_xlim()
    ylim = axes.get_ylim()
    common_lim = (min(xlim[0], ylim[0]), max(xlim[1], ylim[1]))
    axes.set_xlim(*common_lim)
    axes.set_ylim(*common_lim)
    axes.set_aspect('equal')
    pd.rotate_xticks(axes, 20)
    axes.tick_params(axis='x', pad=10)

    if labels is None:
        labels = ['Obs', 'Sim']
    if 'unit' in ts1.props.keys():
        unit_x = ts1.props['unit']
    else:
        unit_x = None
    if unit_x == 'meter':
        unit_x = 'm'
    if unit_x is not None:
        xlabel = labels[0] + " (%s)" % unit_x
        axes.set_xlabel(xlabel)
        ylabel = labels[1] + " (%s)" % unit_x
        axes.set_ylabel(ylabel)


def lanczos_filter(timeseries, cutoff_period='40hour'):
    """ Filter a time series using cosine-Lanczos filter.

        Parameters
        ----------
        timeseries:
        A time series to filter

        Returns
        -------
        A filtered time series
    """
    if timeseries is None:
        return None
    timeseries_fltrd = cosine_lanczos(timeseries,
                                      cutoff_period=cutoff_period)
    timeseries_fltrd.props['filtered'] = 'cosine_lanczos'
    #if np.all(np.isnan(timeseries_fltrd.data)):
    #    timeseries_fltrd = None
    return timeseries_fltrd


def set_xylabels(axes, list_timeseries):
    pd.auto_ylabels(axes, list_timeseries)
    xlim = axes.get_xlim()
    window = xlim[1] - xlim[0]
    if window < 5.:
        pd.set_xaxis_dateformat(axes, date_format="%m/%d/%y %H:%M", pad=10)
        pd.rotate_xticks(axes, angle=20)
    elif window < 15.:
        pd.set_xaxis_day(axes, rotate=20, pad=10)
    elif window > 160.:
        pd.set_xaxis_month(axes, rotate=20, pad=10)
    else:
        pd.set_xaxis_dateformat(axes, date_format="%m/%d/%y", pad=10)
        pd.rotate_xticks(axes, angle=20)


def plot_multiple_timeseries(*args, **kwargs):
    """ Draw multiple time series in one plot without any format

        Parameters
        ----------
        args
            matplotlib axis and the list of vtools time series

        Returns
        -------
        list
            handles of lines
    """
    lines = list()
    if isinstance(args[0], mpl.axes.Axes):
        axes = args[0]
        args = args[1:]
    else:
        axes = mpl.pyplot.subplot(111)
    n_timeseries = 0
    for arg in args:
        if isinstance(arg, TimeSeries) or arg is None:
            n_timeseries += 1
    if 'labels' in kwargs.keys():
        labels = kwargs['labels']
        if labels is not None and len(labels) != n_timeseries:
            raise ValueError("# of labels does not match to # of time series")
        del kwargs['labels']
    else:
        labels = None
    i_line = 0
    for arg in args:
        if arg is None:
            axes.plot(None)
            lines.append(None)
            i_line += 1
        elif isinstance(arg, TimeSeries):
            if np.all(np.isnan(arg.data)):
                axes.plot(None)
                lines.append(None)
            else:
                if labels is not None:
                    line, = axes.plot(arg.times, arg.data,
                                      label=labels[i_line], **kwargs)
                else:
                    line, = axes.plot(arg.times, arg.data, **kwargs)
                lines.append(line)
            i_line += 1
        else:
            raise ValueError("This function supports Vtools time series only.")
    return lines


def calculate_metrics(obs, sim, sim_shifted):
    """ Calculate metric numbers and show it at the bottom of the figure

        Parameters
        ----------
        obs
            Observation time series.
        sim
            Simulated time series. This one has to have the same time interval
            os that of the observation time series
        sim_shifted
            A shifted simulated time series by amount of lag

        Returns
        -------
        list
            RMSE, Bias, NSE, and R
    """
    # Metrics
    rmse_ = rmse(sim, obs)
    bias = median_error(sim_shifted, obs)
    nse = skill_score(sim_shifted, obs)
    corr = corr_coefficient(sim_shifted, obs)
    return [rmse_, bias, nse, corr]


class MetricsPlot(object):
    """ A class to plot metrics plots
    """
    def __init__(self):
        """
        """
        self.fig = plt.gcf()
        self.axes_inst = None
        self.axes_filtered = None
        self.axes_scatter = None
        self.ax_inst = None
        self.ax_fltrd = None
        self.ax_scatter = None
        self.max_shift = minutes(120.)

    def _set_metrics_axes(self):
        """ Create a set of grids for metrics plot
        """
        grids = gridspec.GridSpec(2, 2,
                                  width_ratios=[1.6, 1], height_ratios=[1, 1])
        grid_inst = grids[0, 0:2]
        self.axes_inst = self.fig.add_subplot(grid_inst)
        grid_fltrd = grids[1, 0]
        self.axes_filtered = self.fig.add_subplot(grid_fltrd)
        grid_scatter = grids[1, 1]
        self.axes_scatter = self.fig.add_subplot(grid_scatter)
        grids.update(top=0.93, bottom=0.15, right=0.88,
                     hspace=0.4, wspace=0.8)

    def plot_metrics(self, list_timeseries, window_inst=None, window_avg=None,
                     labels=None, adj=dict()):
        """ Draw a generic metric plot without axis labels and axis treatment.
            It does not put a title or create a graphic file.

            Parameters
            ----------
            list_timeseries: array-like of vtools.data.timeseries
                A list of time series to plot

            window_inst: tuple of datetime.datetime (size 2)
                time window for the instantaneous value plot

            window_avg: tuple of datetime.datetime (size 2)
                time window for the filtered value plot.
                    This window is used to calculate metrics
                    and to plot a scatter plot.
            labels: array-like of str
                A list of names of the time series for legend.
            adj: dict, optional
                Adjustment made in the data
        """
        # Clear up
        self.fig.clf()
        self.fig.set_dpi(300)
        fig_size = pd.metrics_fig_size
        self.fig.set_size_inches(fig_size[0], fig_size[1])
        self._set_metrics_axes()

        # First: Instantaneous plot
        # plot_instantaneous(self.axes_inst, list_timeseries,
        #                    window_inst, labels=labels)
        plot_multiple_timeseries(self.axes_inst, *window_list_timeseries(list_timeseries, window_inst), labels=labels)
        pd.set_dual_axes(self.axes_inst, list_timeseries[0])

        # Second: Average plot
        # plot_filtered(self.axes_filtered, list_timeseries, window_avg)
        # list_timeseries_filtered = filter_list_timeseries(list_timeseries, window_avg)
        list_timeseries_filtered = [lanczos_filter(ts) for ts in window_list_timeseries(list_timeseries, pad=time_interval(days=3))]

        plot_multiple_timeseries(self.axes_filtered,
                                 *window_list_timeseries(list_timeseries_filtered, window_avg))
        ylim_inst = self.axes_inst.get_ylim()
        # ylim_filtered = list(self.axes_filtered.get_ylim())
        # if ylim_filtered[0] > ylim_inst[0]:
        #     ylim_filtered[0] = ylim_inst[0]
        # if ylim_filtered[1] < ylim_inst[1]:
        #     ylim_filtered[1] = ylim_inst[1]
        # self.axes_filtered.set_ylim(ylim_filtered)
        pd.set_dual_axes(self.axes_filtered, list_timeseries_filtered[0])
        # pd.set_nice_tick_intervals(self.axes_filtered, 5)


        # Third: Scatter plot with first two time series
        ts1 = list_timeseries[0]
        ts2 = list_timeseries[1]
        ts2_orig = ts2


        if window_avg is not None:
            pass
            #ts1_window = ts1.window(*window_avg)   # obs
            #ts2_window = ts2.window(*window_avg)   # model
        else:
            window_avg = (ts1.start, ts1.end)
            #ts1_window = ts1_window
            #ts2_window = ts2.window(*window_avg)

        lag = calculate_lag(ts1, ts2, window_avg,
                            max_shift=self.max_shift,
                            period=datetime.timedelta(minutes=12.42*60.))

        if lag == datetime.timedelta(minutes=0):
            ts2_shifted = ts2_orig
        elif lag > self.max_shift or lag < -self.max_shift:
            ts2_shifted = ts2_orig
            lag = None
        else:
            ts2_shifted = ts2_orig.shift(-lag)


        ts1_common, ts2_shifted_common = align_intersect_timeseries(ts1.window(*window_avg), ts2_shifted)

        # The shifted ts2 and the unshifted may have a different common window with ts1.
        # Here we feed in ts_common, and
        # the function align_intersect can only make the window smaller
        ts1_common, ts2_common = align_intersect_timeseries(ts1_common,ts2)
        # Now record the possibly smaller window. ts1_common and ts2_common are already reduced
        # to this window, but ts2_shifted_common may need it
        possibly_smaller_window = (ts1_common.start,ts1_common.end)
        ts2_shifted_common = ts2_shifted_common.window(*possibly_smaller_window)

        plot_scatter(self.axes_scatter, ts1_common, ts2_shifted_common)
        metrics = calculate_metrics(ts1_common, ts2_common, ts2_shifted_common)
        if 'obs_adj' in adj.keys():
            metrics[1] = None
        metrics.append(lag)
        unit = ts1.props['unit']
        str_metrics = self.gen_metrics_string(metrics, unit)
        if len(str_metrics) != 0:
            center = 0.5
            top = -1.7
            self.axes_inst.text(
                center, top, str_metrics,
                horizontalalignment='center',
                verticalalignment='top',
                transform=self.axes_inst.transAxes)

    def set_title(self, title):
        self.axes_inst.set_title(title)

    def save_fig(self, fpath_out):
        plt.savefig(fpath_out)

    @classmethod
    def gen_metrics_string(cls, metrics, unit):
        # Metrics
        metric_text = str()
        if unit is not None:
            metric_text += "RMSE=%.3f %s   " % (metrics[0], unit)
        else:
            metric_text += "RMSE=%.3f      " % (metrics[0])
        lag = metrics[-1]
        if lag is not None:
            seconds_ = lag.total_seconds()
            metric_text += "Lag=%d min   " % (seconds_ / 60)
        else:
            metric_text += "Lag=N/A   "
        if metrics[1] is not None:
            metric_text += r"Bias$_\phi$=%.3f   " % metrics[1]
        else:
            metric_text += r"Bias$_\phi$=N/A    "
        metric_text += r"NSE$_\phi$=%.3f   " % metrics[2]
        metric_text += r"R$_\phi$=%.3f   " % metrics[3]
        return metric_text


def metrics_plot():
    """ Create an instance of MetrcisPlot and return
    """
    mp = MetricsPlot()
    return mp
