# Metric for run72a
import station_db
import obs_links
import plot_default_formats as pd
import schism_postprocess
import unit_conversions
import vtools.data.vtime
import vtools.functions.api
import yaml
import argparse
import re
import read_ts
import station_extractor
import matplotlib.pyplot as plt
import numpy as np
import datetime
import os
import logging


def create_arg_parser():
    parser = argparse.ArgumentParser()
    # Read in the input file
    parser = argparse.ArgumentParser(description="Create metrics plots in a batch mode")
    parser.add_argument(dest='main_inputfile', default=None,
                        help='main input file name')
    return parser


def init_logger():
    logging.basicConfig(level=logging.INFO, filename="metrics.log",filemode='w')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    #%(asctime)s - %(name)s - %(levelname)s
    formatter = logging.Formatter('%(message)s')
    console.setFormatter(formatter)
    filer = logging.FileHandler('metrics_errors.log')
    filer.setLevel(logging.WARNING)
    formatter2 = logging.Formatter('%(name)s - %(message)s')
    filer.setFormatter(formatter2)
    logging.getLogger('').addHandler(filer)
    logging.getLogger('').addHandler(console)


def metrics(params):
    init_logger()
    logger = logging.getLogger("metrics")

    working_dir = params["data dir"]
    time_basis = params["time basis"]

    start_avg = params["start avg"]
    end_avg = params["end avg"]
    start_inst = params["start inst"]
    end_inst = params["end inst"]
    labels_base = params["labels"]

    db_stations = station_db.StationDB(params["stations csv"])

    list_obs_fname = params["obs links csv"]
    db_obs = obs_links.ObsLinks(list_obs_fname)
    pad = datetime.timedelta(days=3)

    max_gap_to_fill = vtools.data.vtime.time_interval(hours=1)
    variable = params["variable"]
    if variable == 'flow':
        station_in = os.path.join(working_dir, "flowlines.yaml")
        staout = station_extractor.flow_extractor(station_in,
                                                  working_dir,
                                                  time_basis)
        staout.flow_fname = "flux.dat"
    else:
        station_in = os.path.join(working_dir, "station.in")
        staout = station_extractor.station_extractor(station_in,
                                                     working_dir,
                                                     time_basis)
    use_alias = False
    pd.set_color_cycle_dark2()
    mp = schism_postprocess.metrics_plot()
    for station in staout.stations:
        station_id = station["name"]
        logger.info("Processing station: %s" % station_id)
        alias = db_stations.alias(station_id)
        station_work = station_id
        if alias is None:
            logger.error("Station is not in the station database: %s" % station_id)
            continue
        long_name = db_stations.name(station_id)

        # Simulation
        if variable == 'flow':
            ts_sim = staout.retrieve_ts(station_id)
            data_expected = db_stations.station_attribute(station_id,"Flow")
        else:
            if variable == 'elev':
                ts_sim = staout.retrieve_ts(variable, name=station_id)
                data_expected = db_stations.station_attribute(station_id,"Stage")
            elif variable in ['salt', ]:
                vert_pos = station["vert_pos"]
                data_expected = db_stations.station_attribute(station_id,"WQ")
                ts_sim = staout.retrieve_ts(variable, name=station_id, depth=vert_pos)
        if ts_sim is None:
            logger.warning("This station is not in staout: %s %s: " %(station_id, long_name))
            continue
        if np.all(np.isnan(ts_sim.data)):
            logger.warning("All simulated values are nan.")
            continue


        # Observation
        if variable == "salt":
            obs_fname = db_obs.filename(station_id, variable, vert_pos=vert_pos)
        else:
            obs_fname = db_obs.filename(station_id, variable)
        if obs_fname is None and use_alias:
            if variable == 'salt':
                obs_fname = db_obs.filename(alias, variable, vert_pos=vert_pos)
            else:
                obs_fname = db_obs.filename(alias, variable)
            station_work = alias

        if alias != station_id:
            alias_use = alias + " " + long_name[:32]
        else:
            alias_use = long_name[:32]

        if obs_fname is None:
            expectstr = "(Omission)" if data_expected else "(Data not expected)"
            level = logging.WARNING if data_expected else logging.DEBUG
            logger.log(level,"%s No %s data link listing for: %s (%s)" % (expectstr,variable,station_id,alias_use))
            continue
        else:
            if not data_expected:
                logger.warning("File link %s found for station %s but station not expected to have data for variable %s" % (obs_fname,station_id,variable))
        logger.info("Observed file for id %s: %s" % (station_id,obs_fname))
        if not os.path.exists(obs_fname):
            logger.error("Observation file path not found on file system: %s" % obs_fname)
            continue

        ts_obs = read_ts.read_ts(obs_fname, start_avg - pad, end_avg + pad)

        if ts_obs is None:
            logger.warning("This obs file does not cover the period: %s" % obs_fname)
            continue

        # Fill gaps
        if np.any(np.isnan(ts_obs.data)):
            logger.debug("Filling gaps ...")
            max_gap = int(max_gap_to_fill.total_seconds() / ts_obs.interval.total_seconds())
            if max_gap == 0:
                max_gap += 1
            ts_obs = vtools.functions.api.interpolate_ts_nan(ts_obs, max_gap=max_gap)

        # agency = db_obs.get_agency(station_id, variable)
        obs_unit = db_obs.unit(station_work, variable)
        logger.debug("obs_unit %s" % obs_unit)
        ts_obs.props["unit"] = obs_unit
        if ts_obs.props["unit"] == 'ft':
            unit_conversions.ft_to_m(ts_obs)
        elif ts_obs.props["unit"] == "cfs":
            logger.debug("Convert cfs to cms....")
            unit_conversions.cfs_to_cms(ts_obs)
        elif ts_obs.props["unit"] == "ec":
            logger.debug("Convert ec to psu...")
            unit_conversions.ec_psu_25c(ts_obs)
        datum_adj = db_obs.datum_adjust(station_work, variable)


        if datum_adj != 0.:
            ts_obs += datum_adj
        datum = db_obs.vdatum(station_work, variable)
        adj = 0.
        if variable == 'elev' and datum == 'STND':
            tss = schism_postprocess.window_list_timeseries((ts_obs, ts_sim))
            adj = np.average(tss[1].data) - np.average(tss[0].data)
            ts_obs += adj

        label_obs = labels_base[0]
        if adj != 0.:
            if adj > 0.:
                label_obs += " + %g" % adj
            else:
                label_obs += u" \u2212 %g" % (-adj)
        labels = [label_obs,]
        labels.extend(labels_base[1:])

        mp.plot_metrics((ts_obs, ts_sim), labels=labels)

        mp.axes_inst.set_xlim(start_inst, end_inst)
        pd.set_xaxis_dateformat(mp.axes_inst, "%m/%d/%y")
        pd.rotate_xticks(mp.axes_inst, 30)

        mp.axes_filtered.set_xlim(start_avg, end_avg)
        pd.set_xaxis_dateformat(mp.axes_filtered, "%m/%d/%y")
        pd.rotate_xticks(mp.axes_filtered, 30)
        mp.axes_inst.legend(prop={'size': 12})
        if variable in ["salt",]:
            if vert_pos == 0:
                title = "%s, Surface, (%s)" % (long_name, alias)
                fout_name = "%s_%s_surface.png" % (variable, alias)
            elif vert_pos == 1:
                title = "%s, Bottom, (%s)" % (long_name, alias)
                fout_name = "%s_%s_bottom.png" % (variable, alias)
        else:
            title = "%s (%s)" % (long_name, alias)
            fout_name = "%s_%s.png" % (variable, alias)
        mp.set_title(title)
        plt.savefig(fout_name, dpi=300)
    print "See metrics_errors.log for errors such as missed stations"


def get_params(fname):
    with open(fname, 'r') as fin:
        params = yaml.load(fin)
        for k, v in params.iteritems():
            if isinstance(v, datetime.date):
                params[k] = datetime.datetime.fromordinal(v.toordinal())
            if isinstance(v, str):
                v = v.split()
                if len(v) == 1:
                    params[k] = v[0]
                else:
                    params[k] = v
    return params


if __name__ == "__main__":
    parser = create_arg_parser()
    args = parser.parse_args()

    params = get_params(args.main_inputfile)
    metrics(params)
