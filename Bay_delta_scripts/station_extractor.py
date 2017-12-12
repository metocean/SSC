"""
Routines to read in SCHISM STAOUT_* and flux.dat files.
"""
#
# Copyright:
# 2014, California Department of Water Resources
# Kijin Nam, Bay-Delta Office
# knam@water.ca.gov

import yaml
from vtools.data.timeseries import *
import numpy as np
import collections
import math
import datetime
import re
import csv
import sys
import os

__all__ = ['station_extractor', 'flow_extractor']

# Ordered Dict YAML
# From http://stackoverflow.com/questions/5121931/in-python-how-can-you-load-yaml-mappings-as-ordereddicts
def dict_representer(dumper, data):
    return dumper.represent_mapping(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, data.iteritems())


def dict_constructor(loader, node):
    return collections.OrderedDict(loader.construct_pairs(node))


class StaoutReader():
    """ SCHISM STAOUT reader
    """
    def __init__(self, working_dir, time_basis):
        """ Constructor

            Parameters
            ----------

            working_dir: A working diretory where STAOUT_* reside

            time_basis: A time basis of the STAOUT_*
        """
        self._working_dir = working_dir
        self._time_basis = time_basis
        self._elev = None
        self._salt = None
        self._salt_fname = 'staout_6'
        self._elev_fname = 'staout_1'
        self._items = None
        self._stations = None

    @property
    def stations(self):
        return self._stations

    def read_station_in(self, fpath):
        """ Read in 'station.in.'
            'station.in' should contains station IDs and description.

            Parameters
            ----------
            fpath: a file path of 'station.in' to read
        """
        stations = list()
        with open(fpath, 'r') as f:
            # First line
            items = map(int, f.readline().split()[:9])
            self._items = items
            n_stations = int(f.readline().split()[0])
            pattern = r"^(\s*\w+\s+[-\w.]+\s+[-\w.]+\s+[-\w.]+)\s*!?\s*(.*)"
            pattern_annotation = r"(\w+)\s+(.*)"
            for i in range(n_stations):
                line = f.readline()
                matched = re.match(pattern, line)
                if matched is None:
                    print "Error in line: ", line
                    raise ValueError("station.in is not well formatted.")
                tkns = line.split()
                coords = tuple(map(float, tkns[1:4]))
                annotation = matched.group(2).strip()
                name = None
                desc = None
                matched = re.match(pattern_annotation, annotation)
                if matched is not None:
                    name = matched.group(1)
                    desc = matched.group(2)
                stations.append({"coords": coords,
                                 "name": name,
                                 "desc": desc})
        self._stations = stations
        self._map_station_ids()
        self._sort_stations_by_depth()
        self._check_integrity_stations()

    def _map_station_ids(self):
        station_ids = dict()
        for index, station in enumerate(self._stations):
            station_id = station["name"]
            if station_id is not None:
                keys = station_ids.keys()
                if station_id in keys:
                    station_ids[station_id] .append(index)
                else:
                    station_ids[station_id] = [index]
        self._station_id_to_index = station_ids

    def _sort_stations_by_depth(self):
        """ Sort stations with the same id.
        """
        for station_id, indices in self._station_id_to_index.iteritems():
            if len(indices) > 1:
                depths = [(self._stations[i]["coords"][2], i) for i in indices]
                depths = sorted(depths, reverse=True)
                self._station_id_to_index[station_id] = zip(*depths)[1]
                for i, index in enumerate(self._station_id_to_index[station_id]):
                    self._stations[index]["vert_pos"] = i
            else:
                self._stations[indices[0]]["vert_pos"] = 0


    def _check_integrity_stations(self):
        # Check station coords per ID
        for station_id, indices in self._station_id_to_index.iteritems():
            if station_id is not None and len(indices) > 1:
                all_coords = [self._stations[i]["coords"][:2] for i in indices]
                counted = collections.Counter(all_coords)
                if len(counted) != 1:
                    print "Warning: Station does not have a unique horizontal position"
                    print "Station:", station_id, ", indices:", indices

        # Check duplicate coords
        all_coords = [station["coords"] for station in self._stations]
        counted = collections.Counter(all_coords)
        for coords, count in counted.iteritems():
            if count > 1:
                stations = list()
                for index, station in enumerate(self._stations):
                    if station["coords"] == coords:
                        stations.append((index, station["name"]))
                print "Warning: Found stations with identical coordinates in the station file."
                print stations

    def retrieve_ts(self, variable, index=None, name=None, depth=None):
        """ Retrieve a time series of an output constituent for a station
            NOTE: For output with extra vertical outputs.

            Parameters
            ----------
            variable:
                A variable to read in. Currently, it should be one of
                elev, salt.
            index:
                index of the station in station.in. Zero-based.
            name:
                station ID
            depth:
                depth of the station

            Returns
            -------
            a time series
        """
        var = variable.lower()
        if index is None:
            if name is None:
                raise ValueError("Coord or station ID must be provided.")
            if not name in self._station_id_to_index.keys():
                # raise ValueError("No matching station for %s" % name)
                print "Warning: No matching station for %s in station.in" % name
                return None
            indices = self._station_id_to_index[name]
            if len(indices) == 1:
                index = indices[0]
            else:
                if var in ['salt',]: # 3D
                    if depth is None:
                        raise ValueError("Coord or station ID and depth must be provided.")
                    else:
                        index = self._station_id_to_index[name][depth]
                else:  # 2D
                    index = indices[0]
        #     hits = self._index_by_coord(coord)
        # if len(hits) == 0:

        if variable == 'elev':
            if self._elev is None:
                print "Reading %s..." % self._elev_fname
                elev_staout_fname = self._elev_fname
                fname = os.path.join(self._working_dir, elev_staout_fname)
                self._elev = self._read_all_staout(fname)
                # Sanity check
                if len(self._elev) != len(self._stations):
                    return ValueError("# of stations in station.in and staout do not correspond.")
                for ts in self._elev:
                    ts.props['timestamp'] = 'INST-VAL'
                    ts.props['aggregation'] = 'INST-VAL'
                    ts.props['unit'] = 'm'
            return self._elev[index]
        elif variable == 'salt':
            if self._salt is None:
                print "Reading %s..." % self._salt_fname
                salt_staout_fname = self._salt_fname
                fname = os.path.join(self._working_dir, salt_staout_fname)
                self._salt = self._read_all_staout(fname)
                for ts in self._salt:
                    ts.props['timestamp'] = 'INST-VAL'
                    ts.props['aggregation'] = 'INST-VAL'
                    ts.props['unit'] = 'PSU'
            return self._salt[index]
        #         salt_staout_fname = 'staout_6_2'
        #         fname = os.path.join(self._working_dir, salt_staout_fname)
        #         if os.path.exists(fname):
        #             self._salt = self._read_all_staout(fname)
        #         else:
        #             sys.stderr.write("Processed salt output ,staout_6_2, cannot be found.\n")
        #             raise ValueError()
        #     found = self._find_a_station(self._stations, station, depth)
        #     if found is None:
        #         sys.stderr.write("Cannot find a station, %s, %.2f\n" % (station, depth))
        #         return None
        #     else:
        #         index = self._stations.index(found)
        #         return self._salt[index]
        else:
            sys.stderr.write("Not one of supported variables: %s\n" % variable)
            return None

    def _read_all_staout(self, fname):
        """ Read the whole staout into the memory.
            Each column will be converted into vtools.data.timeseries.TimeSeries.
        """
        raw = np.loadtxt(fname)
        # Get the first time stamp
        ts_begin = datetime.timedelta(seconds=raw[0, 0])
        # Get dt
        for i in xrange(raw.shape[0]):
            dt = raw[i+1, 0] - raw[i, 0]
            if dt > 0.:
                if i > 0:
                    raw = raw[::i+1]
                dt = datetime.timedelta(seconds=dt)
                break
        # Remove dry spots
        raw[raw < -100.] = np.nan
        # Convert into Vtools time series
        tss = list()
        for i in range(1, raw.shape[1]):
            ts = rts(raw[:, i], self._time_basis + ts_begin, dt)
            tss.append(ts)
        return tss

    # @property
    # def stations(self):
    #     return self._stations

    # @stations.setter
    # def stations(self, value):
    #     self._stations = value


def station_extractor(station_file, working_dir, time_basis):
    """ Create a staout, SCHISM standard output file, extractor

        Parameters
        ----------
        station_file:
            A staout file path to read
        working_dir:
            A directory where staout file resides
        time_basis:
            A time basis of the staout file

        Returns
        -------
        An instance of the station extractor
    """
    sr = StaoutReader(working_dir, time_basis)
    sr.read_station_in(station_file)
    return sr


class FlowReader():
    """ SCHISM flux.dat reader
    """

    # class Station(object):
    #     def __init__(self, name, coord):
    #         self._name = name
    #         self._coord = coord

    #     def __str__(self):
    #         return "%s, (%f, %f ,%f, %f)" % (self._name,
    #             self._coord[0], self._coord[1],
    #             self._coord[2], self._coord[3])

    #     @property
    #     def name(self):
    #         return self._name
    #     @name.setter
    #     def name(self, value):
    #         self._name = value

    #     @property
    #     def coord(self):
    #         return self._coord

    #     @coord.setter
    #     def coord(self, value):
    #         self._coord = value

    def __init__(self, working_dir, time_basis):
        """ Constructor

            Parameters
            ----------
            working_dir: A full path of flux.dat
            time_basis: A time basis of the flux.dat
        """
        self._working_dir = working_dir
        self._time_basis = time_basis
        self._items = None
        self._flow = None
        self._flow_fname = 'flux.dat'
        self._stations = None
        self._station_id_to_index = None

    @property
    def flow_fname(self):
        return self._flow_fname
    @flow_fname.setter
    def flow_fname(self, value):
        self._flow_fname = value

    @property
    def stations(self):
        return self._stations
    # @stations.setter
    # def stations(self, value):
    #     self._stations = value

    # def __iter__(self):
    #     for st in zip(self._names, self._coords):
    #         yield self.Station(*st)

    def read_flow_input(self, fpath):
        """ Read in 'flowlines.yaml'

            Parameters
            ----------
            fpath: a file path of the input
        """
        yaml.add_representer(collections.OrderedDict, dict_representer)
        yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, dict_constructor)
        with open(fpath, 'r') as f:
            flowlines = yaml.load(f)
            stations = list()
            for name, value in flowlines.iteritems():
                stations.append({"name": str(name),
                                 "coord": tuple(map(float, value.split()))})
            self._stations = stations
            self._map_station_ids()
    # def _index_by_name(self, name):
    #     return [i for i, j in enumerate(self._names) if j == name]

    # def _index_by_coord(self, coord):
    #     return [i for i, j in enumerate(self._coords) if j[:len(coord)] == coord]

    def retrieve_ts(self, name=None):
        """ Retrieve a time series of an output constituent for a station

            Parameters
            ----------
            name: str
                station ID in flowlines.yaml

            return
            ------
            vtools.data.timeseries.TimeSeries
                a flow time series
        """
        if not name in self._station_id_to_index.keys():
            # raise ValueError("No matching station for %s" % name)
            print "Warning: No matching station for %s in station.in" % name
            return None
        indices = self._station_id_to_index[name]

        if len(indices) < 1:
            print "More than two outputs with the same name."
            print "Use the first one."

        if self._flow is None:
            print "Reading flow...: ", self.flow_fname
            fname = os.path.join(self._working_dir, self.flow_fname)
            self._flow = self._read_all_flow(fname)
        return self._flow[indices[0]]

    def _read_all_flow(self, fname):
        raw = np.loadtxt(fname)
        # Get the first time stamp
        ts_begin = datetime.timedelta(seconds=np.around(raw[0, 0]*86400))
        # Get dt
        for i in xrange(raw.shape[0]):
            dt = np.around((raw[i+1, 0] - raw[i, 0])*86400.)
            if dt > 0.:
                if i > 0:
                    raw = raw[::i+1]
                dt = datetime.timedelta(seconds=dt)
                break
        # Convert into Vtools time series
        tss = list()
        for i in range(1, raw.shape[1]):
            ts = rts(raw[:, i], self._time_basis + ts_begin, dt)
            ts.props['aggregation'] = 'INST-VAL'
            ts.props['timestamp'] = 'INST-VAL'
            ts.props['unit'] = 'cms'
            tss.append(ts)
        return tss

    def _map_station_ids(self):
        station_ids = dict()
        for index, station in enumerate(self._stations):
            station_id = station["name"]
            if station_id is not None:
                keys = station_ids.keys()
                if station_id in keys:
                    station_ids[station_id] .append(index)
                else:
                    station_ids[station_id] = [index]
        self._station_id_to_index = station_ids


def flow_extractor(flow_input_file, working_dir, time_basis):
    """ Create a flow extractor

        Return
        ------
        FlowReader object
    """
    fr = FlowReader(working_dir, time_basis)
    fr.read_flow_input(flow_input_file)
    return fr
