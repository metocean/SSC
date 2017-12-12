#!/usr/bin/env python2.7
""" Driver module to prepares input files for a SCHISM run.
    Example jobs in this pre-processing are:
      1. Fill up missing land/island boundaries in `hgrid.gr3'.
      2. Re-order the open boundaries.
      3. Create source_sink.in and hydraulics.in
      4. Create spatial information such as elev.ic, rough.gr3 and xlsc.gr3.
      5. Dredge up some part of the domain.
      6. Create an ocean boundary file, elev2D.th
"""

# Author: Kijin Nam, knam@water.ca.gov

import grid_opt
import sms2gr3
import stacked_dem_fill
from schism_setup import *
import yaml
import numpy as np
import collections
import string
import subprocess
from datetime import *
import shutil
import os
import argparse

# Ordered Dict YAML
# From http://stackoverflow.com/questions/5121931/in-python-how-can-you-load-yaml-mappings-as-ordereddicts
def dict_representer(dumper, data):
    return dumper.represent_mapping(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, data.iteritems())


def dict_constructor(loader, node):
    return collections.OrderedDict(loader.construct_pairs(node))


def create_arg_parser():
    # Read in the input file
    parser = argparse.ArgumentParser(description='Prepare SCHISM input files.')
    parser.add_argument(dest='main_inputfile', default=None,
                        help='main input file name')
    return parser


def get_land_boundary(mesh):
    """ Create an array of land boundary information that conforms
        RueenFang's data structure

        Parameters
        ----------
        mesh: schism_mesh


        Returns
        -------
        numpy array
            land baoundaries, [ind, node1, node2, ind_section, ind_in_section]
    """
    new_bound = []
    ind = 0
    ind_bound = 0
    for boundary in mesh.boundaries:
        if boundary.btype >= trimesh.LAND_BOUNDARY:
            nodes = boundary.nodes
            for i in range(len(nodes) - 1):
                item = [ind, nodes[i], nodes[i+1], ind_bound, i]
                new_bound.append(item)
                ind += 1
            ind_bound += 1
    return np.array(new_bound)

def create_hgrid(s, inputs):
    # Preprocessed the hgrid file
    section_name = 'mesh'
    if section_name in inputs.keys():
        section = inputs[section_name]
        # First, boundary building
        option_name = 'open boundary file'
        if option_name in section.keys():
            openboundary_fpath = section[option_name]
            with open(openboundary_fpath, 'r') as f:
                open_boundary_segments = yaml.load(f)
            s.create_open_boundaries(open_boundary_segments)

            # Fill the missing land and island boundary information
            print "Filling missing land and island boundaries..."
            s.mesh.fill_land_and_island_boundarise()
        # Second, Mesh optimization
        option_name = 'depth optimization file'
        if option_name in section.keys():
            opt_param_fpath = os.path.expanduser(section[option_name])
            dem_list_fpath = os.path.expanduser(section['dem list file'])
            with open(opt_param_fpath, 'r') as f:
                opt_params = yaml.load(f)
                opt_params = [opt_params[key] for key in ["damp", "damp_shoreline", "face_coeff", "volume_coeff"]]
                land_boundaries = get_land_boundary(s.mesh)
                new_z = grid_opt.grid_opt(s.mesh.nodes, s.mesh.elems,
                                          land_boundaries,
                                          dem_list_fpath, [opt_params])
                s.mesh.nodes[:,2] = np.negative(new_z[0])

        # Write hgrid.gr3
        option_name = 'gr3 output file'
        if option_name in section.keys():
            print "Write up a new hgrid file..."
            hgrid_out_fpath = os.path.expanduser(section[option_name])
            s.write_hgrid(hgrid_out_fpath)

        # Write hgrid.ll
        option_name = 'll output file'
        if option_name in section.keys():
            print "Create a new hgrid.ll file..."
            hgrid_ll_fpath = os.path.expanduser(section[option_name])
            s.write_hgrid_ll(hgrid_ll_fpath)


def create_source_sink(s, inputs):
    # Create source_sink.in
    section_name = 'source/sink'
    if section_name in inputs.keys():
        section = inputs[section_name]
        option_name = 'input file'
        if option_name in section.keys():
            source_sink_in_fpath = os.path.expanduser(section[option_name])
            with open(source_sink_in_fpath, 'r') as f:
                source_sink = yaml.load(f)
            option_name = 'output file'
            if option_name in section.keys():
                source_sink_out_fpath = os.path.expanduser(section[option_name])
                print "Creating %s..." % source_sink_out_fpath
                s.create_source_sink_in(source_sink, source_sink_out_fpath)


def create_gr3_with_constant(s, inputs):
    # Create GR3 files with one values
    section_name = 'gr3 with constant'
    if section_name in inputs.keys():
        section = inputs[section_name]
        for k, v in section.iteritems():
            gr3_fpath = os.path.expanduser(k)
            # TODO: Improve it
            attr_array = np.empty(s.mesh.n_nodes())
            attr_array.fill(float(v))
            print "Creating %s..." % gr3_fpath
            s.write_hgrid(gr3_fpath, attr_array, False)


def create_gr3_with_polygons(s, inputs):
    # Create GR3 files with polygons
    section_name = 'gr3 with polygons'
    if section_name in inputs.keys():
        section = inputs[section_name]
        for k, v in section.iteritems():
            gr3_fpath = os.path.expanduser(k)
            polygon_fpath = os.path.expanduser(v)
            with open(polygon_fpath, 'r') as f:
                polygons = yaml.load(f)
            print "Creating %s..." % gr3_fpath
            s.create_node_partitioning(gr3_fpath, polygons)


def create_prop_with_polygons(s, inputs):
    # Create prop files with polygons
    section_name = 'prop with polygons'
    if section_name in inputs.keys():
        section = inputs[section_name]
        for k, v in section.iteritems():
            prop_fpath = os.path.expanduser(k)
            polygon_fpath = os.path.expanduser(v)
            with open(polygon_fpath, 'r') as f:
                polygons = yaml.load(f)
            print "Creating %s..." % prop_fpath
            s.create_prop_partitioning(prop_fpath, polygons)


def create_structures(s, inputs):
    # Create a structure file
    section_name = 'hydraulics'
    if section_name in inputs.keys():
        section = inputs[section_name]
        option_name = 'input file'
        if option_name in section.keys():
            structure_in_fpath = os.path.expanduser(section[option_name])
            print "Reading %s..." % structure_in_fpath
            with open(structure_in_fpath, 'r') as f:
                structures = yaml.load(f)
            s.create_structures(structures)
            option_name = 'output file'
            if option_name in section.keys():
                structure_out_fpath = os.path.expanduser(section[option_name])
                print "Creating %s..." % structure_out_fpath
                s.write_structures(structure_out_fpath)


def create_fluxflag(s, inputs):
    # Create fluxflag.gr3
    section_name = 'flow output'
    if section_name in inputs.keys():
        section = inputs[section_name]
        option_name = 'input file'
        if option_name in section.keys():
            flux_in_fpath = os.path.expanduser(section[option_name])
            with open(flux_in_fpath, 'r') as f:
                fluxlines = yaml.load(f)
            option_name = 'output file'
            if option_name in section.keys():
                flux_out_fpath = os.path.expanduser(section[option_name])
                print "Creating %s..." % flux_out_fpath
                s.create_flux_regions(fluxlines, flux_out_fpath)


def update_spatial_inputs(s, inputs):
    """ Create SCHISM grid inputs
        s = schism setup object
        inputs = inputs from an input file
    """
    create_hgrid(s, inputs)
    create_source_sink(s, inputs)
    create_gr3_with_constant(s, inputs)
    create_gr3_with_polygons(s, inputs)
    create_prop_with_polygons(s, inputs)
    create_structures(s, inputs)
    create_fluxflag(s, inputs)


def update_temporal_inputs(s, inputs):
    """ Create temporal inputs. Under development
    """
    # create in interpolated tide file
    sf_tide_out_fpath = os.path.join(output_dir, sf_tide_out_fname)
    s.interpolate_tide(time_start, time_end, dt,
                       sf_tide_in_fpath, sf_tide_out_fpath)
    # Run the FORTRAN code to create elev2D.th
    hgrid_out_fpath = os.path.join(output_dir, hgrid_out_fname)
    webtide_grid_fpath = os.path.join(input_dir, webtide_grid_fname)
    webtide_fpath = os.path.join(input_dir, webtide_fname)
    elev2d_fpath = os.path.join(output_dir, elev2d_fname)
    p = subprocess.Popen(["./gen_elev2D_4_NAVD88", sf_tide_out_fpath,
                      hgrid_out_fpath, webtide_grid_fpath, webtide_fpath,
                      elev2d_fpath],
                     stdout = subprocess.PIPE,
                     stderr = subprocess.PIPE)
    return_code = p.wait()
    if return_code != 0:
        for l in p.stdout:
            print l
        for l in p.stderr:
            print l


def substitute(inputs, env=None):
    if env is None:
        if 'env' in inputs.keys():
            env = inputs['env']
    if env is not None:
        for k, v in inputs.iteritems():
            if isinstance(v, dict):
                v = substitute(v, env)
            else:
                if v is not None:
                    if isinstance(v, str) and '$' in v:
                        temp = string.Template(v)
                        new_v = temp.substitute(**env)
                        inputs[k] = new_v
            if '$' in k:
                temp = string.Template(k)
                new_k = temp.substitute(**env)
                inputs[new_k] = v
                del inputs[k]
        return inputs

def item_exist(inputs, name):
    return True if name in inputs.keys() else False


def main():
    parser = create_arg_parser()
    args = parser.parse_args()
    in_fname = args.main_inputfile

    yaml.add_representer(collections.OrderedDict, dict_representer)
    yaml.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
                         dict_constructor)

    with open(in_fname ,'r') as f:
        inputs = yaml.load(f)
    substitute(inputs)

    out_fname = os.path.splitext(in_fname)[0] \
        + '_echo' + os.path.splitext(in_fname)[1]
    with open(out_fname, 'w') as f:
        f.write(yaml.dump(inputs, default_flow_style=False))

    if item_exist(inputs, 'mesh'):
        if item_exist(inputs['mesh'], 'mesh input file'):
            # Read the grid file to be processed
            mesh_input_fpath = os.path.expanduser(inputs['mesh']['mesh input file'])
            mesh_input_froot, ext = os.path.splitext(mesh_input_fpath)
            if ext == '.2dm':
                # Convert 2dm to gr3 first
                if item_exist(inputs['mesh'], 'dem list file') and \
                   not item_exist(inputs['mesh'], 'depth optimization file'):
                    demlist_fname = inputs['mesh']['dem list file']
                    print "Populate elevation with DEM file...", demlist_fname
                    bak_fname = mesh_input_fpath + '.bak'
                    print "Backing up the 2dm...", mesh_input_fpath
                    shutil.copyfile(mesh_input_fpath, bak_fname)
                    stacked_dem_fill.fill_2dm(bak_fname,
                                              mesh_input_fpath, demlist_fname)
                print "Convert 2dm to gr3..."
                gr3_input_fpath = mesh_input_froot + '.gr3'
                sms2gr3.convert_2dm(mesh_input_fpath, gr3_input_fpath, True)
            else:
                if ext != '.gr3':
                    print "Assume the mesh file is in GR3 format."
                if item_exist(inputs['mesh'], 'dem list file') and \
                   not item_exist(inputs['mesh'], 'depth optimization file'):
                    demlist_fname = inputs['mesh']['dem list file']
                    print "Populate elevation with DEM file...", demlist_fname
                    bak_fname = mesh_input_fpath + '.bak'
                    print "Backing up the gr3...", mesh_input_fpath
                    shutil.copyfile(mesh_input_fpath, bak_fname)
                    stacked_dem_fill.fill_gr3(bak_fname,
                                              mesh_input_fpath, demlist_fname, True)
                gr3_input_fpath = mesh_input_fpath
            s = load_gr3(gr3_input_fpath)
        else:
            raise ValueError("No mesh input file in the mesh section.")
    else:
        raise ValueError("No mesh section in the main input.")
    update_spatial_inputs(s, inputs)
    print "Done."

if __name__ == "__main__":
    main()
