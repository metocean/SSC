""" A package to read a mesh in GR3 format.
"""
## Author: Kijin Nam, knam@water.ca.gov

from schism_mesh import *
from trimesh import *
from base_io import *
import numpy as np
import os

class Gr3IO(BaseIO):
    """ A class that manages I/O of GR3 files
    """

    def __init__(self):
        """ Constructor
        """
        super(Gr3IO, self).__init__()
        self._mesh = None

    def read(self, gr3_fname="hgrid.gr3", mode=0):
        """ Read in a hgrid.gr3 file.
            If mode is 1, it does not read in boundary information.
        """
        if self._verbose > 0:
            print "Reading in a gr3 file...", gr3_fname
        # Create a new mesh
        self._mesh = SchismMesh()
        try:
            f = open(gr3_fname)
        except:
            raise Exception("A grid file is not found.")

        self._linecounter = 0
        # Mesh info
        self._read_header(f)
        if self._verbose > 1:
            print "Reading in nodes..."
        self._read_nodes(f)
        if self._verbose > 1:
            print "Reading in elements..."
        self._read_elems(f)

        # Build grid
        if self._verbose > 1:
            print "Building a mesh..."
        self._mesh.build_edges_from_elems()

        # Boundary cinfo
        if not mode == 1:
            if self._verbose > 1:
                print "Reading in boundaries..."
            self._read_boundaries(f)

        ## Close
        f.close()
        if self._verbose > 0:
            print "Done reading a GR3 file."
        return self._mesh

    def _read_header(self, f):
        ## Header
        # First line: ignored
        tokens, ok = self._read_and_parse_line(f, 0)
        try:
	     if len(tokens) > 0:
                 self._mesh.comment = tokens[0]
	except:
		self._mesh.comment=None
        # Second line: # of elements and # of nodes
        tokens, ok = self._read_and_parse_line(f, 2)
        if not ok:
            raise Exception("The Header information of GR3 file is corrupted.")
        (n_elems, n_nodes) = map(int, tokens[:2])
        if n_nodes <= 0 or n_elems <= 0:
            raise Exception("The Header information of GR3 file is corrupted.")
        self._mesh.allocate(n_elems, n_nodes) # Allocate memory
        self._n_elems = n_elems
        self._n_nodes = n_nodes

    def _read_nodes(self, f):
        node_counter = 0
        for i in range(self._n_nodes):
            tokens, ok = self._read_and_parse_line(f, 4)
            if not ok:
                print "Line: ", self.linecounter
                raise Exception("Node block is corrupt.")
            node_coords = map(float, tokens[1:4])
            self._mesh.set_node(node_counter, node_coords)
            node_counter += 1

    def _read_elems(self, f):
        for elem_i in range(self._n_elems):
            tokens, ok = self._read_and_parse_line(f, 5)
            if not ok:
                print "Line: ", self.linecounter
                raise Exception("Element block is corrupt.")
            type_elem = int(tokens[1])
            if not type_elem == 3:
                raise Exception("Currently only triangular " \
                                "are supported.")
            # Zero-based connectivities
            connectivities = np.subtract(np.array(map(int, tokens[2:5])), 1)
            self._mesh.set_elem(elem_i, connectivities)

    def _read_boundaries(self, f):
        """ Read boundary information
        """
        ## Assuming the order of the different boundary types are consitent
        ## Open boundaries
        # # of open boundaries
        tokens, ok = self._read_and_parse_line(f, 1)
        if not ok:
            print "No boundary information is present?"
            return

        n_open_boundaries = int(tokens[0])
        # total # of open boundary nodes
        tokens, ok = self._read_and_parse_line(f, 1)
        if not ok:
            print "Line: ", self.linecounter
            raise Exception("The number of total open boundary nodes is not"\
                            " correctly provided.")
        n_open_boundary_nodes = int(tokens[0])

        # open boundaries
        for i_boundary in range(n_open_boundaries):
            # # of nodes of this open boundary
            tokens, ok = self._read_and_parse_line(f, 1)
            if not ok:
                print "Line: ", self.linecounter
                raise Exception("The number of nodes for a boundary is not"\
                                " correctly provided.")
            # Comment
            comment = None
            if len(tokens) > 1:
                comment = tokens[1]

            n_nodes = int(tokens[0])
            nodes = []
            for i_node in range(n_nodes):
                tokens, ok = self._read_and_parse_line(f, 1)
                if not ok:
                    print "Line: ", self.linecounter
                    raise Exception("A node for a boundary is not"\
                                    " correctly provided.")

                node = int(tokens[0]) - 1 # Zero based
                nodes.append(node)
            self._mesh.add_boundary(nodes, OPEN_BOUNDARY, comment)

        # land boundaries
        # I found out that there is no distinction between
        # land and island boundaries.
        tokens, ok = self._read_and_parse_line(f, 1)
        if not ok:
            print "No land boundary presented?"
            return
        n_land_boundaries = int(tokens[0])
        # total # of land boundary nodes
        tokens, ok = self._read_and_parse_line(f, 1)
        if not ok:
            print "Line: ", self.linecounter
            raise Exception("The number of total land boundary nodes is not"\
                            " correctly provided.")
        n_land_boundary_nodes = int(tokens[0])
        #import pdb;pdb.set_trace()
	for i_boundary in range(n_land_boundaries):
	        
            # # of nodes of this open boundary
            (tokens, ok) = self._read_and_parse_line(f,1)
            if not ok:
                print "Line: ", self.linecounter
                raise Exception("The number of nodes for a boundary is not"\
                                " correctly provided.")
            # Comment
            comment = None
            #if len(tokens) > 1:
            #    comment = tokens[1]

            n_nodes = int(tokens[0])

#            print i_boundary,tokens[1].partition(' ')[0]
#	    btype = int(tokens[1].partition(' ')[0])  # Obsolete
            stre=tokens[1].partition(' ')[0]
            stre=stre.replace('\t','')
	    stre=stre.replace('=','')
	    btype = int(stre)  # Obsolete
     
	    ### take care of btype if adcirc
	    if btype ==10 or btype==11:
		btype=btype-10
		print "This is ADCIRC BTYPE==> double check"
	    
	    comment= str(btype) + '= Land/Island boundary ' +str(i_boundary+1)
            nodes = []
            for i_node in range(n_nodes):
                (tokens, ok) = self._read_and_parse_line(f, 1)
                if not ok:
                    print "Line: ", self.linecounter
                    raise Exception("A node for a boundary is not correctly provided.")
                node = int(tokens[0]) - 1 # Zero based
                nodes.append(node)
            #self._mesh.add_boundary(nodes, LAND_BOUNDARY, comment)

 
	    self._mesh.add_boundary(nodes, btype, comment)
	    
    def write(self, mesh, fname, node_attr = None, boundary = False):
        """ Write a GR3 format grid.
            mesh = SCHISM mesh (schism_mesh) instance
            fname = output file name
            node_attr = a list of node attribute
            boundary = If true, boundary information will be added. Otherwise,
            it will not be appended.
        """
        print "Writing an hgrid file:", fname
        f = open(fname, 'w')
        # Header
#         if mesh.comment is None:
        buf = "%s \n" % os.path.basename(fname)
#         else:
#             buf = "%s !modified by the preprocessing tool\n" \
#                   % mesh.comment

        f.write(buf)
        n_elems = mesh.n_elems()
        n_nodes = mesh.n_nodes()
        buf = "%d %d ! # of elements and nodes \n" % (n_elems, n_nodes)
        f.write(buf)
	
        # Nodes
        for i in range(n_nodes):
            if not node_attr is None:
                buf = "%d %.7f %.7f %.7g\n" % (i + 1, \
                                     mesh.nodes[i, 0], \
                                     mesh.nodes[i, 1], \
                                     node_attr[i])

            else:
                buf = "%d %.7f %.7f %.3f\n" % (i + 1, \
                                     mesh.nodes[i, 0], \
                                     mesh.nodes[i, 1], \
                                     mesh.nodes[i, 2])
            f.write(buf)

        # Elements
        for i in range(n_elems):
            buf = "%d 3 %d %d %d\n" % (i + 1, \
                                       mesh.elems[i, 0] + 1, \
                                       mesh.elems[i, 1] + 1, \
                                       mesh.elems[i, 2] + 1)
            f.write(buf)

        # Boundaries
        if boundary:
            # Open
            buf = "%d = Number of open boundaries\n" \
                  % mesh.n_boundaries(OPEN_BOUNDARY)
            f.write(buf)
            buf = "%d = Total number of open boundary nodes\n" \
                  % mesh.n_boundary_nodes(OPEN_BOUNDARY)
            f.write(buf)
            openbound_count = 0
            for boundary in mesh.boundaries:
                if boundary.btype == OPEN_BOUNDARY:
                    openbound_count += 1
                    if boundary.comment is None:
                        buf = "%d = Number of nodes for open boundary %d\n" % \
                              (boundary.n_nodes(), openbound_count)
                    else:
                        buf = "%d %s\n" % (boundary.n_nodes(), boundary.comment)
                    f.write(buf)
                    buf = ""
                    for node_i in boundary.nodes:
                        buf += "%d\n" % (node_i + 1)
                    f.write(buf)
#             else:
#                 raise Exception("Unsupported boundary type.")

            # Land
            
	   
            buf = "%d = Number of land boundaries\n" \
                  % (mesh.n_boundaries(LAND_BOUNDARY) + \
                     mesh.n_boundaries(ISLAND_BOUNDARY))
            f.write(buf)
            buf = "%d = Total number of land boundary nodes\n" \
                  % (mesh.n_boundary_nodes(LAND_BOUNDARY) + \
                     mesh.n_boundary_nodes(ISLAND_BOUNDARY))
            f.write(buf)
            landbound_count = 0
            islandbound_count = 0
            for boundary in mesh.boundaries:
                if boundary.btype == LAND_BOUNDARY:
                    landbound_count += 1
                    if boundary.comment is None:
                        buf = "%d = Number of nodes for land boundary %d\n" % \
                          (boundary.n_nodes(), landbound_count)
                    else:
                        buf = "%d %s\n" % (boundary.n_nodes(), boundary.comment)
                    f.write(buf)
                    buf = ""
                    for node_i in boundary.nodes:
                        buf += "%d\n" % (node_i + 1)
                    f.write(buf)
                elif boundary.btype == ISLAND_BOUNDARY:
                    islandbound_count += 1
                    if boundary.comment is None:
                        buf = "%d = Number of nodes for island boundary %d\n" \
                               % (boundary.n_nodes(), islandbound_count)
                    else:
                        buf = "%d %s\n" % (boundary.n_nodes(), boundary.comment)
                    f.write(buf)
                    buf = ""
                    for node_i in boundary.nodes:
                        buf += "%d\n" % (node_i + 1)
                    f.write(buf)
#             else:
#                 raise Exception("Unsupported boundary type.")

        f.flush()
        f.close()

        print "Done writing a GR3 file."
