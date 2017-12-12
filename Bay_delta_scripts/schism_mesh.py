""" schism_mesh module
    SchismMesh class holds mesh and boundary information
"""
## Author: Kijin Nam, knam@water.ca.gov

import trimesh
import numpy as np
import sets

class SchismBoundary(object):
    """ A class describing an individual boundary.
    """
    def __init__(self, nodes, btype, comment = None):
        self._nodes = nodes
        self._btype = btype
        self._comment = comment

    @property
    def btype(self):
        """ Boundary type.

            :getter: Get the type of the boundary.
        """
        return self._btype

    @property
    def nodes(self):
        """ A sequence of node indices in the boundary.

            :getter: Get the sequence of nodes in the boundary.
        """
        return self._nodes

    def n_nodes(self):
        """ Get the total number of nodes in the boundary.

            Return
            ------
            integer
                the total number of nodes in the boundary.
        """
        return len(self._nodes)

    @property
    def comment(self):
        """ A comment or name of the boundary.

            :getter: Get the comment of the boundary.
            :setter: Set the comment of the boundary.
        """
        return self._comment

    @comment.setter
    def comment(self, value):
        self._comment = value


class SchismMesh(trimesh.TriMesh):
    """ Mesh geometry, topology, boundaries and hydraulic structures.
    """
    def __init__(self):
        super(SchismMesh, self).__init__()
    #   self._boundaries = SchismBoundaries.SchismBoundaries()
        self._boundaries = []
        self._comment = None

    @property
    def comment(self):
        """ Comment of the mesh.

            :getter: Get the comment
            :setter: Set the comment
        """
        return self._comment

    @comment.setter
    def comment(self, value):
        self._comment = value

    @property
    def boundaries(self):
        """ An array of the boundary information

            :getter: Get the array of the boundary information
        """
        return self._boundaries

    def add_boundary(self, nodes, btype, comment=None):
        """ Add one boundary.

            Parameter
            ---------
            nodes: an array of integer
                array of inter node indices in the boundary path

            btype: integer
                an integer constant for boundary types.

            comment: optional, string
        """
        nodes = self._rearrange_boundary_nodes_in_ccw(nodes)
        # Call the super class method first
        super(SchismMesh, self).add_boundary(nodes, btype)

        self._boundaries.append(SchismBoundary(nodes, btype, comment))

    def _rearrange_boundary_nodes_in_ccw(self, nodes):
        """ Make sure the boundary nodes are in CCW ordering.
            Assume that nodes are in row at least.
            nodes = the list of node indices
            return = reordered list of the node indices
        """
        # TODO: I may need to consider nodes not in any order...
        edge_i = self._find_edge(nodes[:2], True)
        if edge_i is None:
            new_nodes = []
            for node in reversed(nodes):
                new_nodes.append(node)
            del nodes
            return new_nodes
        else:
            return nodes

    def n_boundaries(self, btype = None):
        """ Get the number of open boundaries. If a boundary type is given,
            it counts only open boundaries with the corresponding type.

            Parameter
            ---------
            btype: integer, optional
                Type of the boundary.
        """
        if btype is None:
            return len(self._boundaries)
        else:
            return sum(1 for b in self._boundaries \
                   if b.btype == btype)

    def n_boundary_nodes(self, btype):
        """ Get the total node boundary of a given type

            Parameter
            ---------
            btype: integer
                Type of the boundary.
        """
        return sum(b.n_nodes() for b in self._boundaries \
               if b.btype == btype)

    def clear_boundary(self):
        """ Delete all the current boundary information
        """
        del self._boundaries[:]
        self._clear_edge_types()


    def _find_intersection(self, s1, s2):
        """ Calculate an intersection point from two line segments
            defined by two end points
        """
        det = (s1[0,0] - s1[1,0]) * (s2[0,1] - s2[1,1]) - \
              (s1[0,1] - s1[1,1]) * (s2[0,0] - s2[1,0])
        if det == 0.:  # parallel
            return None
        x = (s1[0,0] * s1[1,1] - s1[0,1] * s1[1,0]) * \
            (s2[0,0] - s2[1,0]) - \
            (s1[0,0] - s1[1,0]) * \
            (s2[0,0] * s2[1,1] - s2[0,1] * s2[1,0])
        y = (s1[0,0] * s1[1,1] - s1[0,1] * s1[1,0]) * \
            (s2[0,1] - s2[1,1]) - \
            (s1[0,1] - s1[1,1]) * \
            (s2[0,0] * s2[1,1] - s2[0,1] * s2[1,0])
        intersection =  np.array([x, y]) / det
        v1 = s1[0,:] - intersection
        v2 = s1[1,:] - intersection
        sign1 = np.dot(v1, v2)
        v1 = s2[0,:] - intersection
        v2 = s2[1,:] - intersection
        sign2 = np.dot(v1, v2)
        if sign1 > 0. or sign2 > 0.:
            return None
        else:
            return intersection

    def find_two_neighboring_paths(self, line_segment):
        """ Format of the line segment: start_x, start_y, end_x, end_y
            line_segment = a pair of (X, Y)

            Parameter
            ---------
            line_segment: list or array of py:.float with size of 4



            return
            ------
            up_path: list of integers
                upstream side of node paths
            down_path: list of integers
                downstream side of node paths
        """
        # TODO: This routine may return some extra if the line segment
        # does not run through from outside to outside of the domain.
        if len(line_segment) != 4:
            raise ValueError("The size of line_segment array is incorrect."
                            " It must be four.")

        if self._elem_index is None:
            self._build_elem_index()

        x = np.array(line_segment)
        x = x.reshape(2, 2)
        normal = np.array((x[1, 1] - x[0, 1], x[0, 0] - x[1, 0]))
        box = self._box_from_points(x)
        hits = self._elem_index.intersection(box)
        # Test which one is actually intersect
        up_path = sets.Set()
        down_path = sets.Set()
        for hit in hits:
            nodes_i = self._elems[hit]
            nodes = self._nodes[nodes_i][:, :2]
            intersections = []
            for i in range(3):
                edges = np.roll(nodes, -i, axis=0)[:2]
                intersection = self._find_intersection(x, edges)
                intersections.append(intersection)
            for i, intersection in enumerate(intersections):
                if intersection is not None:
                    if np.inner(normal, nodes[i] - intersection) > 0.:
                        up_path.add(nodes_i[i])
                        down_path.add(nodes_i[(i+1)%3])
                    else:
                        down_path.add(nodes_i[i])
                        up_path.add(nodes_i[(i+1)%3])
        up_path = list(up_path)
        down_path = list(down_path)

        # Let's order them up
        up_path = self._order_up_nodes_from_point(up_path, x)
        down_path = self._order_up_nodes_from_point(down_path, x)

        return up_path, down_path

    def _order_up_nodes_from_point(self, nodes, x):
        """ Order up nodes based on the distance from x
        """
        dist = self._distance(nodes, x)
        sorted_indices = np.argsort(dist)
        sorted = list(nodes[i] for i in sorted_indices)
        return sorted

    def _distance(self, nodes_i, x):
        nodes = tuple(self._nodes[i] for i in nodes_i)
        diffs = tuple(np.subtract(nodes[i][:2], x[0,]) \
                     for i in range(len(nodes)))
        dist = tuple(np.linalg.norm(diff) for diff in diffs)
        return dist

    def _check_integrity_of_flowline(self, up_path, down_path):
        """ Check the the integrtity of the flowline which consists of
            a pair of paths.
        """
        if not self._check_if_path_is_continuous(up_path):
            raise Exception('Two nodes are not neighboring')


    def _check_if_path_is_continuous(self, path):
        """ Check if a node string is well constructed continuously
        """
        node_prev = None
        for node_i in path:
            if node_prev is None:
                node_prev = node_i
            else:
                edge_i = self._find_edge([node_prev, node_i])
                if edge_i is None:
                    return False
                node_prev = node_i
        return True

    def create_boundaries(self, open_pairs):
        """ Create boundaries from the given pairs of coordinates for
            open boundaries
        """
        self.clear_boundary()  # Clear the current boundary information first
        # Open boundary first
        for open in open_pairs:
            ns = self._create_boundary(open)
            self.add_boundary(ns, OPEN_BOUNDARY)

        self.fill_land_and_island_boundarise()

    def fill_land_and_island_boundarise(self):
        """ Fill land and island boundaries for boundary edges not assigned
            to any open boundary.
        """
        # Now build land boundaries.
        # Boundary edges not assigned for open boundaries ought to be
        # land boundaries.
        self._fill_land_boundaries()
        print "Land boundaries filled..."
        # Check if there are still boundary edges not assigned for any
        # type of boundary.  If so, assume those are island boundaries.
        not_assigned = self._get_not_assigned_boundary_edges()
        done = False
        i = 0
        while not done:
            i += 1
            if len(not_assigned) > 0:
                self._fill_island_boundaries(not_assigned)
            else:
                done = True

    def _create_boundary(self, pair):
        """ Create a node string for a boundary from the given coordinate pair.
        """
        b = None
        nodes = []
        n1 = self.find_closest_nodes(pair[0], 1, True)
        n2 = self.find_closest_nodes(pair[1], 1, True)
        ns = self._build_boundary_node_string(n1, n2)
        return ns

    def _fill_island_boundaries(self, not_assigned):
        """ This function fills missing island boundaries.
        """
        ns = []
        first_node_i = self._edges[not_assigned[0]][0]
        ns.append(first_node_i)
        last_node_i = first_node_i
        done = False
        while not done:
            next_node = self._get_next_node_on_boundary_and_remove_edge(
                last_node_i, not_assigned)
            ns.append(next_node)
            if next_node == first_node_i:
                done = True
                self.add_boundary(ns, trimesh.ISLAND_BOUNDARY)
            else:
                last_node_i = next_node


    def _get_next_node_on_boundary_and_remove_edge(self, node_i,
                                                   not_assigned,
                                                   ccw = True):
        edges_i = self.get_edges_from_node(node_i)
        for edge_i in edges_i:
            edge = self._edges[edge_i]
            if ccw:
                if (not edge[2] == trimesh.INTERNAL_EDGE) and edge[0] == node_i:
                    not_assigned.remove(edge_i)
                    return edge[1]
            else:
                if (not edge[2] == trimesh.INTERNAL_EDGE) and edge[1] == node_i:
                    not_assigned.remove(edge_i)
                    return edge[0]
        return None

    def _fill_land_boundaries(self):
        """ This function fills the land boundary after all open boundaries
            are created in CCW direction.
        """
        if len(self._boundaries) == 0:  # No open boundary at all?
            # TODO: Create a land boundary
            raise Exception("A case without open boundary "
                            "is not supported yet.")
        else:
            n_open_bounds = self.n_boundaries(trimesh.OPEN_BOUNDARY)
            for i in range(n_open_bounds):
                ns = []
                open_boundary = self._boundaries[i]
#                 print open_boundary.nodes
                last = open_boundary.nodes[-1]
                # Check if there is another boundary right next this
                if not self._check_if_beginning_of_boundary(last):
                    ns.append(last)
                    done = False
                    while not done:
#                         print "last node: ", last
                        next_node = self._get_next_node_on_boundary(last)
                        ns.append(next_node)
                        if self._check_if_beginning_of_boundary(next_node):
                            done = True
                            self.add_boundary(ns, trimesh.LAND_BOUNDARY)
                        else:
                            last = next_node

    def _check_if_beginning_of_boundary(self, node_i):
        """ Check if the given node is any of starting nodes of boundary
            node strings.
        """
        for boundary in self._boundaries:
            if boundary.btype == trimesh.OPEN_BOUNDARY:
                if node_i == boundary.nodes[0]:
                    return True
        return False

    def _get_not_assigned_boundary_edges(self):
        """ get any edge that is not assigned
        """
        notypes = []
        i = 0
        for edge in self._edges:
            if edge[2] == trimesh.BOUNDARY_EDGE:
                notypes.append(i)
            i += 1
        return notypes
