
# Todo:
#  GraphQuadTree
#  centroid and names accumulation at nodes
#  node distance to point
#  add graph node heuristic
#  node attraction heuristic

from . import gqnodes
import pprint
import numpy as np
from numpy.linalg import norm

def qs(index, levels, dimensions=2):
    "format quadtree index as binary string"
    length = levels * dimensions
    formatter = "{0:0>" + str(length) + "b}"
    result = formatter.format(index)
    assert len(result) == length
    return "0b" + result

class GeneralizedQuadtree:

    def __init__(self, origin, sidelength, levels):
        self.root = None
        # minimum position of the volume
        self.origin = np.array(origin)
        # maximum length of each side of the volume
        self.sidelength = sidelength
        self.center = self.origin + self.sidelength * 0.5
        # maximum number of levels in the tree
        self.levels = levels
        # number of dimensions of the volume
        self.dimensions = len(origin)
        # number of voxels on a side in the volume.
        self.int_side = 2 ** levels
        # quadrant fan-out at each node
        self.nquadrants = 2 ** self.dimensions
        self.nquadrants1 = self.nquadrants - 1
        # side length of a voxel
        self.min_side = float(sidelength) / self.int_side

    def quadrant_indices(self, index, level):
        """
        Sequence of quadrant indices at the next level below index.
        """
        dimensions = self.dimensions
        levels = self.levels
        shift = dimensions * (levels - (level + 1))
        assert shift >= 0, "no quadrants below leaf level " + repr((levels, level))
        # Low order bits should be zeros
        assert (index << (shift + dimensions)) & ((1 << (levels * dimensions)) - 1) == 0, (
            "bad index " + repr((index, level, shift, self.qs(index)))
        )
        step = 1 << shift
        limit = (self.nquadrants << shift) + index
        return xrange(index, limit, step)

    def walk(self, callback, data=None):
        "walk reverse breadth first passing (node, tree, data) to callback."
        if self.root is None:
            return   # Do nothing if tree is empty.
        self.root.walk(self, callback, data)

    def adjacent(self, index1, index2, level, pos1=None, pos2=None):
        """
        Test whether quadrants at level containing indices have abs
        common border or vertex.
        """
        shift = self.levels - level
        if pos1 is None:
            pos1 = self.index_to_index_position(index1)
        if pos2 is None:
            pos2 = self.index_to_index_position(index2)
        level_offset = (pos1 >> shift) - (pos2 >> shift)
        max_offset = np.max(np.abs(level_offset))
        return max_offset <= 1

    def adjacency_walk(self, position, callback, data=None):
        """
        Cover the quadtree.
        For each leaf visit the leaf or a node containing the leaf once.
        recurse into quadrants that are adjacent at the same level to
        the quadrant containing position.
        call callback(position, node, tree, data) at non-recursed nodes.
        """
        if self.root is None:
            return None
        iposition = self.int_position(position)
        #print "iposition", iposition
        self.root.adjacency_walk(self, callback, data, position, iposition)

    def index_corner(self, index):
        voxels = int_index_inverse(index, self.levels, self.dimensions)
        return np.array(voxels) * self.min_side + self.origin

    def level_side(self, level):
        return self.sidelength / float(2 ** level)

    def avg_dist_to_quadrant_point(self, index, level, location):
        "Heuristic average L1 distance to a point in the quadrant from location."
        side = self.level_side(level)
        corner = self.index_corner(index)
        ur_corner = corner + side
        location = np.array(location)
        offset = max(0, max(corner - location), max(location - ur_corner))
        return offset + side/2.0

    def avg_dist_between_quadrants(self, index1, level1, index2, level2):
        """
        Heuristic L1 distance between two points in quadrants.
        """
        side1 = self.level_side(level1)
        side2 = self.level_side(level2)
        corner1 = self.index_corner(index1)
        corner2 = self.index_corner(index2)
        ur_corner1 = corner1 + side1
        ur_corner2 = corner2 + side2
        offset = 0  # default
        intersecting = True
        for index in (0, 1):
            d1 = corner1[index] - ur_corner2[index]
            d2 = corner2[index] - ur_corner1[index]
            if d1 >= 0 or d2 >= 0:
                intersecting = False
            offset = max(offset, corner1[index] - ur_corner2[index], corner2[index] - ur_corner1[index])
        if intersecting:
            return 0.5 * max(side1, side2)
        return offset + 0.5 * (side1 + side2)

    def ppnode(self, node):
        pprint.pprint(node.list_dump(self))

    def qs(self, index):
        "Format index for tree parameters as binary string."
        return qs(index, self.levels, self.dimensions)

    def quad_string(self, quadrant):
        "Format quadrant index as binary string"
        return qs(quadrant, 1, self.dimensions)

    def list_dump(self):
        root = self.root
        if root is None:
            return None
        return root.list_dump(self)

    def add(self, at_position, name, info=None):
        if info is None:
            info = {}
        (pos_index, int_pos) = self.index_position(at_position)
        info = info.copy()
        assert "position" not in info, "position dict key is reserved " + repr(info)
        info["position"] = at_position
        leaf = gqnodes.QtLeafNode(pos_index, name, info)
        self.root = self.combine(self.root, leaf)

    def add_at_min_penalty(self, node_penalty_fn, name, info=None, initial_penalty_fn=None, normalize=None):
        
        # p "insert", (name, info)
        root = self.root
        if self.root is None:
            # add node at center
            # p "adding initial node"
            return self.add(self.center, name, info)
        index = 0
        index_to_node = {0: self.root}
        for level in range(self.levels):
            # expand any adjacent nodes at this level
            for node_index in index_to_node.keys():
                # p "looking for nodes to expand for", self.qs(index)
                node = index_to_node[node_index]
                if (isinstance(node, gqnodes.QtInteriorNode) and
                    node.level<=level and
                    self.adjacent(index, node_index, node.level)):
                    # p "    expanding", node.level, self.qs(node.prefix)
                    del index_to_node[node_index]
                    children = node.children
                    for quadrant in children:
                        child = children[quadrant]
                        # p "   expand child", child.level, self.qs(child.prefix)
                        index_to_node[child.prefix] = child
            # find index of best quadrant
            best_penalty = None
            best_corner = None
            for qindex in self.quadrant_indices(index, level):
                voxels = int_index_inverse(qindex, self.levels, self.dimensions)
                corner = self.index_corner(qindex)
                # p "    qindex", self.qs(qindex), voxels, corner, "at level", level
                total_penalty = 0
                if initial_penalty_fn is not None:
                    total_penalty = initial_penalty_fn(qindex, voxels, corner)
                for node_index in index_to_node:
                    node = index_to_node[node_index]
                    node_penalty = node_penalty_fn(node, qindex, voxels, corner)
                    total_penalty += node_penalty
                    # p "        penalty", self.qs(node.prefix), node_penalty, total_penalty, node.get_names()
                penalty = total_penalty
                if normalize:
                    penalty = normalize(total_penalty)
                # p "    total penalty", penalty
                if best_penalty is None or penalty < best_penalty:
                    # p "    choosing", self.qs(qindex)
                    best_penalty = penalty
                    best_corner = corner
                    index = qindex
                elif penalty == best_penalty:
                    # prefer nearer to the center, deterministic choice (for testing)
                    new_norm = norm(self.center - corner)
                    best_norm = norm(self.center - best_corner)
                    if new_norm < best_norm or (new_norm == best_norm and qindex < index):
                        # p "    preferring", self.qs(qindex)
                        best_penalty = penalty
                        best_corner = corner
                        index = qindex
        position = self.index_corner(index)
        # p "final position", self.qs(index), position
        self.add(position, name, info)

    def quadrant(self, index, at_level):
        assert at_level >= 0
        assert at_level <= self.levels
        dimensions = self.dimensions
        nquadrants1 = self.nquadrants1
        shift = (self.levels - at_level) * dimensions
        #return (index >> shift) & (self.nquadrants - 1)
        hi_bits = (index >> shift)
        quadrant = hi_bits & nquadrants1
        prefix = (hi_bits & ~nquadrants1) << shift
        return (prefix, quadrant)

    def combine(self, node, leaf):
        if node is None:
            return leaf
        nprefix = node.prefix
        lprefix = leaf.prefix
        levels = self.levels
        if isinstance(node, gqnodes.QtLeafNode):
            (cprefix, clevel) = self.common_prefix_level(nprefix, lprefix)
            if clevel == levels:
                # Leaf collision: extend leaf data at node.
                node.add_leaf(leaf)
                return node
            # Otherwise create a new parent for the leaves
            result = gqnodes.QtInteriorNode(cprefix, clevel)
            result.add_new_child(node, self)
            result.add_new_child(leaf, self)
            return result
        assert isinstance(node, gqnodes.QtInteriorNode)
        nlevel = node.level
        (cprefix, clevel) = self.common_prefix_level(nprefix, lprefix)
        if clevel >= nlevel:
            # insert below node
            node.add_leaf(leaf, self)
            return node
        # otherwise create a new parent for the leaf and node
        result = gqnodes.QtInteriorNode(cprefix, clevel)
        result.add_new_child(node, self)
        result.add_new_child(leaf, self)
        return result

    def int_position(self, position):
        "Convert a position to integer coordinates relative to the origin"
        assert len(position) == self.dimensions, (
            "Bad dimensions in position " + repr(position)
        )
        origin = self.origin
        min_side = self.min_side
        int_side = self.int_side
        # XXXX This permits "negative positions" that round to zero. Bug?
        result = [int((position[i] - origin_i) / min_side)
                  for (i, origin_i) in enumerate(origin)]
        return np.array(result)

    def index_position(self, position):
        """
        Quadtree index of position, and integer position.
        """
        int_position = self.int_position(position)
        index = int_index(int_position, self.levels)
        return (index, int_position)

    def index(self, position):
        return self.index_position(position)[0]

    def index_to_index_position(self, index):
        """
        Convert index back to int position relative to origin
        """
        return int_index_inverse(
            index, self.levels, self.dimensions)

    def common_prefix_level(self, index1, index2, from_level=0):
        """
        quadrant coordinates which agree between index1 and index2.
        """
        dim = self.dimensions
        levels = self.levels
        shift_level = from_level
        while True:
            shift = shift_level * dim
            prefix1 = index1 >> shift
            prefix2 = index2 >> shift
            if prefix1 == prefix2:
                return (prefix1 << shift, levels - shift_level)
            if shift_level < levels:
                shift_level += 1
            else:
                raise ValueError, ("too many levels "
                    + repr((levels, index1, index2)))

def int_index_inverse(index, levels, dimensions):
    index0 = index  # for diagnositics only
    result = [0] * dimensions
    d_order = range(dimensions)
    #pr "start index {0:b}".format(index)
    for level in range(levels):
        for dim in d_order:
            bit = index & 1
            #property "level, dim, bit, before", level, dim, bit, "{0:b}".format(result[dim])
            index = index >> 1
            result[dim] = result[dim] | (bit << level)
            #pr "after {0:b} {1:b}".format(index, result[dim])
    assert index == 0, "too many bits " + repr((index0, levels, dimensions))
    return np.array(result)

def int_index(position_ints, levels):
    result = 0
    shift = 0
    p = list(position_ints)
    assert min(p) >= 0, "negative entries " + repr(p)
    dimensions = list(range(len(p)))
    shift = 0
    for level in range(levels):
        for d in dimensions:
            p_d = p[d]
            bit = p_d & 1
            p[d] = p_d >> 1
            result = result | (bit << shift)
            shift += 1
    assert max(p) == 0, "unshifted bits " + repr((p, position_ints))
    return result
