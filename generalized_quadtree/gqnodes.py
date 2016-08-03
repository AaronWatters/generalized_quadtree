
#from . import gqtree
import numpy as np

class QtInteriorNode:

    _int_position = None
    _names = None   # names of all descendents

    def __init__(self, prefix, level):
        self.prefix = prefix
        self.level = level
        self.data = {}
        self.children = {}

    def get_names(self):
        names = self._names
        if names is None:
            names = set()
            for child in self.children.values():
                names.update(child.get_names())
            self._names = names
        return names

    def adjacency_walk(self, tree, callback, data, position, iposition):
        level = self.level
        my_ipos = self._int_position
        if my_ipos is None:  
            my_ipos = tree.index_to_index_position(self.prefix)
            self._int_position = my_ipos
        shift = tree.levels - level 
        level_offset = (my_ipos >> shift) - (iposition >> shift)
        max_offset = np.max(np.abs(level_offset))
        #ind = "    " * level
        #print ind, self._int_position, level, shift, "max", max_offset
        if max_offset <= 1:
            #print ind, "expand", level_offset
            # recursively expand node
            for node in self.children.values():
                node.adjacency_walk(tree, callback, data, position, iposition)
        else:
            #print ind, "callback", level_offset
            # visit node without expanding (not adjacent at level)
            callback(position, self, tree, data)

    def walk(self, tree, callback, data):
        "walk reverse breadth first passing (node, tree, data) to callback."
        children = self.children
        for quadrant in children:
            child = children[quadrant]
            child.walk(tree, callback, data)
        callback(self, tree, data)

    def add_new_child(self, node, tree):
        "Add a child in empty quadrant."
        level = self.level
        nprefix = node.prefix
        (remainder, quadrant) = tree.quadrant(nprefix, level + 1)
        assert remainder == self.prefix, ("bad child prefix" +
            repr(tree.qs(self.prefix), tree.qs(remainder), tree.qs(nprefix), level))
        children = self.children
        assert children.get(quadrant) == None, "non-empty quadrant " + repr(quadrant)
        children[quadrant] = node
        names = self._names
        if names is not None:
            names.update(node.get_names())

    def add_leaf(self, leaf, tree):
        level = self.level
        nprefix = self.prefix
        (remainder, quadrant) = tree.quadrant(leaf.prefix, level + 1)
        assert remainder == self.prefix, ("bad leaf prefix" +
            repr(tree.qs(self.prefix), tree.qs(remainder), tree.qs(nprefix), level))
        children = self.children
        old_child = children.get(quadrant)
        new_child = tree.combine(old_child, leaf)
        if isinstance(new_child, QtInteriorNode):
            assert new_child.level > level
        names = self._names
        if names is not None:
            names.update(leaf.get_names())
        children[quadrant] = new_child

    def list_dump(self, tree):
        children_dumped = []
        for (quadrant, child) in sorted(self.children.items()):
            if child:
                dumped = child.list_dump(tree)
                children_dumped.append((tree.quad_string(quadrant), dumped))
        return [
            "node %s LV%s" % (tree.qs(self.prefix), self.level),
            self.data,
            children_dumped]

class QtLeafNode:

    children = {}  # "read only constant"
    level = None  # "read only constant"

    def __init__(self, prefix, name, info):
        self.prefix = prefix
        self.data = {name: info}

    def get_names(self):
        return set(self.data)

    def adjacency_walk(self, tree, callback, data, position, iposition):
        # always visit any leaf that is reached.
        callback(position, self, tree, data)

    def walk(self, tree, callback, data):
        "walk reverse breadth first passing (node, tree, data) to callback."
        callback(self, tree, data)

    def add(self, name, info):
        self.data[name] = info

    def add_leaf(self, leaf):
        for (name, info) in leaf.data.items():
            self.add(name, info)

    def list_dump(self, tree):
        return ("Leaf " + tree.qs(self.prefix), self.data)
