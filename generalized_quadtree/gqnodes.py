
class QtInteriorNode:

    def __init__(self, prefix, level):
        self.prefix = prefix
        self.level = level
        self.data = {}
        self.children = {}

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

    def __init__(self, prefix, name, info):
        self.prefix = prefix
        self.data = {name: info}

    def add(self, name, info):
        self.data[name] = info

    def add_leaf(self, leaf):
        for (name, info) in leaf.data.items():
            self.add(name, info)

    def list_dump(self, tree):
        return ("Leaf " + tree.qs(self.prefix), self.data)
