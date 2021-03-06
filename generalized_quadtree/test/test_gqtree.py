import unittest
from .. import gqtree
from ..gqtree import qs
import pprint

class TestDR(unittest.TestCase):

    def test_qs(self):
        self.assertEqual(gqtree.qs(0b00011000, 4, 2), "0b00011000")

    def test_int_index(self):
        self.assertEqual(gqtree.int_index([0b01, 0b10], 2), 0b1001)
        self.assertEqual(gqtree.int_index([0b01, 0b11, 0b10], 2), 0b110011)
        #Sprint ("{0:b}".format(test))
        self.assertEqual(gqtree.int_index([0b1100, 0b0000], 4), 0b01010000)
        self.assertEqual(
            gqtree.int_index([0b11100, 0b00011, 0b00011], 5), 0b001001001110110)
        self.assertEqual(
            gqtree.int_index([0b11100, 0b11100, 0b11100], 5), 0b111111111000000)
        with self.assertRaises(AssertionError):
            test = gqtree.int_index([0b0111, 0b01], 2)
        with self.assertRaises(AssertionError):
            test = gqtree.int_index([-0b011, 0b01], 2)

    def test_int_position(self):
        gq = gqtree.GeneralizedQuadtree(origin=[1.0, 2.0], sidelength=8.0, levels=2)
        self.assertEqual(list(gq.int_position([1.0, 2.0])), [0, 0])
        self.assertEqual(list(gq.int_position([8, 9])), [3, 3])

    def test_index(self):
        gq = gqtree.GeneralizedQuadtree(origin=[1.0, 2.0], sidelength=8.0, levels=2)
        self.assertEqual(gq.index([1.0, 2.0]), 0x0000)
        self.assertEqual(gq.index([8, 9]), 0b1111)
        #print "neg index", gq.index([0, 1])
        with self.assertRaises(AssertionError):
            test = gq.index([9, 3])
        with self.assertRaises(AssertionError):
            test = gq.index([3, 0])
            #print ("{0:b}".format(test))

    def test_index_corner(self):
        gq = gqtree.GeneralizedQuadtree(origin=[0, 0], sidelength=8.0, levels=2)
        for (index, expected) in [
                (0b0000, [0.0, 0.0]),
                (0b0001, [2.0, 0.0]),
                (0b0010, [0.0, 2.0]),
                (0b0011, [2.0, 2.0]),
                (0b0100, [4.0, 0.0]),
                (0b0101, [6.0, 0.0]),
                (0b0110, [4.0, 2.0]),
                (0b0111, [6.0, 2.0]),
                (0b1000, [0.0, 4.0]),
                (0b1001, [2.0, 4.0]),
                (0b1010, [0.0, 6.0]),
                (0b1011, [2.0, 6.0]),
                (0b1100, [4.0, 4.0]),
                (0b1101, [6.0, 4.0]),
                (0b1110, [4.0, 6.0]),
                (0b1111, [6.0, 6.0]),
            ]:
            self.assertEqual(list(gq.index_corner(index)), expected)

    def test_index_inverse(self):
        def round_trip(p_int, levels):
            dims = len(p_int)
            forward = gqtree.int_index(p_int, levels)
            #pr [gqtree.qs(x, levels, 1) for x in p_int]
            #pr "forward is {0}".format(gqtree.qs(forward, levels, dims))
            back = gqtree.int_index_inverse(forward, levels, dims)
            self.assertEqual(list(p_int), list(back))
        round_trip((0, 1, 0, 0), 1)
        round_trip((0b10, 0b01), 2)
        round_trip((0b1011, 0b0110), 4)
        round_trip((0b10, 0b01, 0b11), 2)
        round_trip((0b101, 0b011, 0b000, 0b111, 0b110), 3)
        with self.assertRaises(AssertionError):
            t = gqtree.int_index_inverse(0b101010, 2, 2)

    def test_common_prefix(self):
        gq = gqtree.GeneralizedQuadtree(origin=[1.0, 2.0], sidelength=8.0, levels=2)
        self.assertEqual((0b0101, 2), gq.common_prefix_level(0b0101, 0b0101))
        self.assertEqual((0b0100, 1), gq.common_prefix_level(0b0101, 0b0111))
        self.assertEqual((0b0000, 0), gq.common_prefix_level(0b0101, 0b1111))
        with self.assertRaises(ValueError):
            test =  gq.common_prefix_level(0b010111, 0b111111)

    def test_common_from_level(self):
        gq = gqtree.GeneralizedQuadtree(origin=[1.0, 2.0], sidelength=8.0, levels=4)
        t = gq.common_prefix_level(0b01010101, 0b01010101, 2)
        self.assertEqual((0b01010000, 2), t)
        t = gq.common_prefix_level(0b01110101, 0b01010101, 2)
        self.assertEqual((0b01000000, 1), t)

    def test_quadrant(self):
        gq = gqtree.GeneralizedQuadtree(origin=[1.0, 2.0], sidelength=8.0, levels=2)
        self.assertEqual(gq.quadrant(0b1001, 1), (0b0000, 0b10))
        self.assertEqual(gq.quadrant(0b1001, 2), (0b1000, 0b01))
        with self.assertRaises(AssertionError):
            test = gq.quadrant(0b1001, 3)
        with self.assertRaises(AssertionError):
            test = gq.quadrant(0b1001, -1)

    def test_index2(self):
        gq = gqtree.GeneralizedQuadtree(origin=[0,0], sidelength=1.1, levels=2)
        lastindex = 0b1100
        for x in (0,1):
            for y in (0,1):
                xy = [x,y]
                index = gq.index(xy)
                # xy, gq.qs(index), gq.qs(lastindex)
                common = gq.common_prefix_level(index, lastindex)
                # "common", common
                self.assertEqual(common, (0, 0))
                lastindex = index

    def test_leaf_collision(self):
        "add 2 leaves descending from root"
        gq = gqtree.GeneralizedQuadtree(origin=[0,0], sidelength=2.1, levels=2)
        gq.add([0,0], "first")
        gq.add([0,0], "second")
        dump = gq.list_dump()
        #print "ADD2"
        #pprint.pprint(dump)
        expected = ('Leaf 0b0000',
                    {'first': {'position': [0, 0]}, 'second': {'position': [0, 0]}})
        self.assertEqual(dump, expected)

    def test_add2(self):
        "add 2 leaves descending from root"
        gq = gqtree.GeneralizedQuadtree(origin=[0,0], sidelength=2.1, levels=2)
        gq.add([0,0], "first")
        gq.add([1,1], "second")
        dump = gq.list_dump()
        #print "ADD2"
        #pprint.pprint(dump)
        expect = [
            'node 0b0000 LV1',
            {},
            [('0b00', ('Leaf 0b0000', {'first': {'position': [0, 0]}})),
             ('0b11', ('Leaf 0b0011', {'second': {'position': [1, 1]}}))]]
        self.assertEqual(dump, expect)

    def test_add3_below(self):
        "add 2 leaves descending from root and a third in the same quadrant as first."
        gq = gqtree.GeneralizedQuadtree(origin=[0,0], sidelength=2.1, levels=2)
        gq.add([0,0], "first")
        gq.add([1,1], "second")
        gq.add([0,0], "third")
        dump = gq.list_dump()
        expect = [
            'node 0b0000 LV1',
            {},
            [('0b00',
              ('Leaf 0b0000',
                {'first': {'position': [0, 0]}, 'third': {'position': [0, 0]}})),
             ('0b11', ('Leaf 0b0011', {'second': {'position': [1, 1]}}))]]
        self.assertEqual(dump, expect)

    def test_add3_above(self):
        "add 2 leaves descending from root and a third at higher level."
        gq = gqtree.GeneralizedQuadtree(origin=[0,0], sidelength=2.1, levels=2)
        gq.add([0,0], "first")
        gq.add([1,1], "second")
        gq.add([2,0], "third")
        dump = gq.list_dump()
        #print "ADD3"
        #pprint.pprint(dump)
        expect = [
             'node 0b0000 LV0',
             {},
             [('0b00',
               ['node 0b0000 LV1',
                {},
                [('0b00', ('Leaf 0b0000', {'first': {'position': [0, 0]}})),
                 ('0b11', ('Leaf 0b0011', {'second': {'position': [1, 1]}}))]]),
              ('0b01', ('Leaf 0b0101', {'third': {'position': [2, 0]}}))]]
        self.assertEqual(dump, expect)

    def test_add(self):
        gq = gqtree.GeneralizedQuadtree(origin=[1.0, 2.0], sidelength=8.0, levels=2)
        dump = gq.list_dump()
        self.assertEqual(dump, None)
        gq.add([7.1, 4.2], "first")
        dump = gq.list_dump()
        #pprint.pprint(dump)
        expect = ('Leaf 0b0111', {'first': {'position': [7.1, 4.2]}})
        self.assertEqual(dump, expect)
        # non-collision
        gq.add([7.1,2.2], "second")
        dump = gq.list_dump()
        #pprint.pprint(dump)
        expected = ['node 0b0100 LV1',
             {},
             [('0b01', ('Leaf 0b0101', {'second': {'position': [7.1, 2.2]}})),
              ('0b11', ('Leaf 0b0111', {'first': {'position': [7.1, 4.2]}}))]]
        self.assertEqual(dump, expected)
        # leaf collision
        gq.add([7.5,0.25], "third")
        dump = gq.list_dump()
        #pprint.pprint(dump)
        expected = ['node 0b0100 LV1',
                     {},
                     [('0b01',
                       ('Leaf 0b0101',
                        {'second': {'position': [7.1, 2.2]},
                         'third': {'position': [7.5, 0.25]}})),
                      ('0b11', ('Leaf 0b0111', {'first': {'position': [7.1, 4.2]}}))]]
        self.assertEqual(expected, dump)
        gq.add([1.5,2.25], "fourth")
        expected = ['0b0100', 1, {},
            [('0b00', ('0b0000', {'fourth': {'position': [1.5, 2.25]}})),
             ('0b01',
              ('0b0101',
               {'second': {'position': [7.1, 2.2]},
                'third': {'position': [7.5, 0.25]}})),
             ('0b11', ('0b0111', {'first': {'position': [7.1, 4.2]}}))]]
        gq.add([7.1, 7.1], "fifth")
        dump = gq.list_dump()
        #pprint.pprint(dump)
        expected = ['node 0b0000 LV0',
                     {},
                     [('0b00', ('Leaf 0b0000', {'fourth': {'position': [1.5, 2.25]}})),
                      ('0b01',
                       ['node 0b0100 LV1',
                        {},
                        [('0b01',
                          ('Leaf 0b0101',
                           {'second': {'position': [7.1, 2.2]},
                            'third': {'position': [7.5, 0.25]}})),
                         ('0b11', ('Leaf 0b0111', {'first': {'position': [7.1, 4.2]}}))]]),
                      ('0b11', ('Leaf 0b1101', {'fifth': {'position': [7.1, 7.1]}}))]]
        self.assertEqual(dump, expected)

    def test_walk(self):
        gq = gqtree.GeneralizedQuadtree(origin=[1.0, 2.0], sidelength=8.0, levels=2)
        gq.add([2,7], "f", {"w": 10})
        gq.add([2,3], "s", {"w": 100})
        gq.add([7,7], "t", {"w": 1000})
        D = {}
        def callback(node, tree, data):
            ltotal = sum(i["w"] for i in node.data.values())
            ctotal = sum(c.sum_w
                         for c in node.children.values())
            node.sum_w = ctotal + ltotal
            D[(node.level, node.prefix)] = (node.sum_w, node.data)
        gq.walk(callback)
        expectD = {(None, 0): (100, {'s': {'position': [2, 3], 'w': 100}}),
                   (None, 8): (10, {'f': {'position': [2, 7], 'w': 10}}),
                   (None, 13): (1000, {'t': {'position': [7, 7], 'w': 1000}}),
                   (0, 0): (1110, {})}
        self.assertEqual(D, expectD)

    def test_adjacency_walk(self):
        gq = gqtree.GeneralizedQuadtree(origin=[1.0, 2.0], sidelength=8.0, levels=5)
        gq.add([2,2], "f", {"w": 10})   # summarize
        gq.add([2,3], "s", {"w": 100})  # summarize
        gq.add([7,7], "t", {"w": 1000})  # visit
        D = {}
        def callback(p, node, t, d):
            D[(node.level, node.prefix)] = node.data
        gq.adjacency_walk((7.1, 7.1), callback)
        expectD = {(None, 864): {'t': {'position': [7, 7], 'w': 1000}}, 
                   (2, 0): {}}
        self.assertEqual(D, expectD)

    def test_get_names(self):
        gq = gqtree.GeneralizedQuadtree(origin=[1.0, 2.0], sidelength=8.0, levels=5)
        gq.add([2,2], "f", {"w": 10})   # summarize
        names = gq.root.get_names()
        self.assertEqual(set(["f"]), names)
        gq.add([2,3], "s", {"w": 100})  # summarize
        names = gq.root.get_names()
        self.assertEqual(set(["f", "s"]), names)
        gq.add([7,7], "t", {"w": 1000})  # visit
        names = gq.root.get_names()
        self.assertEqual(set(["f", "s", "t"]), names)

    def test_avg_dist(self):
        gq = gqtree.GeneralizedQuadtree(origin=[1.0, 2.0], sidelength=2.0, levels=2)
        for loc in [(0, 1), (0, 2), (4, 3), (3, 1)]:
            dist = gq.avg_dist_to_quadrant_point(0, 0, loc)
            self.assertEqual(dist, 2, repr(loc) + " should have dist 2")
        for loc in [(1,2), (3, 2), (3, 4), (2, 3)]:
            dist = gq.avg_dist_to_quadrant_point(0, 0, loc)
            self.assertEqual(dist, 1, repr(loc) + " should have dist 1")
        index = gq.index_position((2,2))[0]
        for loc in [(1, 2), (4, 2), (1, 1), (2, 4)]:
            dist = gq.avg_dist_to_quadrant_point(index, 1, loc)
            self.assertEqual(dist, 1.5, repr(loc) + " should have dist 1.5")
        for loc in [(2,2), (3, 2), (3, 3), (2.5, 2.5)]:
            dist = gq.avg_dist_to_quadrant_point(index, 1, loc)
            self.assertEqual(dist, 0.5, repr((loc, dist)) + " should have dist 0.5")

    def test_quad_dist(self):
        gq = gqtree.GeneralizedQuadtree(origin=[1.0, 2.0], sidelength=4.0, levels=2)
        avg = gq.avg_dist_between_quadrants
        L = [(gq.qs(index), avg(0, 2, index, 2)) for index in range(16)]
        eL = \
            [('0b0000', 0.5),
             ('0b0001', 1.0),
             ('0b0010', 1.0),
             ('0b0011', 1.0),
             ('0b0100', 2.0),
             ('0b0101', 3.0),
             ('0b0110', 2.0),
             ('0b0111', 3.0),
             ('0b1000', 2.0),
             ('0b1001', 2.0),
             ('0b1010', 3.0),
             ('0b1011', 3.0),
             ('0b1100', 2.0),
             ('0b1101', 3.0),
             ('0b1110', 3.0),
             ('0b1111', 3.0)]
        self.assertEqual(L, eL)
        L = [(gq.qs(index), avg(0, 1, index, 2)) for index in range(16)]
        eL = \
           [('0b0000', 1.0),
            ('0b0001', 1.0),
            ('0b0010', 1.0),
            ('0b0011', 1.0),
            ('0b0100', 1.5),
            ('0b0101', 2.5),
            ('0b0110', 1.5),
            ('0b0111', 2.5),
            ('0b1000', 1.5),
            ('0b1001', 1.5),
            ('0b1010', 2.5),
            ('0b1011', 2.5),
            ('0b1100', 1.5),
            ('0b1101', 2.5),
            ('0b1110', 2.5),
            ('0b1111', 2.5)]
        self.assertEqual(L, eL)
        L = [(gq.qs(index), avg(0b1100, 2, index, 2)) for index in range(16)]
        eL = \
           [('0b0000', 2.0),
            ('0b0001', 2.0),
            ('0b0010', 2.0),
            ('0b0011', 1.0),
            ('0b0100', 2.0),
            ('0b0101', 2.0),
            ('0b0110', 1.0),
            ('0b0111', 1.0),
            ('0b1000', 2.0),
            ('0b1001', 1.0),
            ('0b1010', 2.0),
            ('0b1011', 1.0),
            ('0b1100', 0.5),
            ('0b1101', 1.0),
            ('0b1110', 1.0),
            ('0b1111', 1.0)]
        self.assertEqual(L, eL)

    def test_quadrant_indices(self):
        gq = gqtree.GeneralizedQuadtree(origin=[1.0, 2.0], sidelength=2.0, levels=2)
        L = list(gq.quadrant_indices(0b0100, 1))
        sL = [gq.qs(x) for x in L]
        eL = ['0b0100', '0b0101', '0b0110', '0b0111']
        self.assertEqual(sL, eL)
        gq = gqtree.GeneralizedQuadtree(origin=[1.0, 2.0], sidelength=2.0, levels=4)
        L = list(gq.quadrant_indices(0b01100000, 2))
        sL = [gq.qs(x) for x in L]
        eL = ['0b01100000', '0b01100100', '0b01101000', '0b01101100']
        self.assertEqual(sL, eL)
        with self.assertRaises(AssertionError):
            L = list(gq.quadrant_indices(0b01100000, 4))
        with self.assertRaises(AssertionError):
            L = list(gq.quadrant_indices(0b01100010, 2))

    def test_adjacent(self):
        gq = gqtree.GeneralizedQuadtree(origin=[0, 0], sidelength=8.0, levels=8)
        for (loc1, loc2) in [((3, 1), (4.1, 1)), ((3, 3), (4.1, 4.1)), ((6, 3), (6, 4.1))]:
            index1 = gq.index_position(loc1)[0]
            index2 = gq.index_position(loc2)[0]
            for level in range(4):
                self.assertEqual(gq.adjacent(index1, index2, level), True,
                    repr((index1, index2)) + " should be adjacent at " + repr(level))
            for level in range(4, 9):
                self.assertEqual(gq.adjacent(index1, index2, level), False,
                    repr((index1, index2)) + " should not be adjacent at " + repr(level))

    def test_add_min(self):
        def node_penalty_fn(node, qindex, voxels, corner):
            if qindex == node.prefix:
                return len(node.get_names())
            return 0
        gq = gqtree.GeneralizedQuadtree(origin=[0, 0], sidelength=8.0, levels=1)
        for name in list("abcd"):
            gq.add_at_min_penalty(node_penalty_fn, name)
        dump = gq.list_dump()
        expect = \
            ['node 0b00 LV0',
            {},
            [('0b00', ('Leaf 0b00', {'d': {'position': [0.0, 0.0]}})),
            ('0b01', ('Leaf 0b01', {'b': {'position': [4.0, 0.0]}})),
            ('0b10', ('Leaf 0b10', {'c': {'position': [0.0, 4.0]}})),
            ('0b11', ('Leaf 0b11', {'a': {'position': [4.0, 4.0]}}))]]
        self.assertEqual(expect, dump)

    def test_add_min2(self):
        def node_penalty_fn(node, qindex, voxels, corner):
            if qindex == node.prefix:
                return len(node.get_names())
            return 0
        gq = gqtree.GeneralizedQuadtree(origin=[0, 0], sidelength=8.0, levels=2)
        for name in (str(i) for i in range(16)):
            gq.add_at_min_penalty(node_penalty_fn, name)
        dump = gq.list_dump()
        expect = \
            ['node 0b0000 LV0',
            {},
            [('0b00',
            ['node 0b0000 LV1',
                {},
                [('0b00', ('Leaf 0b0000', {'15': {'position': [0.0, 0.0]}})),
                ('0b01', ('Leaf 0b0001', {'6': {'position': [2.0, 0.0]}})),
                ('0b10', ('Leaf 0b0010', {'11': {'position': [0.0, 2.0]}})),
                ('0b11', ('Leaf 0b0011', {'5': {'position': [2.0, 2.0]}}))]]),
            ('0b01',
            ['node 0b0100 LV1',
                {},
                [('0b00', ('Leaf 0b0100', {'9': {'position': [4.0, 0.0]}})),
                ('0b01', ('Leaf 0b0101', {'13': {'position': [6.0, 0.0]}})),
                ('0b10', ('Leaf 0b0110', {'1': {'position': [4.0, 2.0]}})),
                ('0b11', ('Leaf 0b0111', {'2': {'position': [6.0, 2.0]}}))]]),
            ('0b10',
            ['node 0b1000 LV1',
                {},
                [('0b00', ('Leaf 0b1000', {'10': {'position': [0.0, 4.0]}})),
                ('0b01', ('Leaf 0b1001', {'3': {'position': [2.0, 4.0]}})),
                ('0b10', ('Leaf 0b1010', {'14': {'position': [0.0, 6.0]}})),
                ('0b11', ('Leaf 0b1011', {'4': {'position': [2.0, 6.0]}}))]]),
            ('0b11',
            ['node 0b1100 LV1',
                {},
                [('0b00', ('Leaf 0b1100', {'0': {'position': [4.0, 4.0]}})),
                ('0b01', ('Leaf 0b1101', {'7': {'position': [6.0, 4.0]}})),
                ('0b10', ('Leaf 0b1110', {'8': {'position': [4.0, 6.0]}})),
                ('0b11', ('Leaf 0b1111', {'12': {'position': [6.0, 6.0]}}))]])]]
        self.assertEqual(expect, dump)
