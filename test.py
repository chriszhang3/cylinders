"""Unit Tests"""

import unittest
from surface_dynamics import CylinderDiagram
from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
from surface_dynamics import AbelianStratum
from lib import check_pants_condition, list_partitions, \
                  filter_pants_condition, \
                  filter_homologous_condition, filter_leaf_condition, \
                  filter_standard_twist_condition, find_generalized_pants
from Graph import CylinderGraph
from Twist import Twist


class Test(unittest.TestCase):

    def test_find_generalized_pants(self):
        cd = CylinderDiagram("(0)-(2) (1,2,3)-(4,5) (4)-(3) (5)-(0,1)")
        cd_g = CylinderGraph(cd).digraph
        self.assertEqual(find_generalized_pants(cd_g), set([frozenset([1, 2, 3])]))
        cd = CylinderDiagram("(0,3)-(5) (1)-(0) (2,5)-(3,4) (4)-(1,2)")
        cd_g = CylinderGraph(cd).digraph
        self.assertEqual(find_generalized_pants(cd_g), set())
        cd = CylinderDiagram("(0,2,1)-(3,4,5) (3)-(1) (4)-(2) (5)-(0)")
        cd_g = CylinderGraph(cd).digraph
        self.assertEqual(find_generalized_pants(cd_g), set([frozenset([0, 1, 2, 3])]))
        cd = CylinderDiagram("(0)-(2) (1)-(3) (2,4,3)-(5,6) (5)-(4) (6)-(0,1)")
        cd_g = CylinderGraph(cd).digraph
        self.assertEqual(find_generalized_pants(cd_g), set([frozenset({0, 1, 4}), frozenset({2, 3, 4}), frozenset({0, 1, 2, 3})]))

    def test_check_pants_condition(self):
        self.assertTrue(check_pants_condition([{1}, {2}, {3}], [[1, 2, 3]]))
        self.assertFalse(check_pants_condition([{0, 1}, {2, 3}], [[1, 2, 3]]))
        self.assertTrue(check_pants_condition([{0, 1, 2, 3}], [[1, 2, 3]]))
    
    def test_filter_partition(self):
        cd = CylinderDiagram("(0,3)-(5) (1)-(2) (2,5)-(3,4) (4)-(0,1)")
        partitions = list_partitions(4, 2)
        self.assertEqual(len(filter_pants_condition(cd, partitions)), 1)
        cd = CylinderDiagram("(0,3)-(0,5) (1,2)-(1,4) (4)-(3) (5)-(2)")
        partitions = list_partitions(4, 2)
        self.assertEqual(len(filter_pants_condition(cd, partitions)), 7)
        cd = CylinderDiagram("(0,2,1)-(3,4,5) (3)-(1) (4)-(2) (5)-(0)")
        partitions = list_partitions(4, 2)
        self.assertEqual(len(filter_pants_condition(cd, partitions)), 0)

    def test_partitions(self):
        self.assertEqual(len(list_partitions(5, 2, singletons=False)), 10)
        self.assertEqual(len(list_partitions(5, 2, singletons=True)), 15)
        self.assertEqual(len(list_partitions(6, 2, singletons=False)), 25)
        self.assertEqual(len(list_partitions(6, 3, singletons=False)), 15)
    
    def test_pants_relations(self):
        C = CylinderDiagrams()
        H = AbelianStratum(3, 1).components()[0]
        valid_classes = [cd for cd in C.get_iterator(H, 4)
                         if filter_pants_condition(cd, list_partitions(4, 2, False))]
        self.assertFalse(valid_classes)

        H = AbelianStratum(2, 2).components()[1]
        valid_classes = [cd for cd in C.get_iterator(H, 4)
                         if filter_pants_condition(cd, list_partitions(4, 2, False))]
        self.assertEqual(len(valid_classes), 4)

        H = AbelianStratum(2, 1, 1).components()[0]
        valid_classes = [cd for cd in C.get_iterator(H, 4)
                         if filter_pants_condition(cd, list_partitions(4, 2, False))]
        self.assertEqual(len(valid_classes), 9)
    
    def test_find_homologous_cylinders(self):
        cd = CylinderDiagram('(0,2,1)-(4,5) (3,5)-(0,2,1) (4)-(3)')
        tw = Twist(cd)
        self.assertEqual(set(tw.find_homologous_cylinders()[0]), {0, 1})
        cd = CylinderDiagram('(0,3)-(0,5) (1,2)-(1,4) (4,6)-(3,7) (5,7)-(2,6)')
        tw = Twist(cd)
        self.assertEqual(set(tw.find_homologous_cylinders()[0]), {2, 3})
        cd = CylinderDiagram('(0)-(2) (1,2,3)-(4,5) (4)-(3) (5)-(0,1)')
        tw = Twist(cd)
        self.assertEqual(tw.find_homologous_cylinders(), [])
    
    def test_filter_homologous_cylinders(self):
        part_list = list_partitions(3, 2)
        cd = CylinderDiagram('(0,2,1)-(4,5) (3,5)-(0,2,1) (4)-(3)')
        good_part = filter_homologous_condition(cd, part_list)
        self.assertEqual(len(good_part), 1)
        part_list = list_partitions(4, 2)
        cd = CylinderDiagram('(0,3)-(0,5) (1,2)-(1,4) (4,6)-(3,7) (5,7)-(2,6)')
        good_part = filter_homologous_condition(cd, part_list)
        self.assertEqual(len(good_part), 3)
        part_list = list_partitions(5, 2)
        cd = CylinderDiagram('(0)-(2) (1)-(3) (2,4,3)-(5,6) (5)-(4) (6)-(0,1)')
        good_part = filter_homologous_condition(cd, part_list)
        self.assertEqual(len(good_part), 15)
        cd = CylinderDiagram('(0,1)-(0,6) (2)-(5) (3)-(4) (4,5)-(1) (6)-(2,3)')
        good_part = filter_homologous_condition(cd, part_list)
        self.assertEqual(len(good_part), 7)

    def test_filter_leaf_condition(self):
        part_list = list_partitions(4, 2)
        cd = CylinderDiagram("(0,1)-(0,2) (2)-(3) (3,4)-(1,5) (5)-(4)")
        good_part = filter_leaf_condition(cd, part_list)
        self.assertEqual(len(good_part), 4)
        cd = CylinderDiagram("(0)-(1) (1,3,4,2)-(5,6) (5)-(0,4) (6)-(2,3)")
        good_part = filter_leaf_condition(cd, part_list)
        self.assertEqual(len(good_part), 7)
        cd = CylinderDiagram("(0,3)-(5) (1)-(2) (2,5)-(3,4) (4)-(0,1)")
        good_part = filter_leaf_condition(cd, part_list)
        self.assertEqual(len(good_part), 7)
    
    def test_filter_standard_twist_condition(self):
        part_list = list_partitions(4, 3)
        cd = CylinderDiagram("(0,3)-(5) (1)-(0) (2,5)-(3,4) (4)-(1,2)")
        good_part = filter_standard_twist_condition(cd, part_list)
        good_part = [set(part) for part in good_part]
        possible = {frozenset({0, 3}), frozenset({1}), frozenset({2})}
        self.assertTrue(possible in good_part)

        part_list = list_partitions(4, 3)
        cd = CylinderDiagram("(0,2,1)-(3,4,5) (3)-(1) (4)-(2) (5)-(0)")
        good_part = filter_standard_twist_condition(cd, part_list)
        good_part = [set(part) for part in good_part]
        possible = {frozenset({0}), frozenset({1, 3}), frozenset({2})}
        self.assertTrue(possible in good_part)
        possible = {frozenset({0}), frozenset({1, 2}), frozenset({3})}
        self.assertTrue(possible in good_part)
        possible = {frozenset({0}), frozenset({1}), frozenset({2, 3})}
        self.assertTrue(possible in good_part)

        part_list = list_partitions(5, 3)
        cd = CylinderDiagram("(0,2)-(6) (1)-(3) (3,6)-(4,5) (4)-(0) (5)-(1,2)")
        good_part = filter_standard_twist_condition(cd, part_list)
        good_part = [set(part) for part in good_part]
        possible = {frozenset({0, 4}), frozenset({1, 3}), frozenset({2})}
        self.assertTrue(possible in good_part)

        part_list = list_partitions(5, 3)
        cd = CylinderDiagram("(0,6)-(4,5) (1,2)-(3,6) (3)-(2) (4)-(1) (5)-(0)")
        good_part = filter_standard_twist_condition(cd, part_list)
        good_part = [set(part) for part in good_part]
        possible = {frozenset({0, 1}), frozenset({2, 4}), frozenset({3})}
        self.assertTrue(possible in good_part)

if __name__ == "__main__":
    unittest.main()
