import unittest
from surface_dynamics import CylinderDiagram
from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
from surface_dynamics import AbelianStratum
from lib import check_pants_condition, list_partitions, \
                  filter_pants_condition, find_homologous_cylinders, \
                  filter_homologous_condition, filter_leaf_condition


class Test(unittest.TestCase):

    def test_check_pants_condition(self):
        self.assertTrue(check_pants_condition([{1}, {2}, {3}], [[1, 2, 3]]))
        self.assertFalse(check_pants_condition([{0, 1}, {2, 3}], [[1, 2, 3]]))
    
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
        cd2 = CylinderDiagram('(0,3)-(0,5) (1,2)-(1,4) (4,6)-(3,7) (5,7)-(2,6)')
        self.assertEqual(find_homologous_cylinders(cd)[0], [0, 1])
        self.assertEqual(find_homologous_cylinders(cd2)[0], [2, 3])
    
    def test_filter_homologous_cylinders(self):
        part_list = list_partitions(3, 2)
        cd = CylinderDiagram('(0,2,1)-(4,5) (3,5)-(0,2,1) (4)-(3)')
        good_part = filter_homologous_condition(cd, part_list)
        self.assertEqual(len(good_part), 1)
        part_list = list_partitions(4, 2)
        cd = CylinderDiagram('(0,3)-(0,5) (1,2)-(1,4) (4,6)-(3,7) (5,7)-(2,6)')
        good_part = filter_homologous_condition(cd, part_list)
        self.assertEqual(len(good_part), 3)

    def test_filter_leaf_condition(self):
        part_list = list_partitions(4, 2)
        cd = CylinderDiagram("(0,1)-(0,2) (2)-(3) (3,4)-(1,5) (5)-(4)")
        good_part = filter_leaf_condition(cd, part_list)
        self.assertEqual(len(good_part), 4)

if __name__ == "__main__":
    unittest.main()