from collections import defaultdict
import unittest
from surface_dynamics import CylinderDiagram
from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
from surface_dynamics import AbelianStratum
from sage.all import Partitions, SetPartitions, QQ, matrix, vector
from Graph import CylinderGraph

def find_element_in_partition(partition, element):
    """Gives the singularity corresponding to the right vertex of the saddle."""
    for i, part in enumerate(partition):
        if element in part:
            return i

def list_partitions(n, m, singletons=True):
    """List all ways to partition the set [1..n] into m sets.
    
    If singletons==False, do not allow singleton sets."""
    def contains_singleton(l):
        return any([i == 1 for i in l])

    partitions = []
    underlying_part = Partitions(n, length=m).list()
    if not singletons:
        underlying_part = [i for i in underlying_part if not contains_singleton(i)]
    for up in underlying_part:
        partitions.extend(SetPartitions(range(n), up))
    return partitions

def check_pants_condition(partition, pants_list):
    """Check the partition satisfies any homology restrictions coming from 
    the pants in pants_list.
    
    `partition` is a partition of the cylinders of horizontally periodic
    translation surface into M-parallel classes.

    `pants_list` is a list of pants, where each pants is a frozenset containing
    every cylinder in the pants
    
    The cylinders in the pants cannot be contained in exactly two distinct
    sets in `partition`.
    """
    for condition in pants_list:
        partition_sets = map(
            lambda c: find_element_in_partition(partition, c),
            condition
        )
        if len(set(partition_sets)) == 2:
            return False
    return True

def filter_pants_condition(cyl_diag, part_list):
    """Remove cylinder equivalence classes that don't satisfy constraints
    coming from generic pants.
    
    `cyl_diag` is a cylinder diagram
    `part` is a list of partitions to filter through
    returns a sublist  of `part` with unwanted ones removed
    """
    cyl_graph = CylinderGraph(cyl_diag)
    pants_list = cyl_graph.find_generic_pants()
    return [partition for partition in part_list 
                      if check_pants_condition(partition, pants_list)]

def find_homologous_cylinders(cyl_diag):
    """Find any pairs of homologous cylinders.

    `cyl_diag` is a cylinder diagram.
    returns a list of lists of homologous cylinders
    """
    cylinders = cyl_diag.cylinders()
    equations = {}
    relations = []
    homologous_cylinders_partition = defaultdict(list)

    # Convert cylinders into relations.
    for bot, top in cylinders:
        row = [0] * cyl_diag.degree()
        for s in bot:
            row[s] += 1
        for s in top:
            row[s] -= 1
        relations.append(row)
    relations = matrix(QQ, relations)

    # Convert relations into a matrix in reduced row echelon form.
    # Add these relations to `equations`.
    for row in relations.rref():
        for i, one in enumerate(row):
            if one == 1:
                row[i] = 0
                for j in range(len(row)):
                    row[j] = -row[j]
                equations[i] = row
                break
    
    for i, (bot, _) in enumerate(cylinders):
        vec = vector(QQ, [0] * cyl_diag.degree())
        for s in bot:
            vec[s] += 1
        # print(vec)
        for n in equations:
            if vec[n] == 1:
                vec[n] = 0
                vec = vec + equations[n]
        homologous_cylinders_partition[tuple(vec)].append(i)
    output = []
    for v in homologous_cylinders_partition.values():
        if len(v) > 1:
            output.append(v)
    return output

def check_homologous_condition(cyl_diag, partition):
    """Check that homologous cylinders are in the same M-parallel class."""
    classes = find_homologous_cylinders(cyl_diag)
    for homology_class in classes:
        homology_class = frozenset(homology_class)
        for f_set in partition:
            intersection = f_set.intersection(homology_class)
            if intersection != frozenset() and intersection != homology_class:
                return False
    return True

def filter_homologous_condition(cyl_diag, part_list):
    return [part for part in part_list 
                 if check_homologous_condition(cyl_diag, part)]

def check_leaf_condition(cd, partition):
    cylinder_graph = CylinderGraph(cd)
    for leaf, neighbor in cylinder_graph.find_leaves():
        if find_element_in_partition(partition, leaf) == \
           find_element_in_partition(partition, neighbor):
            return False
    return True

def filter_leaf_condition(cd, part_list):
    return [part for part in part_list if check_leaf_condition(cd, part)]

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
