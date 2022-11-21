import unittest
import networkx as nx
from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
from surface_dynamics import AbelianStratum
from sage.combinat.partition import Partitions
from sage.combinat.set_partition import SetPartitions

def find_generic_pants(cyl_diag):
    """Finds cases when n cylinders are all complete attached to the side of a
    single cylinder.
    
    Input: A cylinder diagram.
    
    Output: A list of lists repsenting homology restrictions."""
    cylinders = cyl_diag.cylinders()
    digraph_data = [[None, None] for _ in range(cyl_diag.degree())]
    for i, (bot, top) in enumerate(cyl_diag.cylinders()):
        for separatrix in bot:
            digraph_data[separatrix][0] = i
        for separatrix in top:
            digraph_data[separatrix][1] = i

    digraph = nx.DiGraph()
    digraph.add_nodes_from(range(len(cylinders)))
    for source, dest in digraph_data:
        digraph.add_edge(source, dest)

    output = []
    for n in digraph:
        successors = list(digraph.successors(n))
        if all([list(digraph.predecessors(suc)) == [n] for suc in successors]):
            successors.append(n)
            output.append(successors)

        predecessors = list(digraph.predecessors(n))
        if all([list(digraph.successors(pre)) == [n] for pre in predecessors]):
            predecessors.append(n)
            # check if this relation already exists
            if not set(successors) == set(predecessors):
                output.append(predecessors)
    return output

def find_element_in_partition(partition, element):
    """Gives the singularity corresponding to the right vertex of the saddle."""
    for i, part in enumerate(partition):
        if element in part:
            return i

def check_conditions(partition, condition_list):
    """Check the partition satisfies a property.
    
    condition is a list of integers. 
    Either one set must contain all integers in the list, or they must be split
    among at least 3 of the sets.
    """
    for condition in condition_list:
        partition_sets = map(
            lambda c: find_element_in_partition(partition, c),
            condition
        )
        if len(set(partition_sets)) == 2:
            return False
    return True

def partitions(n, m, singletons=True):
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

def valid_cylinder_equivalence_classes(cyl_diag, num_classes, free_cylinders = True):
    """Cylinders are numbered. An edge is draw if the cylinders are connected."""
    cylinders = cyl_diag.cylinders()   
    relations = find_generic_pants(cyl_diag)
    part = partitions(len(cylinders), num_classes, free_cylinders)
    part = [p for p in part if check_conditions(p, relations)]
    return part

class Test(unittest.TestCase):

    def test_check_conditions(self):
        self.assertTrue(check_conditions([{1}, {2}, {3}], [[1, 2, 3]]))
        self.assertFalse(check_conditions([{0, 1}, {2, 3}], [[1, 2, 3]]))
    
    def test_partitions(self):
        self.assertEqual(len(partitions(5, 2, singletons=False)), 10)
        self.assertEqual(len(partitions(5, 2, singletons=True)), 15)
        self.assertEqual(len(partitions(6, 2, singletons=False)), 25)
        self.assertEqual(len(partitions(6, 3, singletons=False)), 15)
    
    def test_valid(self):
        C = CylinderDiagrams()
        H = AbelianStratum(3, 1).components()[0]
        valid_classes = [cd for cd in C.get_iterator(H, 4) if valid_cylinder_equivalence_classes(cd, 2, False)]
        self.assertFalse(valid_classes)

        H = AbelianStratum(2, 2).components()[1]
        valid_classes = [cd for cd in C.get_iterator(H, 4) if valid_cylinder_equivalence_classes(cd, 2, False)]
        self.assertEqual(len(valid_classes), 4)

        H = AbelianStratum(2, 1, 1).components()[0]
        valid_classes = [cd for cd in C.get_iterator(H, 4) if valid_cylinder_equivalence_classes(cd, 2, False)]
        self.assertEqual(len(valid_classes), 9)

if __name__ == "__main__":
    unittest.main()
