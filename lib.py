from sage.all import Partitions, SetPartitions
from Graph import CylinderGraph
from Twist import Twist

def find_cylinder_in_partition(partition, cylinder):
    """Find the M-parallel class in `partition` than contains `cylinder`.
    
    `partition` is a list of frozen sets.
    `cylinder` is an integer
    Return the index of the element of `partition` that contains `cylinder`."""
    for i, parallel_class in enumerate(partition):
        if cylinder in parallel_class:
            return i
    
    # If cylinder was not found in partition
    return None

def list_partitions(n, m, singletons=True):
    """Return a list of all ways to partition the set [1..n] into m sets.
    
    If singletons==False, do not allow singleton sets.

    Each element of each partition is a frozen set.
    Each partition is an instance of
    `sage.combinat.set_partition.SetPartitions_setparts_with_category.element_class`."""

    partitions = []

    # Partition the integer n into m nonzero integers
    int_parts = Partitions(n, length=m).list()

    # Checks for singleton sets
    if not singletons:
        int_parts = [l for l in int_parts if not any([i == 1 for i in l])]
    
    # Coverts integer partitions into set partitions
    for each_part in int_parts:
        partitions.extend(SetPartitions(range(n), each_part))
    return partitions

def check_pants_condition(partition, pants_list):
    """Check the partition satisfies any homology conditions coming from 
    the pants in pants_list.
    
    `partition` is a partition of the horizontal cylinders into M-parallel
    classes.

    `pants_list` is a list of pants, where each pants is a frozenset containing
    every cylinder in the pants
    
    The homology condition is the following:
    the cylinders in the pants cannot be contained in exactly two distinct sets
    in `partition`.
    """

    for pants in pants_list:
        
        # A list of the cylinder classes that contain a cylinder of pants
        cylinder_classes = map(
            lambda c: find_cylinder_in_partition(partition, c),
            pants
        )
        if len(set(cylinder_classes)) == 2:
            return False
    return True

def filter_pants_condition(cyl_diag, part_list):
    """Filter out the partitions in part_list when check_pants_condition=False.
    """
    cyl_graph = CylinderGraph(cyl_diag)
    pants_list = list(cyl_graph.find_generic_pants())
    return [partition for partition in part_list 
                      if check_pants_condition(partition, pants_list)]

def check_homologous_condition(cyl_diag, partition):
    """Check that homologous cylinders are in the same M-parallel class."""
    tw = Twist(cyl_diag)
    classes = tw.find_homologous_cylinders()
    for homology_class in classes:
        homology_class = frozenset(homology_class)
        for f_set in partition:
            intersection = f_set.intersection(homology_class)
            if intersection != frozenset() and intersection != homology_class:
                return False
    return True

def filter_homologous_condition(cd, part_list):
    """Filter out the partitions in part_list when 
    check_homologous_condition=False."""
    return [part for part in part_list if check_homologous_condition(cd, part)]

def check_leaf_condition(cd, partition):
    """Check for the following condition:
    If a cylinder C is only bordering another cylinder D, then C and D cannot
    be in the same M-parallel class.
    
    Note that this conditions is valid for M a invariant subvariety over Q
    with rank at least 2."""
    cylinder_graph = CylinderGraph(cd)
    for leaf, neighbor in cylinder_graph.find_leaves():
        if find_cylinder_in_partition(partition, leaf) == \
           find_cylinder_in_partition(partition, neighbor):
            return False
    return True

def filter_leaf_condition(cd, part_list):
    """Filter out the partitions in part_list when check_leaf_condition=False.
    """
    return [part for part in part_list if check_leaf_condition(cd, part)]
