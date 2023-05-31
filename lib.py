from sage.all import Partitions, SetPartitions
from Graph import CylinderGraph
from Twist import Twist


### Graph Functions
# TODO: Add unittests for this
def find_generalized_pants(digraph):
    """Finds cases when n cylinders are all only attached to the side
    of a single cylinder. This is like a generalized version of a
    topological pair of pants allows n pant legs.
    
    Input: A cylinder diagram.
    
    Output: A set of pants. Each pants is a frozenset consisting of the 
    cylinders in the pants.
    
    Note: For Python 3.7 and higher, frozensets should maintain insertion
    order. Thus, the first element of each pants is the waist curve and the
    rest are the pant legs. However, we do not use this in our code."""

    pants_set = set()
    for n in digraph:
        # If every cylinder above `C` is only adjacent to `C` along its
        # bottom, this is a generalized pants.
        neighbors_out = list(digraph.neighbors_out(n))
        if all([list(digraph.neighbors_in(suc)) == [n] 
                for suc in neighbors_out]):
            
            pants = frozenset([n] + neighbors_out)
            pants_set.add(pants)

        # If every cylinder below `C` is only adjacent to `C` along its
        # top, this is a generalized pants.
        neighbors_in = list(digraph.neighbors_in(n))
        if all([list(digraph.neighbors_out(pre)) == [n] 
                for pre in neighbors_in]):

            pants = frozenset([n] + neighbors_in)
            pants_set.add(pants)
    return pants_set

def find_leaves(digraph):
    """Return the tuples (leaf, neighbor) for all leaves of `digraph`.
    
    A leaf is a vertex such that it only has one neighbor, where we 
    count a vertex as a neighbor if there is either an edge coming from it
    or an edge going to it.
    
    For a given leaf, `neighbor` is it's unique neighbor."""
    leaf_neighbers = []
    for n in digraph:
        neighbors = set(digraph.neighbors_out(n)) | \
                    set(digraph.neighbors_in(n))
        if len(neighbors) == 1:
            leaf_neighbers.append((n, next(iter(neighbors))))
    return leaf_neighbers



### Partition related
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

### Pants condition
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
    cyl_graph = CylinderGraph(cyl_diag).digraph
    pants_list = list(find_generalized_pants(cyl_graph))
    return [partition for partition in part_list 
                      if check_pants_condition(partition, pants_list)]

### Homologous condition
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

### Leaf condition
def is_simple(cd, index):
    """Check if cylinder number `index` is simple."""
    cylinder = cd.cylinders()[index]
    if len(cylinder[0]) == 1 and len(cylinder[1]) == 1:
        return True
    return False

def check_leaf_condition(cd, partition):
    """Check for the following condition:
    If a cylinder C is only bordering another cylinder D, then C and D cannot
    be in the same M-parallel class.
    
    Note that there are some assumptions for this condition.
    Refer to the paper for these assumptions."""
    cyl_graph = CylinderGraph(cd).digraph
    for leaf, neighbor in find_leaves(cyl_graph):
        if is_simple(cd, leaf):
            if find_cylinder_in_partition(partition, leaf) == \
            find_cylinder_in_partition(partition, neighbor):
                return False
    return True

def filter_leaf_condition(cd, part_list):
    """Filter out the partitions in part_list when check_leaf_condition=False.
    """
    return [part for part in part_list if check_leaf_condition(cd, part)]

### Standard Twist Condition

def filter_standard_twist_condition(cd, part_list):
    tw = Twist(cd)
    return [part for part in part_list 
                 if tw.check_standard_twist_condition(part)]
