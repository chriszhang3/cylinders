"""
A script for listing all partitions and filtering out invalid ones.
Currently uses:

filter_homologous_condition,
filter_pants_condition,
filter_leaf_condition.
"""

from surface_dynamics import AbelianStratum
from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
from lib import list_partitions, filter_homologous_condition, \
                filter_pants_condition, filter_leaf_condition

def filter_partitions(cyl_diag_list, num_classes):
    """For each cylinder diagram in cyl_diag_list, list all partitions with 
    num_classes number of classes and filter out the invalid ones.
    
    Stores the output in a dict of
    (cylinder diagram, list of equivalence classes)"""
    output = {}
    for cd in cyl_diag_list:
        part = list_partitions(len(cd.cylinders()), num_classes)
        part = filter_pants_condition(cd, part)
        part = filter_homologous_condition(cd, part)
        part = filter_leaf_condition(cd, part)
        output[cd] = part
    return output

def list_cylinder_classes(H, num_cylinders, num_classes):
    """Runs filter_partitions on every cylinder diagram with num_cylinders
    number of cylinders in a stratum H."""
    C = CylinderDiagrams()
    cyl_diag_list = C.get_iterator(H, num_cylinders)
    valid = filter_partitions(cyl_diag_list, num_classes)
    for i, (k, v) in enumerate(valid.items()):
        print(f"{i+1}. {k}")
        print(v)

def main():
    # H = AbelianStratum(3, 1).components()[0]
    H = AbelianStratum(2, 2).components()[1]
    # H = AbelianStratum(2, 1, 1).components()[0]
    list_cylinder_classes(H, 4, 3)

if __name__ == '__main__':
    main()
