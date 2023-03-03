"""
A script for listing all partitions and filtering out invalid ones.
Currently uses:

filter_homologous_condition,
filter_pants_condition,
filter_leaf_condition.
"""

import os
import sys
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

def name_of_stratum(H):
    name = 'h'
    for z in H.stratum().zeros():
        name = name + str(z)
    return name

def main():
    directory = 'output'
    if not os.path.exists(directory):
       os.makedirs(directory)

    H = AbelianStratum(2, 1, 1).components()[0]
    for i in range(2, 5):     
        with open(f"{directory}/{name_of_stratum(H)}-c5-{i}.txt", 'w') as sys.stdout:
            list_cylinder_classes(H, 5, i)
    
    strata = [AbelianStratum(3, 1).components()[0], AbelianStratum(2, 2).components()[1], AbelianStratum(2, 1, 1).components()[0]]
    for H in strata:
        for i in range(2, 4):
            with open(f"{directory}/{name_of_stratum(H)}-c4-{i}.txt", 'w') as sys.stdout:
                list_cylinder_classes(H, 4, i)

if __name__ == '__main__':
    main()
