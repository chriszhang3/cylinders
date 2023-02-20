from surface_dynamics import AbelianStratum
from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
from lib import list_partitions, filter_homologous_condition, filter_pants_condition, filter_leaf_condition

def find_valid_partitions(cyl_diag_list, num_cylinders, num_classes):
    output = {}
    for cd in cyl_diag_list:
        part = list_partitions(num_cylinders, num_classes)
        part = filter_pants_condition(cd, part)
        part = filter_homologous_condition(cd, part)
        part = filter_leaf_condition(cd, part)
        output[cd] = part
    return output

def dict_extend(dict1, dict2):
    for k, v in dict2.items():
        if k in dict1:
            dict1[k].extend(v)
        else:
            dict1[k] = v

def list_cylinder_classes(H, num_cylinders, num_classes):
    C = CylinderDiagrams()
    cyl_diag_list = C.get_iterator(H, num_cylinders)
    valid = find_valid_partitions(cyl_diag_list, num_cylinders, num_classes)
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