from surface_dynamics import AbelianStratum
from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
from utils import list_partitions, Twist, contains_pants, filter_homologous_condition, filter_pants_condition

def find_valid_partitions(cyl_diag_list, num_cylinders, num_classes, check_for_pants):
    output = {}
    for cd in cyl_diag_list:
        if check_for_pants and contains_pants(cd):
            continue
        part = list_partitions(num_cylinders, num_classes)
        part = filter_pants_condition(cd, part)
        part = filter_homologous_condition(cd, part)
        output[cd] = part
    return output

def dict_extend(dict1, dict2):
    for k, v in dict2.items():
        if k in dict1:
            dict1[k].extend(v)
        else:
            dict1[k] = v

def check_cylinder_diagrams_in_stratam(H, num_cylinders, check_for_pants=False):
    C = CylinderDiagrams()
    valid = {}
    for i in range(2, num_cylinders):
        cyl_diag_list = C.get_iterator(H, num_cylinders)
        dict_extend(valid, find_valid_partitions(cyl_diag_list, num_cylinders, i, check_for_pants))
    for k, v in valid.items():
        if v:
            print(k)
            print(v)

def main():
    H = AbelianStratum(3, 1).components()[0]
    # H = AbelianStratum(2, 2).components()[1]
    C = CylinderDiagrams()
    for cd in C.get_iterator(H, 4):
        print("\t\item", cd)
    # H = AbelianStratum(2, 1, 1).components()[0]
    # check_cylinder_diagrams_in_stratam(H, 4, False)

if __name__ == '__main__':
    main()