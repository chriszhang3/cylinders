from surface_dynamics import AbelianStratum
from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
from utils import check_pants_relations, check_twist_rel, list_partitions, Twist, check_homologous_cylinders, contains_pants

def find_valid_partitions(cyl_diag_list, num_cylinders, num_classes, check_for_pants):
    output = {}
    for cd in cyl_diag_list:
        if check_for_pants and contains_pants(cd):
            continue
        part = list_partitions(num_cylinders, num_classes)
        part = check_pants_relations(cd, part)
        part = check_twist_rel(Twist(cd), 3, part)
        part = check_homologous_cylinders(cd, part)
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
    H = AbelianStratum(2, 1, 1).components()[0]
    check_cylinder_diagrams_in_stratam(H, 4, True)

if __name__ == '__main__':
    main()