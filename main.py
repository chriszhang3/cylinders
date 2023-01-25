from surface_dynamics import AbelianStratum
from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
from utils import check_pants_relations, check_twist_rel, list_partitions, Twist, check_homologous_cylinders

def find_valid_partitions(H, num_cylinders, num_classes):
    C = CylinderDiagrams()
    output = {}
    for cd in C.get_iterator(H, num_cylinders):
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

def main():
    C = CylinderDiagrams()
    H = AbelianStratum(2, 1, 1).components()[0]

    valid = {}
    for i in range(2, 5):
        dict_extend(valid, find_valid_partitions(H, 5, i))
    for k, v in valid.items():
        if v:
            print(k)
            print(v)

if __name__ == '__main__':
    main()