from surface_dynamics import AbelianStratum
from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
from utils import check_pants_relations, check_twist_rel, list_partitions, Twist, check_homologous_cylinders

def find_valid_partitions(H, num_cylinders, num_classes):
    C = CylinderDiagrams()
    output = {}
    for cd in C.get_iterator(H, num_cylinders):
        part = list_partitions(num_cylinders, num_classes)
        part = check_pants_relations(cd, part)
        part = check_twist_rel(Twist(cd), num_classes, part)
        part = check_homologous_cylinders(cd, part)
        output[cd] = part
    return output
    

def main():
    C = CylinderDiagrams()
    strata = []
    strata.append(AbelianStratum(3, 1).components()[0])
    # strata.append(AbelianStratum(2, 2).components()[1])
    for H in strata:
        valid = find_valid_partitions(H, 4, 3)
        for k, v in valid.items():
            print(k)
            print(v)

if __name__ == '__main__':
    main()