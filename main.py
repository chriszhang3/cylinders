from surface_dynamics import AbelianStratum
from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
from utils import check_pants_relations, check_twist_rel, list_partitions, Twist, check_homologous_cylinders

def main():
    C = CylinderDiagrams()
    strata = []
    strata.append(AbelianStratum(3, 1).components()[0])
    strata.append(AbelianStratum(2, 2).components()[1])
    for H in strata:
        for cd in C.get_iterator(H, 4):
            part = list_partitions(4, 3)
            part = check_pants_relations(cd, part)
            part = check_twist_rel(Twist(cd), 3, part)
            part = check_homologous_cylinders(cd, part)
            print(cd)
            if part:
                print(part)

if __name__ == '__main__':
    main()