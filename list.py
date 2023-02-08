from surface_dynamics import AbelianStratum
from surface_dynamics.databases.flat_surfaces import CylinderDiagrams

def main():
    C = CylinderDiagrams()
    # H = AbelianStratum(3, 1).components()[0]
    # H = AbelianStratum(2, 2).components()[1]
    H = AbelianStratum(2, 1, 1).components()[0]
    for cd in C.get_iterator(H, 5):
        print(cd)

if __name__ == '__main__':
    main()