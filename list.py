"""A script for listing the contents of
https://flatsurf.github.io/surface-dynamics/database.html"""

from surface_dynamics import AbelianStratum
from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
from lib import Twist


def main():
    C = CylinderDiagrams()
    # H = AbelianStratum(3, 1).components()[0]
    # H = AbelianStratum(2, 2).components()[1]
    # H = AbelianStratum(2, 1, 1).components()[0]
    H = AbelianStratum(1, 1, 1, 1).components()[0]
    for i, cd in enumerate(C.get_iterator(H, 5)):
        tw = Twist(cd)
        # print(f"{i+1}. {cd} {tw.twist_dimension()}")
        print(f"{i+1}. {cd}")

if __name__ == '__main__':
    main()