"""Find R, the max norm over all core curves of all horizontally periodic
cylinders."""

from surface_dynamics import AbelianStratum
from surface_dynamics.databases.flat_surfaces import CylinderDiagrams
from Twist import Twist

def max_length(core_curves_list):
    """Find the max norm in a list of vectors."""
    norm = 0
    for curve in core_curves_list:
        norm = max(norm, curve.norm())
    return norm

def main():
    R = 0
    C = CylinderDiagrams()
    
    # Run through each stratum in genus 3.
    strata = [AbelianStratum(3, 1).components()[0], AbelianStratum(2, 2).components()[1], AbelianStratum(2, 1, 1).components()[0], AbelianStratum(1, 1, 1, 1).components()[0]]
    for H in strata:
        
        # Check all cylinder diagrams with at least 4 cylinders.
        for num_cylinders in range(4, 7):
            try:
                cyl_diag_list = C.get_iterator(H, num_cylinders)
                for cd in cyl_diag_list:
                    tw = Twist(cd)
                    R = max(R, max_length(tw.core_curves))
            
            # Skip if no cylinder diagrams with this many cylinders.
            except ValueError:
                continue
    print(R)

if __name__ == "__main__":
    main()
