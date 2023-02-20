from sage.all import QQ, matrix, vector


class Twist:
    """
    A class useful for computing the twist space of a translation surface
    
    input: `cd` which is a surface_dynamics.CylinderDiagram
    """
    
    def __init__(self, cd):
        self.cd = cd
        
        # A list of tuples. Each tuple is a cylinder (bot, top), where
        # `bot` is a tuple containing the bottom saddles of the cylinder
        # `top` is a tuple containing the top saddles of the cylinder
        cylinders = cd.cylinders()

        # `relations` is a matrix. Each row has #saddles entries.
        # Each row gives an equation Σ_i row[i]*saddle[i] = 0
        relations = []

        # Using the equations, write some saddles as sums of other saddles.
        # Subset of the saddles will be a basis for homology.
        equations = {}

        # Write the core curves in terms of the above basis.
        self.core_curves = []

        # All relations come from: top of cylinder = bottom of cylinder
        for bot, top in cylinders:
            row = [0] * cd.degree() # degree() gives the number of saddles
            for s in bot:
                row[s] += 1
            for s in top:
                row[s] -= 1
            relations.append(row)
        relations = matrix(QQ, relations)

        # Derives `equations` from `relations`.
        for row in relations.rref(): # rref() is reduced row echelon form
            for i, elt in enumerate(row):
                if elt == 1:
                    eq = vector([-r for r in row])
                    eq[i] = 0
                    equations[i] = eq
                    break

        # Compute self.core_curves
        for i, (bot, _) in enumerate(cylinders):
            vec = vector(QQ, [0] * cd.degree())
            for s in bot:
                vec[s] += 1
            for n in equations:
                if vec[n] == 1:
                    vec[n] = 0
                    vec = vec + equations[n]
            self.core_curves.append(vec)
