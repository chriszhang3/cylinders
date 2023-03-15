from collections import defaultdict
import numpy as np
from scipy import optimize
from sage.all import QQ, matrix, vector, span


class Twist:
    """
    A class useful for computing the twist space of a cylinder diagram.
    
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

        # Using the equations, write some saddles as linear combinations of
        # other saddles.
        # Subset of the saddles will be linearly independent elements of
        # H^1(M, Sigma).
        # We will be able to write the elements of twist space as linear
        # combinations of these elements.
        equations = {}

        # Write the core curves in terms of the above basis.
        # self.core_curves is a list of vectors
        # A vector is an instance of
        # `sage.modules.vector_rational_dense.Vector_rational_dense`
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
    
    def find_homologous_cylinders(self):
        """List all classes of homologous cylinders.
        
        Output: List of lists. Each list is an equivalence class of > 1
        homologous cylinders."""

        output = []

        # A dictionary whose keys are elements of homology and the values are
        # cylinders whose core curves are that homology class.
        homologous_cylinders_partition = defaultdict(list)
        for i, coefficients in enumerate(self.core_curves):
            homologous_cylinders_partition[tuple(coefficients)].append(i)

        # Only output if an equivalence class contains more than 1 cylinder.
        for homology_class in homologous_cylinders_partition.values():
            if len(homology_class) > 1:
                output.append(homology_class)
        return output
    
    def twist_dimension(self):
        """Compute the dimension of span γ_i.
        
        Compute the dimension of the space spanned by all of the core curves
        of the horizontal cylinders."""
        return span(self.core_curves).dimension()
    
    def to_numpy(self, index):
        return self.core_curves[index].numpy()

    def class_to_numpy(self, parallel_class):
        return np.stack([self.to_numpy(i) for i in parallel_class], axis=0)

    def ordered_partition(self, partition):
        c0 = list(partition[0])
        c1 = list(partition[1])
        c2 = list(partition[2])
        b = self.to_numpy(c0[0])
        A = np.append(self.class_to_numpy(c1), self.class_to_numpy(c2), axis=0)
        if len(c0) > 1:
            A0 = -self.class_to_numpy(c0[1:])
            A = np.append(A0, A, axis=0)
        
        _, error = optimize.nnls(A.T, b)
        if error < 1.0e-8:
            return True
        return False

    def check_standard_twist_condition(self, partition):
        """ If the partition does not have three equivalence classes, return
        true. Otherwise, check that the equation
        Σa_iα_i = Σb_jβ_j + Σc_kγ_k
        has a solution for a_i,b_j,c_k > 0.
        """

        if len(partition) != 3:
            return True
        partition = list(partition)
        order1 = [partition[0], partition[1], partition[2]]
        order2 = [partition[1], partition[0], partition[2]]
        order3 = [partition[2], partition[1], partition[0]]
        return any([self.ordered_partition(order1),
                    self.ordered_partition(order2),
                    self.ordered_partition(order3)])
