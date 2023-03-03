# Code related to my forthcoming thesis work

This code rules out certain cylinder diagrams and M-parallel classes using [Sage](https://www.sagemath.org/) and the [surface-dynamics](https://flatsurf.github.io/surface-dynamics/index.html) package.

Definition: a **cylinder digram** is a way to represent a horizontally periodic translation surface. In this codebase, it is an instance of the class `surface_dynamics.flat_surfaces.separatrix_diagram.CylinderDiagram` see the [documentation](https://flatsurf.github.io/surface-dynamics/surface_topology.html#surface_dynamics.flat_surfaces.separatrix_diagram.CylinderDiagram) for more details.

For a fixed invariant subvariety $\mathcal M$, the horizontal cylinders can be partitioned into $\mathcal M$-parallel classes. There are certain restrictions on which cylinders can be in the same $\mathcal M$-parallel classes. This codebase automatically checks some of these conditions.

## Dependencies

- [Instructions](https://doc.sagemath.org/html/en/installation/index.html) for installing Sage.
- After Sage is installed, you can install the surface-dynamics package with  
`sage -pip install surface-dynamics`

## Run
Run the program using  
`sage -python main.py`  
the output files will be stored in a chosen directory. By default the directory is `output`.