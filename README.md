# 3DLattices

This repository contains a few Mathematica and Matlab functions for the calculation of the homogenized material properties of three-dimensional lattice materials.

The methodology employed here takes as input the model of the uint cell in terms of coordinates of the nodes, topology of the connecting beams and faces, as well as the periodic directions of the lattice, which are the directions along which the unit cell is meant to be replicated, and produces as outputs the stiffness matrix of the homegenized resulting lattice, in addition to show the deformed configuration of the lattice for pure deformation, in the frame of reference of the undeformed unit cell.

More details on the homogenization technique can be found in [Stiffness and strength of tridimensional periodic lattices](https://www.sciencedirect.com/science/article/pii/S0045782512000941)

## Mathematica scripts

the file ```Calc3DMatProp.m``` contains the functions needed to calculate the stiffness matrix of the RVE and to calculate the homegenized stiffness of the lattice. 
The file listed below contain the node coordinate and the connectivity of the lattices analized in the above mentioned  [paper](https://www.sciencedirect.com/science/article/pii/S0045782512000941)
```
- CubicLatticeBCC.nb
- CubicLatticeFCC.nb
- CubicLattice.nb
- Cuboctahedron.nb
- RegularOctet.nb
```

## Matlab scripts

the file ```Find3DMatProp.m``` contatins the function needed to calculate the stiffness matrix of the lattices, the file ```plot3DModel.m``` does the plotting. 
The file listed below contain the node coordinate and the connectivity of the lattices analized in the above mentioned [paper](https://www.sciencedirect.com/science/article/pii/S0045782512000941)
```
- KCubic.m
- KCubOctahedron.m
- KGreatRhombicuboctahedron_a.m
- KGreatRhombicuboctahedron_b.m
- KGreatRhombicuboctahedron_c.m
- KOctet.m
- KSmallRhombicosidodecahedron_a.m
- KSmallRhombicosidodecahedron_b.m
- KSmallRhombicosidodecahedron.m
- KSmallRhombicuboctahedron_a.m
- KSmallRhombicuboctahedron_b.m
- KSmallRhombicuboctahedron_c.m
- KTruncatedCube.m
- KTruncatedDodecahedron_a.m
- KTruncatedOctahedron_a.m
```

as an example the following the figures show the deformed lattice's unit cells and the and the outputs produced for three different choices of the periodic directions for the Small Rhombicuboctahedron lattice. The figures are produced by the ```KSmallRhombicuboctahedron_a.m```, ```KSmallRhombicuboctahedron_b.m``` and ```KSmallRhombicuboctahedron_c.m``` scripts:

```
KGreatRhombicuboctahedron_a
the material volume of the unit cell is : 0.0151
the volume of the unit cell is          : 56.1127
the relative density of the lattice is  : 0.0003
the lattice stiffnes matrix is : 
   1.0e+06*

    7.1256    7.1245    7.1245   -0.0000   -0.0000    0.0000
    7.1245    7.1256    7.1245   -0.0000   -0.0000    0.0000
    7.1245    7.1245    7.1256   -0.0000   -0.0000   -0.0000
   -0.0000   -0.0000    0.0000    0.0005    0.0000   -0.0000
   -0.0000    0.0000   -0.0000    0.0000    0.0005   -0.0000
    0.0000    0.0000    0.0000   -0.0000   -0.0000    0.0005
```
![KGreatRhombicuboctahedron_a](https://github.com/avigliotti/3DLattices/blob/master/figures/KGreatRhombicuboctahedron_a.png)

```
KGreatRhombicuboctahedron_b
the material volume of the unit cell is : 0.0151
the volume of the unit cell is          : 56.2843
the relative density of the lattice is  : 0.0003
the lattice stiffnes matrix is : 
   1.0e+06 *
    7.5557    7.5552    7.5552    0.0000    0.0000    0.0000
    7.5552    7.5557    7.5552    0.0000   -0.0000    0.0000
    7.5552    7.5552    7.5557    0.0000    0.0000   -0.0000
   -0.0000   -0.0000    0.0000    0.0005    0.0000   -0.0000
   -0.0000    0.0000   -0.0000    0.0000    0.0005   -0.0000
    0.0000    0.0000   -0.0000    0.0000   -0.0000    0.0005
```
![KGreatRhombicuboctahedron_b](https://github.com/avigliotti/3DLattices/blob/master/figures/KGreatRhombicuboctahedron_b.png)

```
KGreatRhombicuboctahedron_c
the material volume of the unit cell is : 0.0151
the volume of the unit cell is          : 60.8198
the relative density of the lattice is  : 0.0002
the lattice stiffnes matrix is : 
   1.0e+06*

    7.5783    7.5772    7.5772    0.0000   -0.0000   -0.0000
    7.5772    7.5783    7.5772    0.0000   -0.0000   -0.0000
    7.5772    7.5772    7.5783   -0.0000    0.0000    0.0000
   -0.0000    0.0000   -0.0000    0.0004    0.0000   -0.0000
    0.0000    0.0000    0.0000    0.0000    0.0004    0.0000
   -0.0000   -0.0000   -0.0000   -0.0000    0.0000    0.0004
```   
![KGreatRhombicuboctahedron_c](https://github.com/avigliotti/3DLattices/blob/master/figures/KGreatRhombicuboctahedron_c.png)


If you use these routines, or you find them useful for your work, consider citing the above [study](https://www.sciencedirect.com/science/article/pii/S0045782512000941) as 

```
@ARTICLE{CMAME201227,
   author={Vigliotti, A. and Pasini, D.}, 
   title={Stiffness and strength of tridimensional periodic lattices},
   journal={Computer Methods in Applied Mechanics and Engineering},
   year={2012},
   volume={229-232},
   pages={27-43},
   doi={10.1016/j.cma.2012.03.018},
   url={https://www.sciencedirect.com/science/article/pii/S0045782512000941?via%3Dihub},
   document_type={Article},
}

```
