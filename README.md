# 3DLattices

This repository contains a few Mathematica and Matlab functions for the calculation of the homogenized material properties of three-dimensional lattice materials.

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



If you use these routines, or you find them useful for your work, you are kindly requested to cite the above [study](https://www.sciencedirect.com/science/article/pii/S0045782512000941) as 

```
@ARTICLE{Vigliotti201227,
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
