# MT2D_Escript_example

MT 2D anisotropy forward modelling example using Esys-Escript 

This example contains a comparison of anisotropy forward modelling results between Esys-Escript and published model from Guo et al.(Modular implementation of magnetotelluric 2d forward modeling with general anisotropy, Computers & Geosciences 118 (2018) 27–38. doi:10.1016/j.cageo.
2018.05.004).

For the model from Figure 8 of Guo et al.(2018), a general anisotropy model setup was used for modelling (see guo2018.png).
The data used for compariosn is generated via finite difference (FD) code from Pek and Verner (Finite-difference modelling of magnetotelluric fields in two-dimensional anisotropic media, Geophysical Journal International (128) (1997) 505–521).

Our modelling data is generated using Esys-Escript (https://github.com/esys-escript/esys-escript.github.io) and the finite element (FE) mesh can be generated using Gmsh (https://gitlab.onelab.info/gmsh/gmsh).
