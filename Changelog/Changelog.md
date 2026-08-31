# Changelog
All the changes i make to the code are documented here.
## 2026 
### August
#### 31.08.2026
Extended the dense linear algebra wrapper to also provide a copy version of the inverse matrices, which is used with constexpr.
Implemented a function that precomputes the jacobians and inverse-transpose jacobians. These are then later used for shape_grad(), which i also implemented.
With all these changes, we should now be able to compute everything we need. Next step is to write tests and test it with a simple hypercube poisson equation in 3D.
### July
#### 23.07.2026
Wrote a function that precomputes the tensor-product indices. Implemented the Tensor product for gradients.
Implemented the Jacobian. Fixed some weird unnecessary forward declaration issues.
### June
#### 26.06.2026
Reworked the mesh readIn for Abaqus to support a lexiographic numbering, which makes it later easier to use tensor products.
### April
#### 28.04.2026 
Worked on connecting cells and the mappings for the computation of the jacobian later.
Add wrapper for determinant function for dense matrices.
For next time: Implement TensorproductShapeFunctionGradient and Jacobian.
#### 10.04.2026
Updated the DoFHandler to now properly distribute the DoFs.
Next Step: We need an active cell iterator, which just loops over all the active cells which we can use later for the computation of the jacobian. Once we have that we are good.
#### 02.04.2026
This is more like a collection of everything since i forgot to update everything:
Mesh readIn seems to work. DofHandler implementation was started.
Added a GMRES, MINRES and BiCGStab solver.
### February
#### 22.02.2026
Reworked current and finished implementation of basic QGauss weights and coordinates in 1d. Also added it's constructor.
#### 17.02.2026
Continued with Mesh ReadIn. Number of Volumes is now correct.
However we need to abstract more such that we can use our addVolume functions and so for.
Also i feel like i'm finally at the point where i don't really need ChatGPT anymore. Shit is getting easier. Also when i do things myself it's getting better cause i know what i want.
Next step: We have created the element and the highest volume. Now we need to create the respective edges, faces and so for.
We will not use non-default constructor here due to the nature of how different meshes work.
#### 13.02.2026
Added changes to github.
#### 12.02.2026
Worked again on a bit of mesh stuff, especially the connection between AbsTopological Component and Element. Next step is finishing the element creation
#### 09.02.2026 
Worked a lot on mesh read in for abaqus and in general mesh functions. 
#### 07.02.2026
Added the functions inverse2x2() and inverse 3x3 that provide analytical formulations for the respective inverse matrices.
Also added tests for both.
Started working on reading in some meshes. Not the best experience so far.
#### 03.02.2026 
Added the first couple QGauss weights.
Added a bit of ascii art to the main for no reason :).
Added the Christoffelsymbols for spherical coordinates
#### 02.02.2026 
Started working on the abstract class Quadrature<dim> and a derived class QGauss<dim> which then later will be used for Quadrature.
Added a Test for the Jacobi Solver.
Added a Test for get_shape_function_value().
#### 01.02.2026
Finished the get_QP_grad_grad_degree_i() implementation
Finished the get_QP_third_derivative_degree_i() implementation
Finished the get_QP_fourth_derivative_degree_i() implementation