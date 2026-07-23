/*
  _____   ____   __  __  ______
 / ____| / __ \ |  \/  ||  ____|
| |     | |  | || \  / || |__
| |     | |  | || |\/| ||  __|
| |____ | |__| || |  | || |____
 \_____| \____/ |_|  |_||______|
*/
#include "LinearAlgebra/COME_Linalg_sparse.hpp"
#include "LinearAlgebra/COME_Linalg_vector.hpp"
#include "LinearSolvers/COME_LinearIterativeSolvers.hpp"
#include "LinearSolvers/COME_KrylovSolvers.hpp"
#include "LinearAlgebra/COME_Linalg_dense.hpp"
#include "Grid/COME_Mesh.hpp"
#include "FEM/Quadrature/COME_Quadrature_QGauss.hpp"
#include "FEM/ShapeFunctions/COME_ShapeFunctions_QP.hpp"
#include "DOFs/COME_DofHandler.hpp"
#include "FEM/FEValues/COME_FEValues.hpp"
#include <iostream>
#include <istream>
#include <fstream>

template<int dim>
class ExampleProblem
{
public:
    ExampleProblem();

    void run();
private:
    void readInMesh();
    void setup_system();
    void assemble_system();
    void solve();

    Mesh::Mesh<dim> mesh;
    ShapeFunctions::FE_QP<dim> fe;
    DOFs::DoFHandler<dim> dof_handler;

    LinearAlgebra::SparseMatrix<double> system_matrix;

    LinearAlgebra::Vector<double> solution;
    LinearAlgebra::Vector<double> system_rhs;
};

template <int dim>
ExampleProblem<dim>::ExampleProblem() 
: fe(1)
, dof_handler(mesh)
{}


template <int dim>
void ExampleProblem<dim>::readInMesh()
{
    std::ifstream input_file("../../../Meshes/Abaqus3D1Cube.inp");
    mesh.read_abaqus(input_file);
}

template <int dim>
void ExampleProblem<dim>::setup_system()
{
    dof_handler.distribute_dofs(fe); 
    unsigned int n_dofs = dof_handler.number_of_dofs();
    std::cout << "n dofs" << n_dofs << std::endl;
    system_matrix.reinit(n_dofs);
    system_rhs.reinit(n_dofs);
    solution.reinit(n_dofs);
    
}

template<int dim>
void ExampleProblem<dim>::assemble_system()
{
    Quadrature::QGauss<dim> quadrature_formula(fe.getPolynomialDegree() + 1);
    FEM::FEValues<dim> fe_values(fe, quadrature_formula);

    const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

    LinearAlgebra::FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    LinearAlgebra::Vector<double>     cell_rhs(dofs_per_cell);

    // need to replace with dof_handler
    for (const auto& hex : mesh.getVolumes())
    {
        cell_matrix.reinit(dofs_per_cell);
        cell_rhs.reinit(dofs_per_cell);
        fe_values.reinit(*hex.get());
        fe_values.Jacobian2();

    }

}

template <int dim>
void ExampleProblem<dim>::solve()
{
    LinearAlgebra::Vector<double> x(system_rhs.size());
    LinearSolvers::KrylovSolvers<double> gmres;
    LinearAlgebra::Vector<double>result = gmres.GmresSolve(system_matrix, system_rhs, x);

}

template<int dim>
void ExampleProblem<dim>::run()
{
    readInMesh();
    setup_system();
    assemble_system();

    //solve();
}

int main() {
    ExampleProblem<3> example;
    example.run();
    return 0;
}
