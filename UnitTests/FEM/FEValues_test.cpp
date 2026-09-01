#include "../catch.hpp"
#include "Grid/COME_Mesh.hpp"
#include "FEM/Quadrature/COME_Quadrature_QGauss.hpp"
#include "FEM/ShapeFunctions/COME_ShapeFunctions_QP.hpp"
#include "DOFs/COME_DofHandler.hpp"
#include "FEM/FEValues/COME_FEValues.hpp"
#include <iostream>
#include <istream>
#include <fstream>

TEST_CASE("FEValues testing", "[FEM]")
{
	SECTION("FEValues shape_grad test")
	{
        ShapeFunctions::FE_QP<3> fe(1);
        Mesh::Mesh<3> mesh;

        DOFs::DoFHandler<3> dof_handler(mesh);

        std::ifstream input_file("../../../../Meshes/Abaqus3D1Cube.inp");
        mesh.read_abaqus(input_file);

        Quadrature::QGauss<3> quadrature_formula(fe.getPolynomialDegree() + 1);
        FEM::FEValues<3> fe_values(fe, quadrature_formula);

        dof_handler.distribute_dofs(fe);
        unsigned int n_dofs = dof_handler.number_of_dofs();

        // need to replace with dof_handler
        double check = 0;
        for (const auto& hex : mesh.getVolumes())
        {
            fe_values.reinit(*hex.get());
            for (unsigned int q = 0; q < fe_values.get_number_of_quadrature_points(); q++)
            {   
                for (unsigned int i = 0; i < n_dofs; i++)
                {
                    check += fe_values.shape_grad(i, q).dot(fe_values.shape_grad(i, q))*fe_values.JxW(q);
                }
            }
        }

        REQUIRE(Approx(check) == 8.0/3.0);

	}
}