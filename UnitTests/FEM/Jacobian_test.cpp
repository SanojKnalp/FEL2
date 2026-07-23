#include "../catch.hpp"
#include "Grid/COME_Mesh.hpp"
#include "FEM/Quadrature/COME_Quadrature_QGauss.hpp"
#include "FEM/ShapeFunctions/COME_ShapeFunctions_QP.hpp"
#include "DOFs/COME_DofHandler.hpp"
#include "FEM/FEValues/COME_FEValues.hpp"
#include <iostream>
#include <istream>
#include <fstream>

TEST_CASE("Jacobian testing", "[FEM]")
{
    SECTION("Jacobian 3d test")
    {
        ShapeFunctions::FE_QP<3> fe(1);
        Mesh::Mesh<3> mesh;

        std::ifstream input_file("../../../../Meshes/Abaqus3D1Cube.inp");
        mesh.read_abaqus(input_file);

        Quadrature::QGauss<3> quadrature_formula(fe.getPolynomialDegree() + 1);
        FEM::FEValues<3> fe_values(fe, quadrature_formula);

        double check = 0;
        // need to replace with dof_handler
        for (const auto& hex : mesh.getVolumes())
        {
            for (unsigned int i = 0; i < fe_values.get_number_of_quadrature_points(); i++)
            {
                fe_values.reinit(*hex.get());
                check += fe_values.JxW(i);
            }
        }

        REQUIRE(Approx(check) == 1.0);
    }
}
