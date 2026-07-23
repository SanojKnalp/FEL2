#pragma once
#include "FEM/ShapeFunctions/COME_ShapeFunctions.hpp"
#include "FEM/Quadrature/COME_Quadrature.hpp"
#include "Grid/COME_AbsTopologicalComponent.hpp"
#include "Grid/COME_Node.hpp"
#include "LinearAlgebra/COME_Linalg_dense.hpp"

#include <array>
#include <vector>
#include <iostream>

namespace Mesh {
	template<int dim, int spacedim>
	class AbsTopologicalComponent;
}

namespace FEM
{
	template<int dim, int spacedim = dim>
	class FEValues
	{
	public:
		FEValues(ShapeFunctions::ShapeFunctions<dim, spacedim>& fe, Quadrature::Quadrature<dim>& quadrature, bool update_values = true)
			: fe_(fe)
			, quadrature_(quadrature)
			, update_values_(update_values)
		{
			TensorProductQPointsAndWeights();
		}

		void Jacobian2();

		const double& shape_value(const unsigned int index, const unsigned int q_point) const;
		const double& JxW(const unsigned int q_point) const;
		void reinit(const Mesh::AbsTopologicalComponent<dim,spacedim>& cell);
		std::vector<std::array<double, spacedim>>& get_cell_nodes();

		unsigned int get_number_of_quadrature_points() const;


	private:

		ShapeFunctions::ShapeFunctions<dim, spacedim>& fe_;
		Quadrature::Quadrature<dim>& quadrature_;
		const bool update_values_;
		std::vector<std::array<double, dim>> qPoints_;
		std::vector<double> qWeights_;
		std::vector<std::array<double, spacedim>> cell_nodes;
		std::vector<std::array<unsigned int, dim>> precomputed_indices_;



		double TensorproductShapeFunctionsValue(const unsigned int index, const unsigned int q_point) const;
		std::array<double, dim> TensorproductShapeFunctionGradient(const unsigned int index, const unsigned int q_point) const;
		void TensorProductQPointsAndWeights();
		LinearAlgebra::FullMatrix<double> Jacobian(const unsigned int q_point) const;
		double JacobianDeterminant(const unsigned int q_point) const;

		// precomputes all the indices pairings so that they can be reused later.
		void precompute_indices();

	};

	template <int dim, int spacedim>
	void FEValues<dim, spacedim>::TensorProductQPointsAndWeights()
	{
		std::vector<double> qPoints = quadrature_.getQuadraturePoints();
		std::vector<double> qWeights = quadrature_.getQuadratureWeights();
		int size = qPoints.size();

		if constexpr (dim == 1)
		{
			for (int i = 0; i < size; i++)
			{

				qPoints_.push_back(std::array<double, dim>{ qPoints[i] });
			}
			qWeights_ = qWeights;
		}
		if constexpr (dim == 2)
		{

			for (int i = 0; i < size; i++)
			{
				for (int j = 0; j < size; j++)
				{
					qPoints_.push_back(std::array<double, dim>{ qPoints[i], qPoints[j] });
					qWeights_.push_back(qWeights[i] * qWeights[j]);
				}
			}

		}
		if constexpr (dim == 3)
		{
			for (int i = 0; i < size; i++)
			{
				for (int j = 0; j < size; j++)
				{
					for (int k = 0; k < size; k++)
					{
						qPoints_.push_back( std::array<double, dim>{qPoints[i], qPoints[j], qPoints[k]});
						qWeights_.push_back(qWeights[i] * qWeights[j] * qWeights[k]);
					}
				}
			}

		}
	}

	template <int dim, int spacedim>
	void FEValues<dim, spacedim>::precompute_indices()
	{
		precomputed_indices_.clear();
		for (int i = 0; i < cell_nodes.size(); i++)
		{
			unsigned int temp_index = i;
			const unsigned int n = fe_.getPolynomialDegree() + 1;
			std::array<unsigned int, dim> local_index_nodes;
			for (int i = 0; i < dim; i++)
			{
				unsigned int localIndex = temp_index % n;
				temp_index /= n;
				local_index_nodes[i] = localIndex;
			}
			precomputed_indices_.push_back(local_index_nodes);
		}
	}

	template <int dim, int spacedim>
	void FEValues<dim, spacedim>::reinit(const Mesh::AbsTopologicalComponent<dim, spacedim>& cell)
	{
		cell_nodes.clear();
		for (auto& node : cell.getNodes())
		{
			cell_nodes.push_back(node->getCoordinates());
		}
		precompute_indices();
	}

	template <int dim, int spacedim>
	std::vector<std::array<double, spacedim>>& FEValues<dim, spacedim>::get_cell_nodes()
	{
		return cell_nodes;
	}

	template <int dim, int spacedim>
	void FEValues<dim, spacedim>::Jacobian2()
	{
		std::cout << "jac det" << JacobianDeterminant(0) << std::endl;

	}
}