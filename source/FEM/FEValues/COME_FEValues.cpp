#include "COME_FEValues.hpp"

namespace FEM
{
	template <int dim, int spacedim>
	double FEValues<dim, spacedim>::TensorproductShapeFunctionsValue(const unsigned int index, std::array<double, dim>& localPoints) const
	{
		double value = 1.0;
		unsigned int temp_index = index;
		const unsigned int n = fe_.getPolynomialDegree() + 1;
		for (int i = 0; i < dim; i++)
		{
			int localIndex = temp_index % n;
			temp_index /= n;
			value *= fe_.get_shape_function_value(localIndex, localPoints[i]);
		}
		return value;
	}

	template <int dim, int spacedim>
	const double& FEValues<dim,spacedim>::shape_value(const unsigned int index, const unsigned int q_point) const
	{
		return TensorproductShapeFunctionsValue(index, qPoints_[q_point]);
	}


	template <int dim, int spacedim>
	LinearAlgebra::FullMatrix<double> FEValues<dim, spacedim>::Jacobian(const unsigned int q_point)
	{
		FullMatrix<double> J(dim, dim);
		// we write a generic jacobian function here that computes the jacobian based on the number of nodes inside of it.
		// this is obviously not ideal for equidistant meshes, but works well enough.
		//i need something here so that i can use the gradients of the shape functions. Luckily, we can construct them easily via tensor products.
		for (unsigned int i = 0; i < cell_nodes.size(); i++)
		{

		}
		return J;
	}

	template <int dim, int spacedim>
	double FEValues<dim, spacedim>::JacobianDeterminant(const unsigned int q_point)
	{
		return Jacobian(q_point).det();
	}

	
	template <int dim, int spacedim>
	const double& FEValues<dim, spacedim>::JxW(const unsigned int q_point) const
	{
		return qWeights_[q_point]*JacobianDeterminant(q_point); // we still need the jacobian

	}
}


template class FEM::FEValues<1,1>;
template class FEM::FEValues<2,2>;
template class FEM::FEValues<3,3>;
template class FEM::FEValues<1,2>;
template class FEM::FEValues<1,3>;
template class FEM::FEValues<2,3>;