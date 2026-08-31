#include "COME_FEValues.hpp"

namespace FEM
{
	template <int dim, int spacedim>
	double FEValues<dim, spacedim>::TensorproductShapeFunctionsValue(const unsigned int index, const unsigned int q_point) const
	{
		double value = 1.0;
		const std::array<double, dim>& q = qPoints_[q_point];
		for (int i = 0; i < dim; i++)
		{
			value *= fe_.get_shape_function_value(precomputed_indices_[index][i], q[i]);
		}
		return value;
	}

	template <int dim, int spacedim>
	const double FEValues<dim,spacedim>::shape_value(const unsigned int index, const unsigned int q_point) const
	{
		return TensorproductShapeFunctionsValue(index, q_point);
	}

	template <int dim, int spacedim>
	const LinearAlgebra::Vector<double> FEValues<dim, spacedim>::shape_grad(const unsigned int index, const unsigned int q_point)
	{
		std::array<double, dim> shape_grads = TensorproductShapeFunctionGradient(index, q_point);
		LinearAlgebra::Vector<double> grad(dim);
		for (i = 0; i < dim; i++)
		{
			grad(i) = shape_grads[i];
		}
		return inverse_transpose_jacobians_[index]*grad;
	}


	
	template <int dim, int spacedim>
	std::array<double, dim> FEValues<dim, spacedim>::TensorproductShapeFunctionGradient(const unsigned int index, const unsigned int q_point) const
	{
		std::array<double, dim> shape_grads{};
		const std::array<unsigned int, dim>& indices = precomputed_indices_[index];
		const std::array<double, dim>& q = qPoints_[q_point];


		if constexpr (dim == 1)
		{
			shape_grads[0] = fe_.get_shape_function_grad(indices[0], q[0]);
		}
		else if constexpr (dim == 2)
		{
			shape_grads[0] = fe_.get_shape_function_grad(indices[0], q[0]) * fe_.get_shape_function_value(indices[1], q[1]);
			shape_grads[1] = fe_.get_shape_function_grad(indices[1], q[1]) * fe_.get_shape_function_value(indices[0], q[0]);
		}
		else if constexpr (dim == 3)
		{
			shape_grads[0] = fe_.get_shape_function_grad(indices[0], q[0]) * fe_.get_shape_function_value(indices[1], q[1])* fe_.get_shape_function_value(indices[2], q[2]);
			shape_grads[1] = fe_.get_shape_function_grad(indices[1], q[1]) * fe_.get_shape_function_value(indices[0], q[0])* fe_.get_shape_function_value(indices[2], q[2]);
			shape_grads[2] = fe_.get_shape_function_grad(indices[2], q[2]) * fe_.get_shape_function_value(indices[1], q[1])* fe_.get_shape_function_value(indices[0], q[0]);
		}

		return shape_grads;
	}

	


	template <int dim, int spacedim>
	LinearAlgebra::FullMatrix<double> FEValues<dim, spacedim>::Jacobian(const unsigned int q_point) const
	{
		LinearAlgebra::FullMatrix<double> J(dim, dim);
		// we write a generic jacobian function here that computes the jacobian based on the number of nodes inside of it.
		// this is obviously not ideal for equidistant meshes, but works well enough.
		//i need something here so that i can use the gradients of the shape functions. Luckily, we can construct them easily with tensor products
		for (unsigned int i = 0; i < cell_nodes.size(); i++)
		{
			const auto local_grads = TensorproductShapeFunctionGradient(i, q_point);
			double x = cell_nodes[i][0];
			J(0, 0) += local_grads[0] * x; // dN/dxi  * x

			if constexpr (dim >= 2)
			{
				double y = cell_nodes[i][1];
				J(0, 1) += local_grads[1] * x; // dN/deta * x
				J(1, 0) += local_grads[0] * y; // dN/dxi  * y
				J(1, 1) += local_grads[1] * y; // dN/deta * y
				if constexpr (dim == 3)
				{
					double z = cell_nodes[i][2];
					J(0, 2) += local_grads[2] * x; // dN/dpsi * x
					J(1, 2) += local_grads[2] * y; // dN/dpsi * y
					J(2, 2) += local_grads[2] * z; // dN/dpsi * z
					J(2, 0) += local_grads[0] * z; // dN/dxi  * x
					J(2, 1) += local_grads[1] * z; // dN/dxi  * x
				}
			}
		}
		return J;
	}

	template <int dim, int spacedim>
	double FEValues<dim, spacedim>::JacobianDeterminant(const unsigned int q_point) const
	{
		return jacobians_[q_point].det();
	}

	
	template <int dim, int spacedim>
	const double FEValues<dim, spacedim>::JxW(const unsigned int q_point) const
	{
		return qWeights_[q_point]*JacobianDeterminant(q_point); 

	}

	template <int dim, int spacedim>
	unsigned int FEValues<dim, spacedim>::get_number_of_quadrature_points() const
	{
		return qWeights_.size();
	}
}


template class FEM::FEValues<1,1>;
template class FEM::FEValues<2,2>;
template class FEM::FEValues<3,3>;
template class FEM::FEValues<1,2>;
template class FEM::FEValues<1,3>;
template class FEM::FEValues<2,3>;