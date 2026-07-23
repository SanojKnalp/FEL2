#include "COME_ShapeFunctions_QP.hpp"
#include <cmath>

namespace ShapeFunctions
{
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_shape_function_value(const unsigned int index, const double qPointValue) const
	{
		switch (this->getPolynomialDegree())
		{
		case 1: return get_QP_degree_1(index, qPointValue);
		case 2: return get_QP_degree_2(index, qPointValue);
		case 3: return get_QP_degree_3(index, qPointValue);
		case 4: return get_QP_degree_4(index, qPointValue);
		case 5: return get_QP_degree_5(index, qPointValue);
		case 6: return get_QP_degree_6(index, qPointValue);
		default:
			throw std::out_of_range("A shape function with degree seven or higher was selected. This is not implemented as of right now.");
		}
	}

	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_shape_function_grad(const unsigned int index, const double qPointValue) const
	{
		switch (this->getPolynomialDegree())
		{
		case 1: return get_QP_grad_degree_1(index, qPointValue);
		case 2: return get_QP_grad_degree_2(index, qPointValue);
		case 3: return get_QP_grad_degree_3(index, qPointValue);
		case 4: return get_QP_grad_degree_4(index, qPointValue);
		case 5: return get_QP_grad_degree_5(index, qPointValue);
		case 6: return get_QP_grad_degree_6(index, qPointValue);
		default:
			throw std::out_of_range("A shape function (gradient) with degree seven or higher was selected. This is not implemented as of right now.");
		}
	}

	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_shape_function_grad_grad(const unsigned int index, const double qPointValue) const
	{
		switch (this->getPolynomialDegree())
		{
		case 1: return get_QP_grad_grad_degree_1(index, qPointValue);
		case 2: return get_QP_grad_grad_degree_2(index, qPointValue);
		case 3: return get_QP_grad_grad_degree_3(index, qPointValue);
		case 4: return get_QP_grad_grad_degree_4(index, qPointValue);
		case 5: return get_QP_grad_grad_degree_5(index, qPointValue);
		case 6: return get_QP_grad_grad_degree_6(index, qPointValue);
		default:
			throw std::out_of_range("A shape function second derivative with degree seven or higher was selected. This is not implemented as of right now.");
		}
	}

	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_shape_function_third_derivative(const unsigned int index, const double qPointValue) const
	{
		switch (this->getPolynomialDegree())
		{
		case 1: return get_QP_third_derivative_degree_1(index, qPointValue);
		case 2: return get_QP_third_derivative_degree_2(index, qPointValue);
		case 3: return get_QP_third_derivative_degree_3(index, qPointValue);
		case 4: return get_QP_third_derivative_degree_4(index, qPointValue);
		case 5: return get_QP_third_derivative_degree_5(index, qPointValue);
		case 6: return get_QP_third_derivative_degree_6(index, qPointValue);
		default:
			throw std::out_of_range("A shape function third derivative with degree seven or higher was selected. This is not implemented as of right now.");
		}
	}

	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_shape_function_fourth_derivative(const unsigned int index, const double qPointValue) const
	{
		switch (this->getPolynomialDegree())
		{
		case 1: return get_QP_fourth_derivative_degree_1(index, qPointValue);
		case 2: return get_QP_fourth_derivative_degree_2(index, qPointValue);
		case 3: return get_QP_fourth_derivative_degree_3(index, qPointValue);
		case 4: return get_QP_fourth_derivative_degree_4(index, qPointValue);
		case 5: return get_QP_fourth_derivative_degree_5(index, qPointValue);
		case 6: return get_QP_fourth_derivative_degree_6(index, qPointValue);
		default:
			throw std::out_of_range("A shape function fourth derivative with degree seven or higher was selected. This is not implemented as of right now.");
		}
	}
	//-------------------------------------------
	//------- FE shape values -------------------
	//-------------------------------------------

	// ---------------- degree 1 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_degree_1(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;

		switch (index)
		{
		case 0: return 0.5 * (1.0 - x);
		case 1: return 0.5 * (1.0 + x);
		default:
			throw std::out_of_range("Invalid shape function index for Q1.");
		}
	}

	// ---------------- degree 2 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_degree_2(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		const double x2 = x * x;

		switch (index)
		{
		case 0: return 0.5 * x * (x - 1.0);
		case 1: return 1.0 - x2;
		case 2: return 0.5 * x * (x + 1.0);
		default:
			throw std::out_of_range("Invalid shape function index for Q2.");
		}
	}

	// ---------------- degree 3 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_degree_3(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		const double x2 = x * x;
		const double x3 = x2 * x;

		switch (index)
		{
		case 0:
			return 0.5625 * (x2 - x3) + 0.0625 * (x - 1.0);

		case 1:
			return 1.6875 * (x3 - x) + 0.5625 * (1.0 - x2);

		case 2:
			return 1.6875 * (x - x3) + 0.5625 * (x2 - 1.0);

		case 3:
			return 0.5625 * (x3 + x2) - 0.0625 * (x + 1.0);

		default:
			throw std::out_of_range("Invalid shape function index for Q3.");
		}
	}

	// ---------------- degree 4 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_degree_4(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		const double x2 = x * x;
		const double x3 = x2 * x;
		const double x4 = x2 * x2;

		switch (index)
		{
		case 0:
			return (2.0 / 3.0) * (x4 - x3) + (1.0 / 6.0) * (x - x2);

		case 1:
			return (8.0 / 3.0) * (x2 - x4) + (4.0 / 3.0) * (x3 - x);

		case 2:
			return 4.0 * x4 - 5.0 * x2 + 1.0;

		case 3:
			return (8.0 / 3.0) * (x2 - x4) + (4.0 / 3.0) * (x - x3);

		case 4:
			return (2.0 / 3.0) * (x4 + x3) - (1.0 / 6.0) * (x + x2);

		default:
			throw std::out_of_range("Invalid shape function index for Q4.");
		}
	}

	// ---------------- degree 5 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_degree_5(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		const double x2 = x * x;
		const double x3 = x2 * x;
		const double x4 = x2 * x2;
		const double x5 = x4 * x;

		switch (index)
		{
		case 0:
			return 0.8138020833333333 * (x4 - x5)
				+ 0.3255208333333333 * (x3 - x2)
				+ 0.01171875 * (1.0 - x);

		case 1:
			return  4.069010416666667 * x5
				- 2.44140625 * x4
				- 4.231770833333333 * x3
				+ 2.5390625 * x2
				+ 0.162760416666667 * x
				- 0.09765625;

		case 2:
			return -8.138020833333333 * x5
				+ 1.627604166666667 * x4
				+ 11.067708333333333 * x3
				- 2.213541666666667 * x2
				- 2.9296875 * x
				+ 0.5859375;

		case 3:
			return  8.138020833333333 * x5
				+ 1.627604166666667 * x4
				- 11.067708333333333 * x3
				- 2.213541666666667 * x2
				+ 2.9296875 * x
				+ 0.5859375;

		case 4:
			return -4.069010416666667 * x5
				- 2.44140625 * x4
				+ 4.231770833333333 * x3
				+ 2.5390625 * x2
				- 0.162760416666667 * x
				- 0.09765625;

		case 5:
			return 0.8138020833333333 * (x4 + x5)
				- 0.3255208333333333 * (x3 + x2)
				+ 0.01171875 * (1.0 + x);

		default:
			throw std::out_of_range("Invalid shape function index for Q5.");
		}
	}

	// ---------------- degree 6 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_degree_6(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		const double x2 = x * x;
		const double x3 = x2 * x;
		const double x4 = x2 * x2;
		const double x5 = x4 * x;
		const double x6 = x3 * x3;
		switch (index)
		{
		case 0: return 	1.0125 * (x6 - x5)
			+ 0.5625 * (x3 - x4)
			+ 0.05 * (x2 - x);

		case 1:return -6.075 * x6
			+ 4.05 * x5
			+ 6.75 * x4
			- 4.5 * x3
			- 0.675 * x2
			+ 0.45 * x;

		case 2:return 15.1875 * x6
			- 5.0625 * x5
			- 21.9375 * x4
			+ 7.3125 * x3
			+ 6.75 * x2
			- 2.25 * x;

		case 3:return -20.25 * x6
			+ 31.5 * x4
			- 12.25 * x2
			+ 1;
		case 4:return 15.1875 * x6
			+ 5.0625 * x5
			- 21.9375 * x4
			- 7.3125 * x3
			+ 6.75 * x2
			+ 2.25 * x;
		case 5:return -6.075 * x6
			- 4.05 * x5
			+ 6.75 * x4
			+ 4.5 * x3
			- 0.675 * x2
			- 0.45 * x;
		case 6:return 1.0125 * (x6 + x5)
			- 0.5625 * (x3 + x4)
			+ 0.05 * (x2 + x);
		default:
			throw std::out_of_range("Invalid shape function index for Q6.");
		}
	}



	//-------------------------------------------
	//------- FE shape grads --------------------
	//-------------------------------------------

	// ---------------- degree 1 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_grad_degree_1(const unsigned int index,
		const double qPointValue) const
	{
		switch (index)
		{
		case 0: return -0.5;
		case 1: return  0.5 ;
		default:
			throw std::out_of_range("Invalid shape function index for Q1.");
		}
	}

	// ---------------- degree 2 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_grad_degree_2(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		switch (index)
		{
		case 0: return x-0.5;
		case 1: return  -2.0*x;
		case 2: return x+0.5;
		default:
			throw std::out_of_range("Invalid shape function index for Q2.");
		}
	}


	// ---------------- degree 3 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_grad_degree_3(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		const double x2 = x*x;
		switch (index)
		{
		case 0: return 1.125*x-1.6875*x2+0.0625;
		case 1: return 5.0625*x2-1.125*x-1.6875 ;
		case 2: return -5.0625 * x2 + 1.125 * x + 1.6875;
		case 3: return 1.6875*x2+1.125*x-0.0625;
		default:
			throw std::out_of_range("Invalid shape function index for Q3.");
		}
	}

	// ---------------- degree 4 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_grad_degree_4(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		const double x2 = x * x;
		const double x3 = x2 * x;
		switch (index)
		{
		case 0: return 8.0 / 3.0 * x3 - 2.0 * x2 - 1.0 / 3.0 * x;
		case 1: return 16.0 / 3.0 * x - 32.0 / 3.0 * x3 + 4.0 * x2 - 4.0 / 3.0;
		case 2:	return 16.0 * x3 - 10.0 * x;
		case 3:	return 16.0 / 3.0 * x - 32.0 / 3.0 * x3 + 4.0 / 3.0 - 4.0 * x2;
		case 4:	return 8.0 / 3.0 * x3 + 2.0 * x2 - 1.0 / 6.0 - 1.0 / 3.0 * x;
		default:
			throw std::out_of_range("Invalid shape function index for Q4.");
		}

	}

	// ---------------- degree 5 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_grad_degree_5(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		const double x2 = x * x;
		const double x3 = x2 * x;
		const double x4 = x2 * x2;
		switch (index)
		{
		case 0:	return 625.0 / 192.0 * x3 - 3125.0 / 768.0 * x4 + 125.0 / 128.0 * x2 - 125.0 / 192.0 * x - 0.01171875;
		case 1:	return 15625.0 / 768.0 * x4 - 625.0 / 64.0 * x3 - 1625.0 / 128.0 * x2 + 5.078125 * x + 0.162760416666667;
		case 2:	return -15625.0 / 384.0 * x4 + 625.0 / 96.0 * x3 + 2125.0 / 64.0 * x2 - 425.0 / 96.0 * x - 2.9296875;
		case 3:	return 15625.0 / 384.0 * x4 + 625.0 / 96.0 * x3 - 2125.0 / 64.0 * x2 - 425.0 / 96.0 * x + 2.9296875;
		case 4:	return -15625.0 / 768.0 * x4 - 625.0 / 64.0 * x3 + 1625.0 / 128.0 * x2 + 5.078125 * x - 0.162760416666667;
		case 5:	return 625.0 / 192.0 * x3 + 3125.0 / 768.0 * x4 - 125.0 / 128.0 * x2 - 125.0 / 192.0 * x + 0.01171875;
		default:
			throw std::out_of_range("Invalid shape function index for Q5.");
		}

	}

	// ---------------- degree 6 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_grad_degree_6(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		const double x2 = x * x;
		const double x3 = x2 * x;
		const double x4 = x2 * x2;
		const double x5 = x3 * x2;
		switch (index)
		{
		case 0:	return 6.075 * x5 - 5.0625 * x4 + 1.6875 * x2 - 2.25 * x3 + 0.1 * x - 0.05;
		case 1:	return -36.45 * x5 + 20.25 * x4 + 27.0 * x3 - 13.5 * x2 - 1.35 * x + 0.45;
		case 2:	return 191.125 * x5 - 25.3125 * x4 - 87.75 * x3 + 21.9375 * x2 + 13.5 * x - 2.25;
		case 3:	return -121.5 * x5 + 126.0 * x3 - 24.5 * x;
		case 4:	return 191.125 * x5 + 25.3125 * x4 - 87.75 * x3 - 21.9375 * x2 + 13.5 * x + 2.25;
		case 5:	return -36.45 * x5 - 20.25 * x4 + 27.0 * x3 + 13.5 * x2 - 1.35 * x - 0.45;
		case 6:	return 6.075 * x5 + 5.0625 * x4 - 1.6875 * x2 - 2.25 * x3 + 0.1 * x + 0.05;
		default:
			throw std::out_of_range("Invalid shape function index for Q6.");
		}

	}



	//-------------------------------------------
	//------- FE shape grad grads ---------------
	//-------------------------------------------

	// ---------------- degree 1 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_grad_grad_degree_1(const unsigned int index,
		const double qPointValue) const
	{
		if (index >1)
		{
			throw std::out_of_range("Invalid shape function index for Q1.");
		}
		return 0.0;
	}

	// ---------------- degree 2 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_grad_grad_degree_2(const unsigned int index,
		const double qPointValue) const
	{
		switch (index)
		{
		case 0: return 1.0;
		case 1: return -2.0;
		case 2:	return 1.0;
		default:
			throw std::out_of_range("Invalid shape function index for Q2.");

		}

	}

	// ---------------- degree 3 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_grad_grad_degree_3(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		switch (index)
		{
		case 0: return 1.125 - 3.375 * x;
		case 1: return 10.125 * x - 1.125;
		case 2: return -10.125 * x + 1.125;
		case 3: return 3.375 * x + 1.125;
		default:
			throw std::out_of_range("Invalid shape function index for Q3.");
		}

	}

	// ---------------- degree 4 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_grad_grad_degree_4(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		const double x2 = x * x;
		switch (index)
		{
		case 0: return 8.0 * x2 - 4.0 * x - 1.0 / 3.0;
		case 1: return 16.0 / 3.0 - 32.0 * x2 + 8.0 * x;
		case 2:	return 48.0 * x2 - 10.0;
		case 3:	return 16.0 / 3.0 - 32.0 * x2 - 8.0 * x;
		case 4:	return 8.0 * x2 + 4.0 * x - 1.0 / 3.0;
		default:
			throw std::out_of_range("Invalid shape function index for Q4.");
		}

	}

	// ---------------- degree 5 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_grad_grad_degree_5(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		const double x2 = x * x;
		const double x3 = x * x2;
		switch (index)
		{
		case 0:	return 625.0 / 64.0 * x2 - 3125.0 / 192.0 * x3 + 125.0 / 64.0 * x - 125.0 / 192.0;
		case 1:	return 15625.0 / 192.0 * x3 - 1875.0 / 64.0 * x2 - 1625.0 / 64.0 * x + 5.078125;
		case 2:	return -15625.0 / 96.0 * x3 + 625.0 / 32.0 * x2 + 2125.0 / 32.0 * x - 425.0 / 96.0;
		case 3:	return 15625.0 / 96.0 * x3 + 625.0 / 32.0 * x2 - 2125.0 / 32.0 * x - 425.0 / 96.0;
		case 4:	return -15625.0 / 192.0 * x3 - 1875.0 / 64.0 * x2 + 1625.0 / 64.0 * x + 5.078125;
		case 5:	return 625.0 / 64.0 * x2 + 3125.0 / 192.0 * x3 - 125.0 / 64.0 * x - 125.0 / 192.0;
		default:
			throw std::out_of_range("Invalid shape function index for Q5.");
		}

	}

	// ---------------- degree 6 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_grad_grad_degree_6(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		const double x2 = x * x;
		const double x3 = x2 * x;
		const double x4 = x2 * x2;
		switch (index)
		{
		case 0:	return 30.375 * x4 - 20.25 * x3 + 3.375 * x - 6.75 * x2 + 0.1;
		case 1:	return -182.25 * x4 + 81.0 * x3 + 81.0 * x2 - 27.0 * x - 1.35;
		case 2:	return 955.625 * x4 - 101.25 * x3 - 263.25 * x2 + 43.875 * x + 13.5;
		case 3:	return -607.5 * x4 + 378.0 * x2 - 24.5;
		case 4:	return 955.625 * x4 + 101.25 * x3 - 263.25 * x2 - 43.875 * x + 13.5;
		case 5:	return -182.25 * x4 - 81.0 * x3 + 81.0 * x2 + 27.0 * x - 1.35;
		case 6:	return 30.375 * x4 + 20.25 * x3 - 3.375 * x - 6.75 * x2 + 0.1;
		default:
			throw std::out_of_range("Invalid shape function index for Q6.");
		}

	}



	//-------------------------------------------
	//------- FE shape third derivative ---------
	//-------------------------------------------

	// ---------------- degree 1 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_third_derivative_degree_1(const unsigned int index,
		const double qPointValue) const
	{
		if (index >1)
		{
			throw std::out_of_range("Invalid shape function index for Q1.");
		}
		return 0.0;

	}

	// ---------------- degree 2 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_third_derivative_degree_2(const unsigned int index,
		const double qPointValue) const
	{
		if (index >2)
		{
			throw std::out_of_range("Invalid shape function index for Q2.");
		}
		return 0.0;

	}

	// ---------------- degree 3 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_third_derivative_degree_3(const unsigned int index,
		const double qPointValue) const
	{
		switch (index)
		{
		case 0: return 3.375;
		case 1: return 10.125;
		case 2: return 10.125;
		case 3: return 3.375;
		default:
			throw std::out_of_range("Invalid shape function index for Q3.");
		}

	}

	// ---------------- degree 4 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_third_derivative_degree_4(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		switch (index)
		{
		case 0: return 16.0 * x - 4.0;
		case 1: return - 64.0 * x + 8.0;
		case 2:	return 96.0 * x;
		case 3:	return - 64.0 * x - 8.0;
		case 4:	return 16.0 * x + 4.0;
		default:
			throw std::out_of_range("Invalid shape function index for Q4.");
		}

	}

	// ---------------- degree 5 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_third_derivative_degree_5(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		const double x2 = x * x;
		switch (index)
		{
		case 0:	return 625.0 / 32.0 * x - 3125.0 / 64.0 * x2 + 125.0 / 64.0;
		case 1:	return 15625.0 / 64.0 * x2 - 1875.0 / 32.0 * x - 1625.0 / 64.0;
		case 2:	return -15625.0 / 32.0 * x2 + 625.0 / 16.0 * x + 2125.0 / 32.0;
		case 3:	return 15625.0 / 32.0 * x2 + 625.0 / 16.0 * x - 2125.0 / 32.0;
		case 4:	return -15625.0 / 64.0 * x2 - 1875.0 / 32.0 * x + 1625.0 / 64.0;
		case 5:	return 625.0 / 32.0 * x + 3125.0 / 64.0 * x2 - 125.0 / 64.0;
		default:
			throw std::out_of_range("Invalid shape function index for Q5.");
		}

	}

	// ---------------- degree 6 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_third_derivative_degree_6(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		const double x2 = x * x;
		const double x3 = x2 * x;
		switch (index)
		{
		case 0:	return 121.5 * x3 - 60.75 * x2 + 3.375 - 13.5 * x;
		case 1:	return -729.0 * x3 + 243.0 * x2 + 162.0 * x - 27.0;
		case 2:	return 3822.5 * x3 - 303.75 * x2 - 526.5 * x + 43.875;
		case 3:	return -2430.0 * x3 + 756.0 * x;
		case 4:	return 3822.5 * x3 + 303.75 * x2 - 526.5 * x - 43.875;
		case 5:	return -729.0 * x3 - 243.0 * x2 + 162.0 * x + 27.0;
		case 6:	return 121.5 * x3 + 60.75 * x2 - 3.375 - 13.5 * x;
		default:
			throw std::out_of_range("Invalid shape function index for Q6.");
		}

	}

	//-------------------------------------------
	//------- FE shape fourth derivative --------
	//-------------------------------------------

	// ---------------- degree 1 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_fourth_derivative_degree_1(const unsigned int index,
		const double qPointValue) const
	{
		if (index >1)
		{
			throw std::out_of_range("Invalid shape function index for Q1.");
		}
		return 0.0;

	}

	// ---------------- degree 2 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_fourth_derivative_degree_2(const unsigned int index,
		const double qPointValue) const
	{
		if (index >2)
		{
			throw std::out_of_range("Invalid shape function index for Q2.");
		}
		return 0.0;

	}

	// ---------------- degree 3 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_fourth_derivative_degree_3(const unsigned int index,
		const double qPointValue) const
	{
		if (index >3)
		{
			throw std::out_of_range("Invalid shape function index for Q3.");
		}
		return 0.0;

	}

	// ---------------- degree 4 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_fourth_derivative_degree_4(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		switch (index)
		{
		case 0: return 16.0;
		case 1: return -64.0;
		case 2:	return 96.0;
		case 3:	return -64.0;
		case 4:	return 16.0;
		default:
			throw std::out_of_range("Invalid shape function index for Q4.");
		}

	}
	// ---------------- degree 5 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_fourth_derivative_degree_5(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		switch (index)
		{
		case 0:	return 625.0 / 32.0 - 3125.0 / 32.0 * x;
		case 1:	return 15625.0 / 32.0 * x - 1875.0 / 32.0;
		case 2:	return -15625.0 / 16.0 * x + 625.0 / 16.0;
		case 3:	return 15625.0 / 16.0 * x + 625.0 / 16.0;
		case 4:	return -15625.0 / 32.0 * x - 1875.0 / 32.0;
		case 5:	return 625.0 / 32.0 + 3125.0 / 32.0 * x;
		default:
			throw std::out_of_range("Invalid shape function index for Q5.");
		}

	}
	// ---------------- degree 6 ----------------
	template <int dim, int spacedim>
	double FE_QP<dim, spacedim>::get_QP_fourth_derivative_degree_6(const unsigned int index,
		const double qPointValue) const
	{
		const double x = qPointValue;
		const double x2 = x * x;
		switch (index)
		{
		case 0:	return 364.5 * x2 - 121.5 * x - 13.5;
		case 1:	return -2187.0 * x2 + 486.0 * x + 162.0;
		case 2:	return 11467.5 * x2 - 607.5 * x - 526.5;
		case 3:	return -7290.0 * x2 + 756.0;
		case 4:	return 11467.5 * x2 + 607.5 * x - 526.5;
		case 5:	return -2187.0 * x2 - 486.0 * x + 162.0;
		case 6:	return 364.5 * x2 + 121.5 * x - 13.5;
		default:
			throw std::out_of_range("Invalid shape function index for Q6.");
		}

	}

}

template class ShapeFunctions::FE_QP<1>;
template class ShapeFunctions::FE_QP<2>;
template class ShapeFunctions::FE_QP<3>;
template class ShapeFunctions::FE_QP<1, 2>;
template class ShapeFunctions::FE_QP<1, 3>;
template class ShapeFunctions::FE_QP<2, 3>;