#include "COME_Quadrature_QGauss.hpp"


namespace Quadrature
{
	template <int dim>
	std::vector<double> QGauss<dim>::getQuadratureWeights() const
	{
		switch (this->getQuadratureOrder())
		{
		case 1:	return get_QGauss_weight_degree_1();
		case 2:	return get_QGauss_weight_degree_2();
		case 3:	return get_QGauss_weight_degree_3();
		case 4: return get_QGauss_weight_degree_4();
		case 5: return get_QGauss_weight_degree_5();
		case 6: return get_QGauss_weight_degree_6();
		case 7: return get_QGauss_weight_degree_7();
		default:
			throw std::out_of_range("You tried to access a quadrature order of degree 8 or higher. This is not implemented as of right now.");
		}
	}

	template<int dim>
	std::vector<double> QGauss<dim>::getQuadraturePoints() const
	{
		switch (this->getQuadratureOrder())
		{
		case 1:	return get_QGauss_points_degree_1();
		case 2:	return get_QGauss_points_degree_2();
		case 3:	return get_QGauss_points_degree_2();
		case 4: return get_QGauss_points_degree_2();
		case 5: return get_QGauss_points_degree_5();
		case 6: return get_QGauss_points_degree_6();
		case 7: return get_QGauss_points_degree_7();
		default:
			throw std::out_of_range("You tried to access a quadrature order of degree 8 or higher. This is not implemented as of right now.");
		}
	}

	//---------------------------------------------
	//----------- Quadrature weights --------------
	//---------------------------------------------

	//--------------- degree 1---------------------
	template <int dim>
	std::vector<double> QGauss<dim>::get_QGauss_weight_degree_1() const
	{
		return std::vector<double>{2.0};
	}

	//--------------- degree 2---------------------
	template <int dim>
	std::vector<double> QGauss<dim>::get_QGauss_weight_degree_2() const
	{
		return std::vector<double>{1.0, 1.0};
	}

	//--------------- degree 3---------------------
	template <int dim>
	std::vector<double> QGauss<dim>::get_QGauss_weight_degree_3() const
	{
		return std::vector<double>{5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
	}

	//--------------- degree 4---------------------
	template <int dim>
	std::vector<double> QGauss<dim>::get_QGauss_weight_degree_4() const
	{
		return std::vector<double>{(18.0 - std::sqrt(30.0)) / 36.0, (18.0 + std::sqrt(30.0)) / 36.0, (18.0 + std::sqrt(30.0)) / 36.0, (18.0 - std::sqrt(30.0)) / 36.0};
	}

	//--------------- degree 5---------------------
	template <int dim>
	std::vector<double> QGauss<dim>::get_QGauss_weight_degree_5() const
	{
		return std::vector<double>{(322.0 - 13.0 * std::sqrt(70.0)) / 900.0, (322.0 + 13.0 * std::sqrt(70.0)) / 900.0, 128.0 / 225.0, (322.0 + 13.0 * std::sqrt(70.0)) / 900.0, (322.0 - 13.0 * std::sqrt(70.0)) / 900.0};
	}

	//--------------- degree 6---------------------
	template <int dim>
	std::vector<double> QGauss<dim>::get_QGauss_weight_degree_6() const
	{
		return std::vector<double>{0.17132449237917034504029614217273, 0.36076157304813860756983351383772, 0.46791393457269104738987034398955, 0.46791393457269104738987034398955, 0.36076157304813860756983351383772, 0.17132449237917034504029614217273};
	}

	//--------------- degree 7---------------------
	template <int dim>
	std::vector<double> QGauss<dim>::get_QGauss_weight_degree_7() const
	{
		return std::vector<double>{0.12948496616886969327061143267908, 0.27970539148927666790146777142378, 0.38183005050511894495036977548898, 0.41795918367346938775510204081633, 0.38183005050511894495036977548898, 0.27970539148927666790146777142378, 0.12948496616886969327061143267908};
	}

	//---------------------------------------------
	//----------- Quadrature points  --------------
	//---------------------------------------------

	//--------------- degree 1---------------------
	template <int dim>
	std::vector<double> QGauss<dim>::get_QGauss_points_degree_1() const
	{
		return std::vector<double>{0.0};
	}

	//--------------- degree 2---------------------
	template <int dim>
	std::vector<double> QGauss<dim>::get_QGauss_points_degree_2() const
	{
		return std::vector<double>{-0.57735026918962576450914878050196, 0.57735026918962576450914878050196};
	}

	//--------------- degree 3---------------------
	template <int dim>
	std::vector<double> QGauss<dim>::get_QGauss_points_degree_3() const
	{
		return std::vector<double>{-0.77459666924148337703585307995648, 0.0, 0.77459666924148337703585307995648};
	}

	//--------------- degree 4---------------------
	template <int dim>
	std::vector<double> QGauss<dim>::get_QGauss_points_degree_4() const
	{
		return std::vector<double>{-0.86113631159405257522394648889281, -0.33998104358485626480266575910324, 0.33998104358485626480266575910324,0.86113631159405257522394648889281};
	}

	//--------------- degree 5---------------------
	template <int dim>
	std::vector<double> QGauss<dim>::get_QGauss_points_degree_5() const
	{
		return std::vector<double>{-0.90617984593866399279762687829939, -0.53846931010568309103631442070021, 0.0, 0.53846931010568309103631442070021,0.90617984593866399279762687829939};
	}

	//--------------- degree 6---------------------
	template <int dim>
	std::vector<double> QGauss<dim>::get_QGauss_points_degree_6() const
	{
		return std::vector<double>{-0.93246951420315202781230155449399, -0.66120938646626451366139959501991, -0.23861918608319690863050172168071, 0.23861918608319690863050172168071,0.66120938646626451366139959501991,0.93246951420315202781230155449399};
	}

	//--------------- degree 7---------------------
	template <int dim>
	std::vector<double> QGauss<dim>::get_QGauss_points_degree_7() const
	{
		return std::vector<double>{-0.94910791234275852452618968404785, -0.74153118559939443986386477328079, -0.40584515137739716690660641207696, 0.0, 0.40584515137739716690660641207696,0.74153118559939443986386477328079,0.94910791234275852452618968404785};
	}


}

template class Quadrature::QGauss<1>;
template class Quadrature::QGauss<2>;
template class Quadrature::QGauss<3>;