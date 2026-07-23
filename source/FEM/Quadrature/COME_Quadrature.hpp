#pragma once
#include <iostream>
#include <vector>
#include <cmath>

namespace Quadrature
{
	template<int dim>
	class Quadrature
	{
	public:
		Quadrature(int order) 
			: quadOrder_(order) {};
		int getQuadratureOrder() const;
		virtual std::vector<double>  getQuadratureWeights() const = 0;
		virtual std::vector<double> getQuadraturePoints() const = 0;

	protected:
		int quadOrder_;
	};

	template<int dim>
	int Quadrature<dim>::getQuadratureOrder() const
	{
		return quadOrder_;
	}



}