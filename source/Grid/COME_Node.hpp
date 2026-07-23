#pragma once
#include "COME_AbsTopologicalComponent.hpp"
#include <iostream>
#include <array>

namespace Mesh
{
	template<int dim, int spacedim = dim>
	class Node : public AbsTopologicalComponent<dim,spacedim>
	{
	public:
		Node() = default;
		Node(int id, std::array<double, spacedim> position)
		: position_(position)
		{}
		const std::array<double, spacedim> getCoordinates() const ;
		std::vector<Node<dim, spacedim>*> getNodes() const;
	private:
		std::array<double,spacedim> position_;
		Node* subNode_ = nullptr;

	};

	
}