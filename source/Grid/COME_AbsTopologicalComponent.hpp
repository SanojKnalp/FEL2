#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <stdexcept>

//forward declaration:
namespace FEM {
	template <int dim, int spacedim>
	class FEValues;
}


namespace Mesh
{
	template <int dim, int spacedim>
	class Element;

	template <int dim, int spacedim>
	class Node;

	template <int dim, int spacedim = dim>
	class AbsTopologicalComponent
	{
	public:
		virtual ~AbsTopologicalComponent() = default;

		void activate() { isActive_ = true; };
		void addAdjacentElement(Element<dim, spacedim>* el);
		const std::vector<Element<dim, spacedim>*>& getAdjacentElements() const;
		void add_dof(unsigned int id);
		const std::vector<unsigned int>& get_dofs() const;
		void clear_dofs();

		virtual std::vector<Node<dim, spacedim>*> getNodes() const = 0;

	protected:
		unsigned int polynomialDegree_ = 1;
		bool isActive_ = false;
		std::vector<Element<dim,spacedim>*> listOfAdjacentElements_;
		unsigned int refinementLevel_ = 0;
		std::vector<unsigned int> dof_ids_;

		friend class FEM::FEValues<dim,spacedim>;

	};

	template<int dim, int spacedim>
	void AbsTopologicalComponent<dim, spacedim>::addAdjacentElement(Element<dim, spacedim>* el)
	{
		listOfAdjacentElements_.push_back(el);
	}
	template<int dim, int spacedim>
	const std::vector<Element<dim, spacedim>*>& AbsTopologicalComponent<dim, spacedim>::getAdjacentElements() const
	{
		return listOfAdjacentElements_;
	}

	template<int dim, int spacedim>
	void AbsTopologicalComponent<dim, spacedim>::add_dof(unsigned int id)
	{
		dof_ids_.push_back(id);
	}

	template<int dim, int spacedim>
	const std::vector<unsigned int>& AbsTopologicalComponent<dim, spacedim>::get_dofs() const
	{
		return dof_ids_;
	}

	template<int dim, int spacedim>
	void AbsTopologicalComponent<dim, spacedim>::clear_dofs()
	{
		dof_ids_.clear();
	}

}