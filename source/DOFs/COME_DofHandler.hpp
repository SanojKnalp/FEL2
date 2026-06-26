#pragma once
#include "Grid/COME_Mesh.hpp"
#include "FEM/ShapeFunctions/COME_ShapeFunctions.hpp"
#include "COME_DofUtils.hpp"

#include <vector>
#include <memory>

namespace DOFs
{
	enum class EntityType { Node, Edge, Face, Volume };

	struct LocalDoFMetadata {
		EntityType type;
		unsigned int entity_index;  // e.g., which edge (0, 1, 2, 3)
		unsigned int internal_index; // index inside that specific edge/face
	};

	template <int dim, int spacedim = dim>
	class DoFHandler
	{
	public:
		DoFHandler(Mesh::Mesh<dim, spacedim>& mesh) : mesh_(mesh) {}

		void distribute_dofs(const ShapeFunctions::ShapeFunctions<dim, spacedim>& fe);
		void get_dof_indices(const auto& cell, std::vector<unsigned int>& local_to_global_map) const;
		unsigned int number_of_dofs() const;
	private:
		Mesh::Mesh<dim, spacedim>& mesh_;

		std::vector<std::unique_ptr<DOF<dim, spacedim>>> DoFList_;

		/* here we need to implement a structure which gets the correct local numbering for each cell
		Essentially what we need to extract find a way to get the local and global numbering working for each cell.
		So ideally we need something like a "Reference cell" structure that 
		*/ 
		std::vector<LocalDoFMetadata> local_dof_pattern_;



	};
	
}