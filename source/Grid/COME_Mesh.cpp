#include "COME_Mesh.hpp"

namespace Mesh
{

	template<int dim, int spacedim>
	std::vector<std::string> Mesh<dim,spacedim>::split_csv(const std::string& line) const{
		std::vector<std::string> tokens;
		std::stringstream ss(line);		//extract values from the line
		std::string token;				

		while (std::getline(ss, token, ',')) {
			token.erase(0, token.find_first_not_of(" \t"));			//remove all leading whitespaces
			token.erase(token.find_last_not_of(" \t") + 1);			//remove all trailing whitespaces
			tokens.push_back(token);								
		}
		return tokens;
	}

	template <int dim, int spacedim>
	void Mesh<dim, spacedim>::addEdgeIfMissing(int nodeA, int nodeB)
	{
		auto [minNode, maxNode] = std::minmax(nodeA, nodeB);
		std::pair key{ minNode, maxNode };

		auto [it, inserted] = edgeMap_.try_emplace(key, nullptr);
		if (inserted)
		{
			listOfEdges_.emplace_back(std::make_unique<Edge<dim, spacedim>>());
			Edge<dim, spacedim>* ptr = listOfEdges_.back().get();
			listOfEdges_.back()->linkNodes(std::array{ minNode,maxNode }, nodeIdMap_);
			edgeMap_[key] = ptr;
		}
	}
	template <int dim, int spacedim>
	void Mesh<dim, spacedim>::addFaceIfMissing(int nodeA, int nodeB, int nodeC, int nodeD)
	{
		// 1. Keep the original order for connectivity
		std::array<int, 4> connectivity{ nodeA, nodeB, nodeC, nodeD };

		// 2. Create a sorted version for the Map Key
		if (nodeA > nodeB) std::swap(nodeA, nodeB);
		if (nodeC > nodeD) std::swap(nodeC, nodeD);
		if (nodeA > nodeC) std::swap(nodeA, nodeC);
		if (nodeB > nodeD) std::swap(nodeB, nodeD);
		if (nodeB > nodeC) std::swap(nodeB, nodeC);

		std::array<int, 4> key{ nodeA, nodeB, nodeC, nodeD };

		auto [it, inserted] = faceMap_.try_emplace(key, nullptr);

		if (inserted)
		{
			auto newFace = std::make_unique<Face<dim, spacedim>>();

			// Link nodes using original order so the face knows its shape
			newFace->linkNodes(connectivity, nodeIdMap_);

			// Link edges using original order so it finds (A,B), (B,C), etc.
			newFace->linkEdges(connectivity, edgeMap_);

			it->second = newFace.get(); // Don't forget to update the map pointer!
			listOfFaces_.push_back(std::move(newFace));
		}

	}

	template<int dim, int spacedim>
	void Mesh<dim, spacedim>::addEdge(const std::array<int, 1 << dim>& nodes)
	{
		if constexpr (dim == 2) {
			for (int i = 0; i < 4; ++i) {
				int n1 = nodes[i];
				int n2 = nodes[(i + 1) % 4];

				// CRITICAL: Must match the lookup sorting!
				auto [minN, maxN] = std::minmax(n1, n2);

				if (edgeMap_.find({ minN, maxN }) == edgeMap_.end()) {
					auto edge = std::make_unique<Edge<dim, spacedim>>();
					// setup your edge ...
					edgeMap_[{minN, maxN}] = edge.get();
					listOfEdges_.push_back(std::move(edge));
				}
			}
		}

		else if constexpr (dim == 3)
		{
			// A 3D Hexahedron (Abaqus C3D8) has 12 edges.
			// Nodes 0,1,2,3 are the bottom face. Nodes 4,5,6,7 are the top face.
			std::array<std::pair<int, int>, 12> hex_edge_indices = { {
					// Bottom face edges
					{0, 1}, {1, 2}, {2, 3}, {3, 0},
					// Top face edges
					{4, 5}, {5, 6}, {6, 7}, {7, 4},
					// Vertical edges connecting bottom and top
					{0, 4}, {1, 5}, {2, 6}, {3, 7}
				} };

			for (const auto& indices : hex_edge_indices) {
				int n1 = nodes[indices.first];
				int n2 = nodes[indices.second];

				// CRITICAL: Must match the lookup sorting!
				auto [minN, maxN] = std::minmax(n1, n2);

				if (edgeMap_.find({ minN, maxN }) == edgeMap_.end()) {
					auto edge = std::make_unique<Edge<dim, spacedim>>();
					// setup your edge ...
					edgeMap_[{minN, maxN}] = edge.get();
					listOfEdges_.push_back(std::move(edge));
				}
			}
		}
	}

	template<int dim, int spacedim>
	void Mesh<dim, spacedim>::addNode(int index, const std::array<double, spacedim>& positions)
	{
		listOfNodes_.emplace_back(std::make_unique<Node<dim, spacedim>>(index, positions));
		nodeIdMap_[index] = listOfNodes_.back().get();
	}


	template <int dim, int spacedim>
	void Mesh<dim, spacedim>::addFace(const std::array<int, 1 << dim>& nodes)
	{
		// 3D Hexahedron (Abaqus C3D8)
		if constexpr (dim == 3)
		{
			// Define the 6 faces using STRICT perimeter loop winding.
			// This ensures Face::linkEdges never creates a diagonal.
			std::array<std::array<int, 4>, 6> hex_faces = { {
				{nodes[0], nodes[1], nodes[2], nodes[3]}, // Bottom face
				{nodes[4], nodes[5], nodes[6], nodes[7]}, // Top face
				{nodes[0], nodes[1], nodes[5], nodes[4]}, // Front face
				{nodes[1], nodes[2], nodes[6], nodes[5]}, // Right face
				{nodes[2], nodes[3], nodes[7], nodes[6]}, // Back face
				{nodes[3], nodes[0], nodes[4], nodes[7]}  // Left face
			} };

			for (const auto& face_nodes : hex_faces)
			{
				// 1. Create a sorted copy of the nodes to use strictly as the Map Key
				std::array<int, 4> key = face_nodes;
				std::sort(key.begin(), key.end());

				// 2. Check if the face already exists
				if (faceMap_.find(key) == faceMap_.end())
				{
					auto face = std::make_unique<Face<dim, spacedim>>();

					// CRITICAL: Pass the original PERIMETER-wound array to the link functions, 
					// NOT the sorted key!
					face->linkNodes(face_nodes, nodeIdMap_);
					face->linkEdges(face_nodes, edgeMap_);

					// 3. Store in maps using the universally sorted key
					faceMap_[key] = face.get();
					listOfFaces_.push_back(std::move(face));
				}
			}
		}
	}

	template<int dim, int spacedim>
	void Mesh<dim, spacedim>::addElement(int elementnumber, const std::array<int, 1 << dim>& nodes)
	{
		auto element = std::make_unique<Element<dim, spacedim>>();
		Element<dim, spacedim>* elPtr = element.get();

		if constexpr (dim == 1) {
			auto edge = std::make_unique<Edge<dim, spacedim>>();
			elPtr->setTopologicalSupport(edge.get());
			listOfEdges_.push_back(std::move(edge));
		}
		else if constexpr (dim == 2) {
			auto face = std::make_unique<Face<dim, spacedim>>();
			elPtr->setTopologicalSupport(face.get());
			listOfFaces_.push_back(std::move(face));
		}
		else if constexpr (dim == 3) {
			auto volume = std::make_unique<Volume<dim, spacedim>>();
			elPtr->setTopologicalSupport(volume.get());
			listOfVolumes_.push_back(std::move(volume));
		}

		// Build standard perimeter-wound topology tracking
		if constexpr (dim == 2) {
			addEdge(nodes); // Native Abaqus loop
		}
		else if constexpr (dim == 3) {
			addEdge(nodes);
			addFace(nodes);
		}

		// Link steps get the unaltered native array
		if constexpr (dim == 1)
		{
			listOfEdges_.back()->linkNodes(nodes, nodeIdMap_);
		}
		else if constexpr (dim == 2)
		{
			listOfFaces_.back()->linkNodes(nodes, nodeIdMap_);
			listOfFaces_.back()->linkEdges(nodes, edgeMap_); // Success! No more diagonal lookups.
		}
		else if constexpr (dim == 3)
		{
			listOfVolumes_.back()->linkNodes(nodes, nodeIdMap_);
			listOfVolumes_.back()->linkEdges(nodes, edgeMap_);
			listOfVolumes_.back()->linkFaces(nodes, faceMap_);
		}

		listOfElements_.push_back(std::move(element));
	}


	/*
	Generally for all read_mesh functions, we make the following assumption:
	The DOFs are filled from x to z. So e.g. when we have a 2d mesh, we have x and y and not x and z.
	This makes my life significantly easier cause i don't need to make sure i account for all the different directions my mesh could possibly be oriented towards.
	It also makes sense cause one can always rotate the mesh properly
	We also for right now assume that we don't need to do line skips and stuff like that, so i guess we start with very basic meshes for now.
	*/
	template<int dim, int spacedim>
	void Mesh<dim, spacedim>::read_abaqus(std::istream& in)
	{
		
		if (!in) {
			std::cerr << "File could not be opened. Current directory: "
				<< std::filesystem::current_path() << "\n";
			throw std::runtime_error("File could not be opened");
		}


		std::string line;

		AbaqusMeshSection section = AbaqusMeshSection::NONE;
		//Loop through all lines
		while (std::getline(in, line)) {
			
			if (line.empty()) continue; //ignore empty lines

			if (line[0] == '*') {						//now if we actually hit a new type, we need to change the respective type we read in
				if (line.find("*NODE") == 0) {
					section = AbaqusMeshSection::NODE;
					
				}
				else if (line.find("*ELEMENT") == 0) {
					section = AbaqusMeshSection::ELEMENT;
				}
				else {
					section = AbaqusMeshSection::NONE;
				}
				continue;
			}

			auto fields = split_csv(line);			//split the line

			//If we are reading in a Node
			if (section == AbaqusMeshSection::NODE)
			{
				std::array<double, spacedim> position;
				for (int i = 0; i < spacedim; i++)
				{
					position[i] = std::stod(fields[i + 1]);		//convert all positions from strings to integers
				}
				this->addNode(std::stoi(fields[0]),position);
			}
			if (section == AbaqusMeshSection::ELEMENT)
			{
				std::array<int, 1 << dim> raw_nodes;
				for (int i = 0; i < (1 << dim); i++)
				{
					raw_nodes[i] = std::stoi(fields[i + 1]);
				}

				// Convert Abaqus winding to Lexicographic tensor-product ordering
				std::array<int, 1 << dim> lex_nodes;

				if constexpr (dim == 1) {
					lex_nodes = raw_nodes;
				}
				else if constexpr (dim == 2) {
					// Abaqus: 0(BottomLeft), 1(BottomRight), 2(TopRight), 3(TopLeft)
					// Lexicographic: 0(BL), 1(BR), 2(TL), 3(TR)
					lex_nodes = { raw_nodes[0], raw_nodes[1], raw_nodes[3], raw_nodes[2] };
				}
				else if constexpr (dim == 3) {
					// Abaqus bottom face CCW (0-3), top face CCW (4-7)
					// Lexicographic (x fastest, then y, then z)
					lex_nodes = {
						raw_nodes[0], raw_nodes[1], raw_nodes[3], raw_nodes[2], // z = 0
						raw_nodes[4], raw_nodes[5], raw_nodes[7], raw_nodes[6]  // z = 1
					};
				}

				// Now all downstream operations use tensor-product ordering
				this->addElement(std::stoi(fields[0]), lex_nodes);
			}

		}
	}

	template<int dim, int spacedim>
	void Mesh<dim, spacedim>::read_msh(std::istream& in)
	{

	}

	template<int dim, int spacedim>
	const std::unordered_map<int, Node<dim, spacedim>*>& Mesh<dim, spacedim>::getNodeMap() const
	{
		return nodeIdMap_;
	}

	template <int dim, int spacedim>
	const std::vector<std::unique_ptr<Node<dim, spacedim>>>& Mesh<dim, spacedim>::getNodes() const
	{
		return listOfNodes_;
	}

	template <int dim, int spacedim>
	const std::vector<std::unique_ptr<Edge<dim, spacedim>>>& Mesh<dim, spacedim>::getEdges() const
	{
		return listOfEdges_;
	}

	template <int dim, int spacedim>
	const std::vector<std::unique_ptr<Face<dim, spacedim>>>& Mesh<dim, spacedim>::getFaces() const
	{
		return listOfFaces_;
	}

	template <int dim, int spacedim>
	const std::vector<std::unique_ptr<Volume<dim, spacedim>>>& Mesh<dim, spacedim>::getVolumes() const
	{
		return listOfVolumes_;
	}

	template <int dim, int spacedim>
	const std::vector<std::unique_ptr<Element<dim, spacedim>>>& Mesh<dim,spacedim>::getElements() const
	{
		return listOfElements_;
	}

}

template class Mesh::Mesh<1>;
template class Mesh::Mesh<2>;
template class Mesh::Mesh<3>;
template class Mesh::Mesh<1, 2>;
template class Mesh::Mesh<1, 3>;
template class Mesh::Mesh<2, 3>;