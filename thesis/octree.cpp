#include "stdafx.h"

#include <vector>
#include <stdexcept>

#include "octree.h"
#include "constants.h"


Octree::TreeNode::TreeNode() :
	is_leaf(true),
	xyz(NULL), xyZ(NULL),
	xYz(NULL), xYZ(NULL),
	Xyz(NULL), XyZ(NULL),
	XYz(NULL), XYZ(NULL)
{
}

Octree::TreeNode::TreeNode(const Vec3f &set_center, float set_dimension_half) : TreeNode() {
	center = set_center;
	dimension_half = set_dimension_half;
}

Octree::TreeNode* Octree::TreeNode::get_child(const Vec3f &position) const {
	// Get the child of the node that is in the direction of position relative
	// to the center of the TreeNode. Does not check position to actually be
	// inside the bounds of the child.

	if (is_leaf) {
		throw std::logic_error("Leaf nodes do not contain children");
	}

	if (position.x < center.x) {
		if (position.y < center.y) {
			if (position.z < center.z) {
				return xyz;
			}
			else {
				return xyZ;
			}
		}
		else {
			if (position.z < center.z) {
				return xYz;
			}
			else {
				return xYZ;
			}
		}
	}
	else {
		if (position.y < center.y) {
			if (position.z < center.z) {
				return Xyz;
			}
			else {
				return XyZ;
			}
		}
		else {
			if (position.z < center.z) {
				return XYz;
			}
			else {
				return XYZ;
			}
		}
	}
}

void Octree::TreeNode::divide() {
	// Create children for the node that divide the space into 8 partitions

	if (!is_leaf) {
		throw std::logic_error("Cannot divide non-leaf node");
	}

	is_leaf = false;
	float child_dimension_half = dimension_half/2.0f;
	Vec3f
		dx(child_dimension_half, 0.0f, 0.0f),
		dy(0.0f, child_dimension_half, 0.0f),
		dz(0.0f, 0.0f, child_dimension_half);

	xyz = new TreeNode(center-dx-dy-dz, child_dimension_half);
	xyZ = new TreeNode(center-dx-dy+dz, child_dimension_half);
	xYz = new TreeNode(center-dx+dy-dz, child_dimension_half);
	xYZ = new TreeNode(center-dx+dy+dz, child_dimension_half);
	Xyz = new TreeNode(center+dx-dy-dz, child_dimension_half);
	XyZ = new TreeNode(center+dx-dy+dz, child_dimension_half);
	XYz = new TreeNode(center+dx+dy-dz, child_dimension_half);
	XYZ = new TreeNode(center+dx+dy+dz, child_dimension_half);
}

bool Octree::TreeNode::intersects(Vec3f poi, float distance) const {
	// Checks whether the tree node volume intersects with a cube with center poi and sides 2*distance.

	if (
		(poi.x + distance > this->center.x - this->dimension_half ||
		 poi.x - distance < this->center.x + this->dimension_half) &&
		(poi.y + distance > this->center.y - this->dimension_half ||
		 poi.y - distance < this->center.y + this->dimension_half) &&
		(poi.z + distance > this->center.z - this->dimension_half ||
		 poi.z - distance < this->center.z + this->dimension_half))
		return true;
	else
		return false;
}


void Octree::add_particle(TreeNode *r, const ParticleSystem::Particle *p) {
	// Add particle to tree with root r.

	const size_t leaf_capacity = 10;		// Max number of particles in leaf node

	if (r->is_leaf) {
		// If given particle is not a leaf node, add the particle to the sub-tree whose
		// root is the appropriate child of r.
		add_particle(r->get_child(p->position), p);
	}
	else {
		// If r is a leaf node, check whether it is full.
		if (r->particles.size() <= leaf_capacity) {
			// If not, just add the particle to the node.
			r->particles.push_back(p);
		}
		else {
			// Otherwise divede the node into 8 children and move all particles to
			// the appropriate child node.
			r->divide();
			for (size_t i = 0; i < leaf_capacity; i++) {
				add_particle(r->get_child(r->particles[i]->position), r->particles[i]);
			}
			add_particle(r->get_child(p->position), p);
		}
	}
}

void Octree::append_neighbours(const ParticleSystem::Particle &p, float distance, TreeNode* const r, std::vector<ParticleSystem::Particle*> &neighbours) const {
	// Append all particles in sub-tree r that are possibly within distance of p to neighbours.

	// First check if the space around p intersects with the space represented by the tree root.
	// If not, the entire sub-tree contains particles that are too far away.
	if (r->intersects(p.position, distance)) {
		if (r->is_leaf) {
			// If r is a leaf node, consider all particles possible neighbours of p
			neighbours.insert(neighbours.end(), r->particles.begin(), r->particles.end());
		}
		else {
			// If non-leaf, then recursively append_neighbours of each sub-tree
			append_neighbours(p, distance, r->xyz, neighbours);
			append_neighbours(p, distance, r->xyZ, neighbours);
			append_neighbours(p, distance, r->xYz, neighbours);
			append_neighbours(p, distance, r->xYZ, neighbours);
			append_neighbours(p, distance, r->Xyz, neighbours);
			append_neighbours(p, distance, r->XyZ, neighbours);
			append_neighbours(p, distance, r->XYz, neighbours);
			append_neighbours(p, distance, r->XYZ, neighbours);
		}
	}
}

void Octree::destroy_tree(TreeNode* &r) {
	// Recursively destroys all sub-trees of given tree.
	// Needs the reference to the pointer of r so as to change the pointer to NULL afterwards!

	if (!r) {
		// Technically safe, but maybe add something here to debug code that tries to destroy non-existent tree.
	}

	if (!r->is_leaf) {
		destroy_tree(r->xyz);
		destroy_tree(r->xyZ);
		destroy_tree(r->xYz);
		destroy_tree(r->xYZ);
		destroy_tree(r->Xyz);
		destroy_tree(r->XyZ);
		destroy_tree(r->XYz);
		destroy_tree(r->XYZ);
	}
	delete r;
	r = NULL;
}

Octree::Octree() : root(NULL) {
}

Octree::Octree(const ParticleSystem &ps) {
	construct_tree(ps);
}

Octree::~Octree() {
	destroy_tree(root);
}

void Octree::construct_tree(const ParticleSystem &ps) {
	// Construct tree for given ps. Destroy any previously stored tree.

	destroy_tree(root);
	for (size_t i; i < num_of_particles; i++) {
		add_particle(root, &(ps.particles[i]));
	}
}

std::vector<ParticleSystem::Particle*> Octree::find_neighbours(const ParticleSystem::Particle &p, float distance) const {
	// Returns a vector of all particles that are possibly within distance of p.

	std::vector<ParticleSystem::Particle*> neighbours;
	append_neighbours(p, distance, root, neighbours);
	return neighbours;
}