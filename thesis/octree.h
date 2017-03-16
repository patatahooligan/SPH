#pragma once

#include <vector>

#include "vec3f.h"

class ParticleSystem;
class Particle;

class Octree {
	private:
		struct TreeNode {
			public:
				TreeNode
					*xyz, *xyZ,
					*xYz, *xYZ,
					*Xyz, *XyZ,
					*XYz, *XYZ;

				bool is_leaf;
				Vec3f center;
				float dimension_half;

				std::vector<const Particle *> particles;

				TreeNode();
				TreeNode(const Vec3f &set_center, float set_dimension_half);

				TreeNode* get_child(const Vec3f &position) const;

				void divide();

				bool intersects(Vec3f poi, float distance) const;
		};

		TreeNode *root;

		void add_particle(TreeNode *r, const Particle *p);

		void append_neighbours(const Particle &p, float distance, TreeNode* const r, std::vector<const Particle*> &neighbours) const;

		void destroy_tree(TreeNode* &r);

	public:
		Octree();

		Octree(const ParticleSystem &ps);

		~Octree();

		void construct_tree(const ParticleSystem &ps);

		std::vector<const Particle*> find_neighbours(const Particle &p, float distance) const;
};