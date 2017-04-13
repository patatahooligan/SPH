========================================================================
							PROJECT STRUCTURE
========================================================================

thesis.cpp
	Holds the global ParticleSystem object and the main() function of the project

physics.cpp & physics.h
	Definition of class ParticleSystem which serves as a context for the whole
	simulation. It holds all the Particles (defined in particles.h) and offers
	all necessary functions for the physics simulation.

particle.h
	Definition of Particle class.

render.cpp & render.h
	Functions for rendering the current state of ParticleSystem using freeglut.

vec3f.cpp & vec3f.h
	Definition of a 3D float vector that is used for the physics emulation.

octree.cpp and octree.h
	An octree implementation to be used for neighbour searching to reduce
	the complexity of some methods in physics.h.

video.cpp and video.h
	Definition of Video class that can be used to save the output from render.cpp
	to a video file.

constants.h
	Used to hold compile time constants that are needed across multiple files.

stdafx.cpp & stafx.h
	These files are used to build a precompiled header (PCH) file
	named thesis.pch and a precompiled types file named StdAfx.obj.