#pragma once
#include "Nucleon.h"
#include <vector>

class Particle_Manager
{
public:
	Particle_Manager();
	~Particle_Manager();
public:
	std::vector<Nucleon> Particles;
	int width;
	int height;

private:
	std::vector<std::vector<std::vector<Nucleon>>> Space

};

