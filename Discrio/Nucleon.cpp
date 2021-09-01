#include "Nucleon.h"



Nucleon::Nucleon(int species, int x_coord, int y_coord)
{
	type = species;
	x = x_coord;
	y = y_coord;
	mass = 1;

}


Nucleon::~Nucleon()
{
}
