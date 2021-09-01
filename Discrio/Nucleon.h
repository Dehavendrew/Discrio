#pragma once
class Nucleon
{
public:
	Nucleon(int species, int x_coord, int y_coord);
	~Nucleon();

public:
	int type;
	int x;
	int y;
	int mass;
};

