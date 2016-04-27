#ifndef PARTICLE_H
#define PARTICLE_H

#include "Utils.h"

#include <valarray>

class Particle
{
public:
	std::valarray<double> x;
	std::valarray<double> U2;
	std::valarray<double> U3;
	double sigma;
	double epsilon;
	double lambda;
	double mass;
	int id;

	Particle();
	Particle(double _mass, double _sigma, double _epsilon, double _lambda);
	double Get_F();

	bool Particle::operator!= (const Particle &_a);


};


#endif