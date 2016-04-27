#include "Particle.h"

Particle::Particle(double _mass, double _sigma, double _epsilon, double _lambda) :
x(0.0, DIM), U2(0.0, DIM), U3(0.0, DIM)
{
	mass = _mass;
	sigma = _sigma;
	epsilon = _epsilon;
	lambda = _lambda;
}

bool Particle::operator!= (const Particle &_a){
	if (id != _a.id) return true;
	else return false;
}

