
#ifndef N_BODY_POT
#define N_BODY_POT

#include <math.h>
#include <list>
#include <vector>
#include "Particle.h"
#include "Utils.h"

class N_Body_pot{

	std::list<pot_sample> Values;

public:
	N_Body_pot();

	void S_W_Pot(std::list<Particle> p_list, double _A, double _B, double _Gamma, double _a, double _p, double _q);

	void Sort_Res();

	void print_res();

	void print_VTK();

	void print_VTK_Struct_Points(double x, double y, double z, double offx, double offy, double offz);

	void print_VTK_Struc_Grid(double x, double y, double z, double offx, double offy, double offz);

	void print_VTK_Rect_Grid(double x, double y, double z, double offx, double offy, double offz, double step);

};




#endif