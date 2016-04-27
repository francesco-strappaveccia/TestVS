#include "nbody_pot.h"
#include "Particle.h"
#include <list>
#include <random>

void Test_Pot();

void Cube_test();

double Distance_m(Particle& i, Particle& j);



int main(){

	//Test_Pot();
	Cube_test();

	return 0;
}

void Test_Pot(){

	double max_dist = 2;

	// Random Generator
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution1(9.5, 10.5);
	std::uniform_real_distribution<double> distribution2(9.5, 10.5);
	std::uniform_real_distribution<double> distribution3(9.5, 10.5);

	std::list<Particle> p_list;
	Particle a(28.085, 0.20951, 2.1678, 21.0); 
	a.id = 0;
	a.x[0] = 10.0; a.x[1] = 10.0; a.x[2] = 10.0;
	p_list.push_back(a);

	Particle b(28.085, 0.20951, 2.1678, 21.0);
	b.id = 1;
	b.x[0] = 15.0; b.x[1] = 15.0; b.x[2] = 15.0;
	p_list.push_back(b);


	Particle c(28.085, 0.20951, 2.1678, 21.0);
	c.id = 2;

	N_Body_pot Pot;

	// Moving particle
	int samples = 5000;
	int s = 0;
	while (s < samples){

		c.x[0] = distribution1(generator); c.x[1] = distribution2(generator); c.x[2] = distribution3(generator);
		p_list.push_back(c);

		Pot.S_W_Pot(p_list, 7.049556277, 0.6022245584, 1.2, 1.8, 4.0, 0.0);

		p_list.pop_back();

		s++;

	}

	Pot.Sort_Res();

	Pot.print_res();

	Pot.print_VTK();

}

void Cube_test(){

	// Set Cube Dimensions and step

	std::list<Particle> p_list;
	Particle a(28.085, 0.20951, 2.1678, 21.0);
	a.id = 0;
	a.x[0] = 20.0; a.x[1] = 20.0; a.x[2] = 20.0;
	p_list.push_back(a);
	
	Particle b(28.085, 0.20951, 2.1678, 21.0);
	b.id = 1;
	b.x[0] = 20.0; b.x[1] = 20.36; b.x[2] = 20.0; //19.880000 19.980000 20.140000
	p_list.push_back(b);

	Particle b1(28.085, 0.20951, 2.1678, 21.0);
	b1.id = 2;
	b1.x[0] = 20.0; b1.x[1] = 19.70; b1.x[2] = 20.0; //19.880000 19.980000 20.140000
	p_list.push_back(b1);

	//double dist = Distance_m(a, b);

	// Moving Particle
	Particle c(28.085, 0.20951, 2.1678, 21.0);
	c.id = 3;

	// Set moves in 3d space
	double off_x, off_y, off_z; // test volume offsets from central particle
	off_x = off_y = off_z = 0.4;
	double step = 0.015; // grid spacing

	// Center coordinates
	double c_x, c_y, c_z;
	c_x = a.x[0]; c_y = a.x[1];  c_z = a.x[2];

	N_Body_pot Potential;

	// test Potential
	for (double z = c_z - off_z; z < c_z + off_z; z += step){
		for (double y = c_y - off_y; y < c_y + off_y; y += step){
			for (double x = c_x - off_x; x < c_x + off_x; x += step){

				c.x[0] = x; c.x[1] = y; c.x[2] = z;
				p_list.push_back(c);

				Potential.S_W_Pot(p_list, 7.049556277, 0.6022245584, 1.2, 1.8, 4.0, 0.0);

				p_list.pop_back();

			}
		}
	}
	 
	// Print VTK File
	//Potential.print_VTK();
	Potential.print_VTK_Struct_Points(c_x, c_y, c_z, off_x, off_y, off_z); // Structures Points
	//Potential.print_VTK_Struc_Grid(c_x, c_y, c_z, off_x, off_y, off_z); // Structured Grid
	//Potential.print_VTK_Rect_Grid(c_x, c_y, c_z, off_x, off_y, off_z, step); // Rectilinear Grid

	Potential.Sort_Res();
	Potential.print_res();
	

}

double Distance_m(Particle& i, Particle& j){

	double dist2 = 0.0; double dist = 0.0;

	for (int d = 0; d < DIM; d++){
		dist2 += sqr(i.x[d] - j.x[d]);
	}

	dist = sqrt(dist2);

	return dist;


}