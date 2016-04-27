#include <iostream>
#include <string>
#include "nbody_pot.h"

double GetTheta(Particle& i, Particle& j, Particle& k); // i centre of angle

double GetLambda(Particle& i, Particle& j, Particle& k);
double GetLambda(Particle& i, Particle& j);

double GetSigma(Particle& i, Particle& j);

double GetEps(Particle& i, Particle& j, Particle& k);
double GetEps(Particle& i, Particle& j);

double Distance(Particle& i, Particle& j);

double U2(Particle& i, Particle& j, double A, double B, double p, double q, double a);

double U3(Particle& i, Particle& j, Particle& k, double Gamma, double a);

double h(double s_ab, double s_bc, double theta_b, double lambda, double gamma, double a); // b angle vertex s = r_ab/sigma_ab 

bool sort_sample(const pot_sample &a, const pot_sample &b){
	return a.d1 < b.d1;
}

N_Body_pot::N_Body_pot(){

}

void N_Body_pot::S_W_Pot(std::list<Particle> p_list, double _A, double _B, double _Gamma, double _a, double _p, double _q){

	double A = _A;
	double B = _B;
	double Gamma = _Gamma;
	double a = _a;
	double p = _p;
	double q = _q;

	// Main Loop

	// Sample the potential following the moving particle
	pot_sample sample{ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	auto first = p_list.begin();
	auto part = --p_list.end();

	sample.x1 = (*part).x[0];
	sample.x2 = (*part).x[1];
	sample.x3 = (*part).x[2];
	//distance from first particle
	sample.d1 = Distance((*part), (*first));

	double U = 0.0;

	for (auto i = p_list.begin(); i != p_list.end(); ++i){
		//auto next_j = i;
		for (auto j = i; j != p_list.end(); ++j){ // pair interaction
			if ((*i) != (*j)){
				/*std::cout << (*i).id;
				std::cout << (*j).id;
				std::cout << "\n";*/
				double U2_val = U2((*i), (*j), A, B, p, q, a);

				U += U2_val;
				for (auto k = j; k != p_list.end(); ++k){ // multi-body interaction
					if ((*k) != (*i) && (*k) != (*j)){
						/*std::cout << (*i).id;
						std::cout << (*j).id;
						std::cout << (*k).id;
						std::cout << "\n";*/
						double U3_val = U3((*i), (*j), (*k), Gamma, a);

						U += U3_val;
					}
				}
			}
		}
		
	}

	sample.value = U;
	
	Values.push_back(sample);
}

void N_Body_pot::Sort_Res(){

	Values.sort(&sort_sample);

}

void N_Body_pot::print_res(){

	double min, max;
	for (auto it = Values.begin(); it != Values.end(); it++){

	}

	FILE *f;
	f = fopen("results.dat", "w");

	// First Line
	fprintf(f, "#	X	Y1	Y2\n");

	for (auto it = Values.begin(); it != Values.end(); ++it){

		fprintf(f, "%f	%f	%f %f %f \n", (*it).d1, (*it).value, (*it).x1, (*it).x2, (*it).x3);

	}

	fclose(f);

}

void N_Body_pot::print_VTK(){

	FILE *f;
	std::string ext = ".vtk";
	std::string folder = "results\\";
	std::string step = std::to_string(0);
	std::string filename = "data" + step;
	std::string path = folder + filename + ext;

	fopen_s(&f, path.c_str(), "w");

	// File First Part
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "MD TEST\n");
	fprintf(f, "ASCII\n");
	// Type of data
	fprintf(f, "DATASET POLYDATA\n");

	int N = Values.size();

	// Data Points
	fprintf(f, "POINTS %d double\n", N);
	//for (int i = 0; i < N; i++){
	for (auto it = Values.begin(); it != Values.end(); ++it){
		fprintf(f, "%.20f ", (*it).x1);
		fprintf(f, "%.20f ", (*it).x2);
		fprintf(f, "%.20f ", (*it).x3);
		/*for (int c = 0; c < DIM; c++){
			fprintf(f, "%.20f ", (*it).x[c]);
		}*/
		//fprintf(f, "0.0");
		fprintf(f, "\n");
	}

	fprintf(f, "\n");

	// Vertices
	fprintf(f, "VERTICES %d %d\n", N, 2*N);
	for (int i = 0; i < N; i++){
		fprintf(f, "1 %d\n", i);
	}
	fprintf(f, "\n");

	// POINTS DATA
	fprintf(f, "POINT_DATA %d\n", N);
	// Scalar quantities
	fprintf(f, "SCALARS PotValue double\n");
	fprintf(f, "LOOKUP_TABLE default\n");
	//for (int i = 0; i < N; i++){
	for (auto it = Values.begin(); it != Values.end(); ++it){

		double val = (*it).value;

		if (val > 0.0){
			if (val > 1.0e09)
				val = 1.0e09;
		}
		else{
			if (val < -10000000.00)
				val = -10000000.00;
		}
		fprintf(f, "%.20f\n", val);
	}

	fprintf(f, "\n");

	// Vectors Quantities

	//// FORCE
	//fprintf(f, "VECTORS force double\n");
	////for (int i = 0; i < N; i++){
	//for (auto it = Values.begin(); it != Values.end(); ++it){
	//	/*for (int c = 0; c < DIM; c++){
	//		fprintf(f, "%.20f ", particles[i].F[c]);
	//	}*/
	//	fprintf(f, "%.20f ", (*it).value);
	//	fprintf(f, "0.0");
	//	fprintf(f, "\n");
	//}
	//fprintf(f, "\n");

	//// VELOCITY
	//fprintf(f, "VECTORS velocity double\n");
	//for (int i = 0; i < N; i++){
	//	for (int c = 0; c < DIM; c++){
	//		fprintf(f, "%.20f ", particles[i].v[c]);
	//	}
	//	fprintf(f, "0.0");
	//	fprintf(f, "\n");
	//}

	fclose(f);

	return;


}

void N_Body_pot::print_VTK_Struct_Points(double x, double y, double z, double offx, double offy, double offz){

	FILE *f;
	std::string ext = ".vtk";
	std::string folder = "results\\";
	std::string filename = "data_S_P";
	std::string path = folder + filename + ext;

	fopen_s(&f, path.c_str(), "w");

	int N = Values.size();

	printf("Print STRUCTURED POINTS\n");

	// File First Part
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Potential TEST\n");
	fprintf(f, "ASCII\n");
	// Type of data
	fprintf(f, "DATASET STRUCTURED_POINTS\n");
	// Volume dimensions
	double v_dim = cbrt(N);
	fprintf(f, "DIMENSIONS %d %d %d\n", (int)v_dim, (int)v_dim, (int)v_dim);
	// Aspect Ratio
	fprintf(f, "ASPECT_RATIO %d %d %d \n", 1, 1, 1);
	// Cube origin
	fprintf(f, "ORIGIN %d %d %d\n", 0, 0, 0);
	// Points
	fprintf(f, "POINT_DATA %d\n", N);
	// Scalars
	fprintf(f, "SCALARS volume_scalars float 1\n");
	// Look up table
	fprintf(f, "LOOKUP_TABLE default\n");

	// Print Values
	int count = 0;
	for (auto it = Values.begin(); it != Values.end(); ++it, count++){
		double val = (*it).value;

		if (val > 0.0){
			if (val > 1.0e09)
				val = 1.0e09;
		}
		else{
			if (val < -10000000.00)
				val = -10000000.00;
		}
			

		fprintf(f, "%e ", val);
		if (count%10 == 0 && count > 0)
			fprintf(f, "\n");

	}

	

}

void N_Body_pot::print_VTK_Struc_Grid(double x, double y, double z, double offx, double offy, double offz){

	FILE *f;
	std::string ext = ".vtk";
	std::string folder = "results\\";
	std::string step = std::to_string(0);
	std::string filename = "data_S_G";
	std::string path = folder + filename + ext;

	fopen_s(&f, path.c_str(), "w");

	int N = Values.size();

	// File First Part
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Potential TEST\n");
	fprintf(f, "ASCII\n");
	// Type of data
	fprintf(f, "DATASET STRUCTURED_GRID\n");
	// Volume dimensions
	double v_dim = cbrt(N);
	fprintf(f, "DIMENSIONS %d %d %d\n", (int)v_dim, (int)v_dim, (int)v_dim);
	// Points
	fprintf(f, "POINTS %d float\n", N);
	// Point Coords
	for (auto it = Values.begin(); it != Values.end(); ++it){
		fprintf(f, "%f %f %f\n", (*it).x1, (*it).x2, (*it).x3);
	}

	fprintf(f, "\n");


	fprintf(f, "POINT_DATA %d\n", N);
	// Scalars
	fprintf(f, "SCALARS Pot_Values float 1\n");
	// Look up table
	fprintf(f, "LOOKUP_TABLE default\n");

	// Print Values
	int count = 0;
	for (auto it = Values.begin(); it != Values.end(); ++it, count++){
		double val = (*it).value;

		if (val > 0.0){
			if (val > 1.0e09)
				val = 1.0e09;
		}
		else{
			if (val < -10000000.00)
				val = -10000000.00;
		}


		fprintf(f, "%f ", val);
		if (count % 10 == 0 && count > 0)
			fprintf(f, "\n");

	}
}

void  N_Body_pot::print_VTK_Rect_Grid(double x, double y, double z, double offx, double offy, double offz, double step){

	FILE *f;
	std::string ext = ".vtk";
	std::string folder = "results\\";
	std::string filename = "data_R_G";
	std::string path = folder + filename + ext;

	fopen_s(&f, path.c_str(), "w");

	int N = Values.size();

	printf("Printing Rectangular Grid\n");

	// File First Part
	fprintf(f, "# vtk DataFile Version 2.0\n");
	fprintf(f, "Potential TEST\n");
	fprintf(f, "ASCII\n");
	// Type of data
	fprintf(f, "DATASET RECTILINEAR_GRID\n");
	// Volume dimensions
	double v_dim = cbrt(N);
	fprintf(f, "DIMENSIONS %d %d %d\n", (int)v_dim, (int)v_dim, (int)v_dim);

	// Volume grid steps
	// X
	fprintf(f, "X_COORDINATES %d float\n", (int)v_dim);
	for (double v = x - offx; v < x + offx; v+=step){
		fprintf(f, "%f ", v);
	}
	fprintf(f, "\n");
	// Y
	fprintf(f, "Y_COORDINATES %d float\n", (int)v_dim);
	for (double v = x - offx; v < x + offx; v += step){
		fprintf(f, "%f ", v);
	}
	fprintf(f, "\n");
	// Z
	fprintf(f, "Z_COORDINATES %d float\n", (int)v_dim);
	for (double v = x - offx; v < x + offx; v += step){
		fprintf(f, "%f ", v);
	}
	fprintf(f, "\n");
	// Points
	fprintf(f, "POINT_DATA %d\n", N);

	// Scalars
	fprintf(f, "SCALARS potential float 1\n"); 

	// Lookup Table
	fprintf(f, "LOOKUP_TABLE default\n"); 

	for (auto it = Values.begin(); it != Values.end(); ++it){
		double val = (*it).value;

		if (val > 0.0){
			if (val > 1.0e09)
				val = 1.0e09;
		}
		else{
			if (val < -1000000000.00)
				val = -1000000000.00;
		}


		fprintf(f, "%f\n", val);
		
	}

	return;

}

double GetTheta(Particle& i, Particle& j, Particle& k){

	double theta = 0.0;
	// Distances
	std::valarray<double> ij = j.x - i.x;
	std::valarray<double> ik = k.x - i.x;

	// Dot product distances
	double dot_p = 0.0;
	for (int d = 0; d < DIM; d++)
		dot_p += ij[d] * ik[d];

	double norm2_ij = 0.0, norm2_ik = 0.0;

	for (int d = 0; d < DIM; d++){

		//norm2_ij += sqr(ij[d]);
		norm2_ij += ij[d] * ij[d];
		//norm2_ik += sqr(ik[d]);
		norm2_ik += ik[d] * ik[d];
	}

	// Norms and squared norms
	double norm_ij = 0.0, norm_ik = 0.0;

	norm_ij = sqrt(norm2_ij);
	norm_ik = sqrt(norm2_ik);


	double res = dot_p / (norm_ij*norm_ik);

	/*for (int d = 0; d < DIM; d++){

		ij[d] /= norm_ij;
		ik[d] /= norm_ik;
		}

		double res = 0.0;

		for (int d = 0; d < DIM; d++){
		res += ij[d] * ik[d];
		}*/

	// Angle value in radiants
	if ((res < -1.0) && ((res + 1.0) < TOLL))
		res = -1.0;
	if ((res > 1.0) && ((res - 1.0) < TOLL))
		res = 1.0;


	theta = acos(res);

	if (theta * (180 / PI) > 180)
		printf("ERROR, INVALID THETA: %f, THETA:%f\n", theta * 180.0 / PI, theta);

	return theta;
}

double GetLambda(Particle& i, Particle& j, Particle& k){
	double lambda = 0.0;

	double lambda_ij = GetLambda(i, j);
	double lambda_jk = GetLambda(j, k);

	lambda = sqrt(lambda_ij*lambda_jk);

	return lambda;
}

double GetLambda(Particle& i, Particle& j){

	double lambda = 0.0;

	lambda = sqrt(i.lambda*j.lambda);

	return lambda;

}
double GetSigma(Particle& i, Particle& j){

	double sigma = 0.0;

	sigma = (i.sigma+j.sigma)*0.5;

	return sigma;

}

double GetEps(Particle& i, Particle& j, Particle& k){

	double eps_ij = GetEps(i, j);
	double eps_jk = GetEps(j, k);
	
	return sqrt(eps_ij*eps_jk);

}

double GetEps(Particle& i, Particle& j){

	double eps = 0.0;

	eps = sqrt(i.epsilon * j.epsilon);

	return eps;


}

double Distance(Particle& i, Particle& j){

	double dist2 = 0.0; double dist = 0.0;

	for (int d = 0; d < DIM; d++){
		dist2 += sqr(i.x[d] - j.x[d]);
	}

	dist = sqrt(dist2);

	return dist;


}

double U2(Particle& i, Particle& j, double A, double B, double p, double q, double a){

	double res = 0.0;

	double eps_ij = GetEps(i, j);
	double sigma_ij = GetSigma(i, j);
	double r_ij = Distance(i, j);

	double s = r_ij / sigma_ij;
	if (s < a){

		double s1 = pow(s, -p);
		double s2 = pow(s, -q);

		double exp_arg = pow((s - a), -1.0);
		double sigma_diff = B*s1 - s2;

		res = eps_ij*A*(sigma_diff)*exp(exp_arg);

		return res;
	}
	else
		return 0.0;

}

double U3(Particle& i, Particle& j, Particle& k, double Gamma, double a){

	double res = 0.0;

	//20.280000 20.240000 20.080000

	// Epsilon
	double eijk = GetEps(i, j, k);
	double ejik = GetEps(j, i, k);
	double ekij = GetEps(k, i, j);

	// Lambda
	double lambdaijk = GetLambda(i, j, k);
	double lambdajik = GetLambda(j, i, k);
	double lambdakij = GetLambda(k, i, j);

	// dist / sigma
	double dij = Distance(i, j);
	double dik = Distance(i, k);
	double sigmaij = GetSigma(i, j);
	double sigmaik = GetSigma(i, k);
	double sij = dij / sigmaij;
	double sik = dik / sigmaik;

	double dji = Distance(j, i);
	double djk = Distance(j, k);
	double sigmaji = GetSigma(j, i);
	double sigmajk = GetSigma(j, k);
	double sji = dji / sigmaji;
	double sjk = djk / sigmajk;

	double dki = Distance(k, i);
	double dkj = Distance(k, j);
	double sigmaki = GetSigma(k, i);
	double sigmakj = GetSigma(k, j);
	double ski = dki / sigmaki;
	double skj = dkj / sigmakj;

	// Theta Angle (first index = centre of angle)
	double theta_i, theta_j, theta_k;
	theta_i = GetTheta(i, j, k);
	theta_j = GetTheta(j, i, k);
	theta_k = GetTheta(k, i, j);

	// Sum terms
	double hijk, hjik, hkij;
	(sij < a && sik < a) ? hijk = h(sij, sik, theta_i, lambdaijk, Gamma, a) : hijk = 0.0;
	(sji < a && sjk < a) ? hjik = h(sji, sjk, theta_j, lambdajik, Gamma, a) : hjik = 0.0;
	(ski < a && skj < a) ? hkij = h(ski, skj, theta_k, lambdakij, Gamma, a) : hkij = 0.0;

	// Potential Value
	res = eijk*hijk + ejik*hjik + ekij*hkij;
	//res = eijk*hijk; 
	return res;
}

double h(double s_ab, double s_bc, double theta_b, double lambda, double gamma, double a){

	double h = 0.0;

	double const_angl = -0.333333333;

	double cos_term = cos(theta_b) + const_angl;

	double exp_arg = gamma*(pow((s_ab - a), -1.0) + pow((s_bc - a), -1.0));
	
	h = lambda*exp(exp_arg)*sqr(cos_term);

	return h;

}