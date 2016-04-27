#ifndef UTILS_H
#define UTILS_H

typedef struct{
	int i;
	int j;
	int k;
	double theta;
	double val;

} u3_comp;

typedef struct{
	int i;
	int j;
	double val;

} u2_comp;

typedef struct{
	double x1;
	double x2;
	double x3;
	double d1;
	double d2;
	double value;

}pot_sample;



#define DIM 3

#define PI 3.14159265

#define sqr(x) ((x)*(x))
const double TOLL = 0.0000001;

#endif