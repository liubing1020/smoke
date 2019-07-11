#ifndef _ODEFUNDEF
#define _ODEFUNDEF
#include <iostream>
#include <cmath>
#include <nvector/nvector_serial.h>
#define N_SPECIES 17
#define N_PARAMS 27
#define N_ASSIGNS 3
#define N_CTS 2

#define MAX_CONV_FAIL 1000000
#define MAX_STEPS 10000
#define MAX_ERRFAILS 15
#define MIN_STEPSIZE 0.000000000000010000000000000000
#define MAX_STEPSIZE 100000000000000.00000

double MM2_1_(double param_0, double modif_0, double sub_0, double sub_1, double param_1, double param_2) 	//MM2 [1]
{return  param_0*modif_0*sub_0*sub_1/(param_1*param_2+sub_0*param_1+sub_1*param_2+sub_0*sub_1);} 
double Henri_Michaelis_Menten_irreversible_(double sub_0, double param_0, double param_1) 	//Henri-Michaelis-Menten (irreversible)
{return  param_1*sub_0/(param_0+sub_0);} 
double MM1(double param_0, double modif_0, double sub_0, double param_1) 	//MM1
{return  param_0*modif_0*sub_0/(param_1+sub_0);} 
double ConstantFlux_irreversible_(double param_0) 	//Constant flux (irreversible)
{return  param_0;} 

int odefun(double t, N_Vector x_in, N_Vector dx_in, void *f_data){
	double *p = (double*) f_data;
	double *x, *dx, *y, *ct;
	x = NV_DATA_S(x_in);
	dx = NV_DATA_S(dx_in);
	y = (double *) malloc(N_ASSIGNS * sizeof(double));
	ct = (double *) malloc(N_CTS * sizeof(double));

	ct[0] = 2.0899999999999995e-06;	//ct[0] conserved total for 'GSH'
	ct[1] = 3.9999999999999991e-12;	//ct[1] conserved total for '15LOXFe2'

	y[0] = ct[0]-x[2];	//metabolite 'GSH': reactions
	y[1] = ct[1]-x[3];	//metabolite '15LOXFe2': reactions
	y[2] = x[16]*9.24027/(T+1.00000000000000000) * p[0];	//model entity 'measurement_death':assignment

	dx[0] = (p[22] * x[9]) *p[0]-(p[23] * x[0]) *p[0]-(p[2] * x[6] * x[0]) *p[0]-(p[2] * x[1] * x[0]) *p[0]+(p[25] * x[14]) *p[0];
	dx[1] = -MM2_1_(p[3], x[10], y[0], x[1], p[4], p[5])*p[0]+MM1(p[8], x[3], x[7], p[9])*p[0]-(p[24] * x[1]) *p[0]-(p[2] * x[1] * x[0]) *p[0];
	dx[2] = MM2_1_(p[3], x[10], y[0], x[1], p[4], p[5])*p[0]-Henri_Michaelis_Menten_irreversible_(x[2], p[6], p[7])*p[0]+MM1(p[11], x[10], y[0], p[12])*p[0];
	dx[3] = Henri_Michaelis_Menten_irreversible_(y[1], p[16], p[17])*p[0]-MM1(p[18], x[0], x[3], p[19])*p[0];
	dx[4] = ConstantFlux_irreversible_(p[13])*p[0]-MM1(p[14], x[13], x[4], p[15])*p[0];
	dx[5] = MM2_1_(p[3], x[10], y[0], x[1], p[4], p[5])*p[0]-(p[10] * x[5]) *p[0];
	dx[6] = (p[24] * x[1]) *p[0]-(p[2] * x[6] * x[0]) *p[0];
	dx[7] = -MM1(p[8], x[3], x[7], p[9])*p[0]+MM1(p[14], x[13], x[4], p[15])*p[0];
	dx[8] = -(p[21] * x[8]) *p[0];
	dx[9] = -(p[26] * x[9] * x[15]) *p[0];
	dx[10] = -(p[20] * x[8] * x[10]) *p[0];
	dx[11] = (p[2] * x[1] * x[0]) *p[0];
	dx[12] = (p[2] * x[6] * x[0]) *p[0];
	dx[13] = 0;
	dx[14] = 0;
	dx[15] = 0;
	dx[16] = x[6]+x[1];	//model entity 'death':ode

	free(y);
	free(ct);

	return 0;
};

#endif
