#ifndef _ODEFUNDEF
#define _ODEFUNDEF
#include <iostream>
#include <cmath>
#include <nvector/nvector_serial.h>
#define N_SPECIES 6
#define N_PARAMS 21

#define MAX_CONV_FAIL 1000000
#define MAX_STEPS 1000000000
#define MAX_ERRFAILS 15
#define MIN_STEPSIZE 0.000000000000010000000000000000
#define MAX_STEPSIZE 100000000000000.00000

int odefun(double t, N_Vector x_in, N_Vector dx_in, void *f_data){
	double *p = (double*) f_data;
	double *x, *dx;
	x = NV_DATA_S(x_in);
	dx = NV_DATA_S(dx_in);

	dx[0] = -p[0]*x[0] + p[6]/(p[7]+p[8]*pow(x[4],p[9]));
	dx[1] = -p[1]*x[1] + p[10]/(p[11]+p[12]*pow(x[5],p[13]));
	dx[2] = -p[2]*x[2] + p[14]/(p[15]+p[16]*pow(x[3],p[17]));
	dx[3] = p[3]*x[0] -p[18]*x[3];
	dx[4] = p[4]*x[1] -p[19]*x[4];
	dx[5] = p[5]*x[2] -p[20]*x[5];
	
	return 0;
};

#endif