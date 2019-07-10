#ifndef _ODEFUNDEF
#define _ODEFUNDEF
#include <iostream>
#include <nvector/nvector_serial.h>
#define N_SPECIES 16
#define N_PARAMS 71

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

	dx[0] = +(1*p[59]*p[1]) -(p[59]*1*p[2]*x[0]/(p[0]+x[0])) -(p[59]*1*p[56]*x[0]*pow(p[49],p[57])/(pow(p[49],p[57])+pow(x[4],p[57]))) ;
	dx[1] = +(p[59]*1*p[56]*x[0]*pow(p[49],p[57])/(pow(p[49],p[57])+pow(x[4],p[57]))) -(p[59]*1*p[4]*x[1]/(p[3]+x[1])) -(p[59]*1*(p[5]*x[1]-p[6]*x[2])) ;
	dx[2] = +(p[59]*1*(p[5]*x[1]-p[6]*x[2])) -(p[59]*1*p[8]*x[2]/(p[7]+x[2])) ;
	dx[3] = +(p[59]*1*p[12]*(p[10]/(p[10]+x[8]))*pow(x[2],p[58])/(pow(p[14],p[58])+pow(x[2],p[58]))) -(p[59]*1*p[13]*x[3]/(p[9]+x[3])) ;
	dx[4] = +(p[59]*1*p[17]*x[3]) -(p[59]*1*p[16]*x[4]/(p[15]+x[4])) ;
	dx[5] = +(p[60]*1*p[64]*p[50]/(p[50]+2)*x[9]/(p[35]+x[9])*(3-x[8])/3) -(p[60]*1*p[65]*x[5]/(p[36]+x[5])) -(p[60]*1*p[20]*x[5]) ;
	dx[6] = -(p[60]*1*(p[38]*x[6]-p[37]*x[9])) ;
	dx[7] = +(p[60]*1*(p[33]*(3-x[8])-p[34]*x[7]*x[8])) +(p[60]*1*p[30]*x[10]) -(p[60]*1*p[31]*x[7]/(p[32]+x[7])) ;
	dx[8] = +(p[60]*1*(p[33]*(3-x[8])-p[34]*x[7]*x[8])) ;
	dx[9] = +(p[60]*1*p[19]) -(p[60]*1*p[18]*x[9]) -(p[60]*1*p[64]*p[50]/(p[50]+2)*x[9]/(p[35]+x[9])*(3-x[8])/3) +(p[60]*1*p[65]*x[5]/(p[36]+x[5])) +(p[60]*1*(p[38]*x[6]-p[37]*x[9])) ;
	dx[10] = +(p[60]*1*p[21]) +(p[60]*1*(p[22]*pow(x[6],p[27])/(pow(p[24],p[27])+pow(x[6],p[27])))) +(p[60]*1*(p[29]*pow(x[13],p[28])/(pow(p[25],p[28])+pow(x[13],p[28])))) -(p[60]*1*p[23]*x[10]/(p[26]+x[10])) ;
	dx[11] = +(p[61]*1*p[66]*pow(1,p[55])/(pow(p[43],p[55])+pow(1,p[55]))*(2-x[11])/(p[44]+(2-x[11]))) -(p[61]*1*p[67]*x[11]/(p[45]+x[11])) ;
	dx[12] = +(p[61]*1*p[68]*x[11]/2*(2-x[12])/(p[47]+(2-x[12]))) -(p[61]*1*p[42]*x[15]*x[12]/(p[46]+x[12])) ;
	dx[13] = +(p[61]*1*p[69]*x[12]/2*(2-x[13])/(p[48]+(2-x[13]))) -(p[61]*1*p[70]*x[13]/(p[51]+x[13])) ;
	dx[14] = +(p[61]*1*p[62]*pow(x[13],p[54])/(pow(p[52],p[54])+pow(x[13],p[54]))) -(p[61]*1*p[63]*x[14]/(p[53]+x[14])) ;
	dx[15] = +(p[61]*1*p[39]*x[14]) -(p[61]*1*p[40]*x[15]/(p[41]+x[15])) ;
	
	return 0;
};

#endif