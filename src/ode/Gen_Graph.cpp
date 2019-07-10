#include "mex.h"
#include "matrix.h"
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

#ifdef MODEL_THROMBIN
	#include "../../models/thrombin/thrombin_ode.h"
#elif MODEL_SEGCLOCK
	#include "../../models/segmentationclock/seg_clock_ode.h"
#elif MODEL_EGFNGF
	#include "../../models/egf-ngf/egf_ngf_ode.h"
#elif MODEL_REPRESS
	#include "../../models/repressilator/repressilator_ode.h"
#elif MODEL_TLR
	#include "../../models/TLRpathway/TLR_ode.h"
#endif

#include "cvode_sim.c"

using namespace std;

extern void _main();

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	double *p_in;
	double *x_in;
	double *x_pre,*x_post;
	double dt;
	double pathway_identifier;
	int n_x_in ,n_formulas_in,i,n,n_p_in,n_elements_limits,nTimes;
	const mwSize *dims;
	mwSize ndim;

	// the argument list is in the order - parameters,variables,nTimes,the formulas to be checked, the discretization of the variables, dt
	p_in = mxGetPr(prhs[0]);
	x_in = mxGetPr(prhs[1]);
	n_p_in = mxGetN(prhs[0]);
	n_x_in = mxGetN(prhs[1]);
	nTimes= (int)mxGetScalar(prhs[2]);
	dt = mxGetScalar(prhs[3]);

	//convert the initial x values to vectors
	x_pre = new double[n_x_in];
	x_post = new double[2*n_x_in];
	
	memcpy(x_pre,x_in,n_x_in*sizeof(double));

	double ** Elements_Limit_Holder = new double*[nTimes];
	for (int i=0;i<nTimes;i++)
	{
		Elements_Limit_Holder[i] = new double[n_x_in];
	}
	
	double ts[2]={0.0,dt};
	double* lastState = x_pre; 
	for(int i=0;i<nTimes;i++)
	{

		for(int j=0;j< n_x_in;j++){
			Elements_Limit_Holder[i][j] = lastState[j];
		}
		int nStepsTaken = cvode_sim(ts,2,x_pre,p_in,odefun,x_post);
		lastState = x_post+n_x_in;
		memcpy(x_pre,lastState,n_x_in*sizeof(double));
		
	}

	plhs[0] = mxCreateDoubleMatrix(nTimes, n_x_in, mxREAL);
	double *candi = mxGetPr(plhs[0]);
	for(int i=0; i<nTimes; i++)
	{
			memcpy(&candi[n_x_in*i], Elements_Limit_Holder[i], sizeof(double)*n_x_in);
	}

	for (int i=0;i<nTimes;i++)
	{
		delete(Elements_Limit_Holder[i]);
	}
	
	delete(Elements_Limit_Holder);
	delete(x_pre);
	delete(x_post);


	return;
}
