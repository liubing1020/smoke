#ifndef ENZSUBS_ODE_H
#define ENZSUBS_ODE_H

#include <boost/array.hpp>

typedef boost::array<double,3> param_type;
typedef boost::array<double,4> state_type;

struct enzsubs_ode
{
	const param_type &p;
	
	enzsubs_ode(const param_type &par) : p(par){}
	
	void operator()(const state_type &x, state_type &dxdt, const double ) const {
		dxdt[0] = -p[0]*x[0]*x[2]+(p[1]+p[2])*x[1];
		dxdt[1] = +p[0]*x[0]*x[2]-(p[1]+p[2])*x[1];
		dxdt[2] = -p[0]*x[0]*x[2] + p[1]*x[1];
		dxdt[3] = p[2]*x[1];
		}
	
};

#endif