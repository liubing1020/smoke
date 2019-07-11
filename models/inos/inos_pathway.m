function inos = inos_pathway()
inos = pathway();
inos.data = [];
inos.data_times = [];
inos.meas_fun = @(x) x;
inos.x_names= {'NO' 'OOH-AA-PE' 'GSSG' '15LOXFe3' 'AA' 'OH-AA-PE' 'PEoxtr' 'AA-PE' 'RSL3' 'iNOS' 'GPX4' 'NO-PE' 'NO-PEoxtr' 'ACSL4/LPCAT3' 'SNAP' 'L-NIL' 'death'   };
inos.x_init = [0 0 1.045e-06 4e-12 1e-12 0 0 1e-15 5e-10 1e-11 6e-11 1e-12 0 6.8e-13 0 0 0 ];
inos.x_init_upper = 1000*ones(1,length(inos.x_init));
inos.x_init_lower = 0*ones(1,length(inos.x_init));
inos.x_init_estimate = [];
inos.p_nominal = [1e-12 9.24027 0.001 756.429 8599.98 1.85041e+06 17748.5 115306 0.1161 9514.37 0.1 1146.67 12978.2 1 10 1000 52.3 0.098 0.1 50 6.58678e-05 0.0012335 0.01 0.01 0.1 0.1 0.1 ];
inos.p_lim_upper = 2*inos.p_nominal;
inos.p_lim_lower = 0.5*inos.p_nominal;

inos.p_estimate = [13 14 15 16 17 18 19 20 21 24 25 34];
inos.sim_reltol = 1e-2;
inos.sim_abstol = 1e-5;
inos.sim_fun = @inos_ode_mex;
end
