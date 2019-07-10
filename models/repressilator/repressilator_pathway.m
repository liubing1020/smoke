function represillator = represillator_pathway()
represillator = pathway();
represillator.data = [];
represillator.data_times = [];
represillator.meas_fun = @(x) x;
%represillator.x_names = { 'LacI protein' 'TetR protein' 'cI protein' 'LacI mRNA' 'TetR mRNA' 'cI mRNA'};
represillator.x_names = { 'm1' 'm2' 'm3' 'p1' 'p2' 'p3'};
represillator.x_init = [0 2 0 0 2 0];
represillator.x_init_estimate = [];
represillator.p_nominal = [ 50 50 50 1 1 1 100 1 8 3 100 1 8 3 100 1 8 3 1 1 1];
%use for parameter estimation
represillator.p_estimate = [1 2 3 7 11 15 9 13 17];
represillator.sim_reltol = 1e-2;
represillator.sim_abstol = 1e-5;
end