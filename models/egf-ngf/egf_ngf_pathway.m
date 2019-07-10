function egf_ngf = egf_ngf_pathway()
egf_ngf = pathway();
egf_ngf.data = [];
egf_ngf.data_times = [];
egf_ngf.meas_fun = @(x) x;
egf_ngf.x_names = { 'EGF' 'NGF' 'freeEGFReceptor' 'boundEGFReceptor' 'freeNGFReceptor' 'boundNGFReceptor' 'SosInactive' 'SosActive' 'P90RskInactive' 'P90RskActive' 'RasInactive' 'RasActive' 'RasGapActive' 'Raf1Inactive' 'Raf1Active' 'BRafInactive' 'BRafActive' 'MekInactive' 'MekActive' 'ErkInactive' 'ErkActive' 'PI3KInactive' 'PI3KActive' 'AktInactive' 'AktActive' 'C3GInactive' 'C3GActive' 'Rap1Inactive' 'Rap1Active' 'RapGapActive' 'PP2AActive' 'Raf1PPtase'};
egf_ngf.x_init = [ 10002000 456000 80000 0 10000 0 120000 0 120000 0 120000 0 120000 120000 0 120000 0 600000 0 600000 0 120000 0 120000 0 120000 0 120000 0 120000 120000 120000];
egf_ngf.x_init_upper = 1000*ones(1,length(egf_ngf.x_init));
egf_ngf.x_init_lower = 0*ones(1,length(egf_ngf.x_init));
egf_ngf.x_init_estimate = [];
egf_ngf.p_nominal = [ 2.185e-005 0.012101 1.3821e-007 0.0072381 694.731 6086070 389.428 2112.66 1611.97 896896 32.344 35954.3 1509.36 1432410 0.8841 62464.6 185.759 4768350 125.089 157948 2.8324 518753 9.8537 1007340 8.8912 3496490 0.02137 763523 10.6737 184912 0.077107 272056 0.056628 653951 15.1212 119355 146.912 12876.2 1.4015 10965.6 27.265 295990 2.21 1025460 0.12633 1061.71 441.287 10879500];
%egf_ngf.p_lim_upper = [0.0000438	0.0242016	0.000000276	0.01447622	1389.462	12180000	778.856	4225.32	3223.94	1793792	64.688	71908.6	3018.72	2860000	1.768192	124929.2	371.518	9540000	250.178	315896	5.66486	1037506	19.70734	2020000	17.7824	7000000	0.0427394	1527046	21.3474	369824	0.1542134	544112	0.1132558	1307902	30.2424	238710	293.824	25752.4	2.8029	21931.2	54.53	591980	4.4199	2060000	0.252658	2123.42	882.574	21800000];
%egf_ngf.p_lim_lower = 0*ones(1,length(egf_ngf.p_nominal));

egf_ngf.p_lim_upper = 2*egf_ngf.p_nominal;
egf_ngf.p_lim_lower = 0.0*egf_ngf.p_nominal;

egf_ngf.p_estimate = [1 2 3 4 11 12 15 17 23 27 28 29 33 34 37 38 39 41 43 44];
egf_ngf.sim_reltol = 1e-2;
egf_ngf.sim_abstol = 1e-5;
egf_ngf.sim_fun = @egf_ngf_ode_mex;
end