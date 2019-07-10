function seg_clock = seg_clock_pathway()
seg_clock = pathway();
seg_clock.data = [];
seg_clock.data_times = [];
seg_clock.meas_fun = @(x) x;
seg_clock.x_names = { 'Notch protein' 'cytosolic NicD' 'nuclear NicD' 'Lunatic fringe mRNA' 'Lunatic Fringe protein' 'phosph. beta-catenin' 'nuclear beta-catenin' 'Axin2 protein' 'Gsk3' 'beta-catenin' 'Axin2 mRNA' 'active Ras' 'active ERK' 'active TF X' 'Dusp6 mRNA' 'Dusp6 protein'};
seg_clock.x_init = [ 0.5 0.2 0 0.1 0.001 0.1 0.001 0.1 3 0.1 0.1 0.5 0.2 0.1 0.1 0.1];
seg_clock.x_init_upper = [0.5	2	0.12	4	4	0.1	1	1	3.2	2	20	2.7	2.7	2.7	10	16];
seg_clock.x_init_lower = 0*ones(1,length(seg_clock.x_init));
seg_clock.x_init_estimate = [];
seg_clock.p_nominal = [ 1.4 0.23 2.82 0.001 0.01 0.1 0.1 0.001 0.1 0.768 2.5 0.0000 3 1.92 0.05 0.37 0.39 0.3 0 0.087 7.062 0.06 1.64 0.8 0.7 0.05 0.48 2 2 0.5 0.02 0.6 0.63 0.1 1.8 0.28 0.03 0.7 1.5 0.5 2 0.5 1.35 0.5 0.103 0.1 0.05 0.05 0.05 0.5 0.5 0.05 0.5 0.5 2 2 3.45 2 2 0.3 1.5 0.3 0.9 0.5 5.08 1 4.968 0.41 3.3 1.6 0.5];
%parameter ranges as used in LiuBings TCS Paper 0.9-1.1*nomial_value of parameters
%seg_clock.p_lim_upper = [1.54	0.253	3.102	0.001	0.01	0.11	0.11	0.0011	0.11	0.8448	2.75	0	3.3	2.112	0.05	0.407	0.429	0.33	0	0.087	7.7682	0.06	1.804	0.8	0.77	0.05	0.48	2	2	0.55	0.022	0.66	0.693	0.1	1.8	0.28	0.03	0.77	1.65	0.55	2.2	0.55	1.485	0.55	0.1133	0.11	0.05	0.05	0.05	0.5	0.5	0.05	0.55	0.55	2	2	3.45	2	2	0.3	1.5	0.3	0.99	0.55	5.08	1	5.4648	0.451	3.63	1.76	0.55];
%seg_clock.p_lim_lower = [1.26	0.207	2.538	0.001	0.01	0.09	0.09	0.0009	0.09	0.6912	2.25	0	2.7	1.728	0.05	0.333	0.351	0.27	0	0.087	6.3558	0.06	1.476	0.8	0.63	0.05	0.48	2	2	0.45	0.018	0.54	0.567	0.1	1.8	0.28	0.03	0.63	1.35	0.45	1.8	0.45	1.215	0.45	0.0927	0.09	0.05	0.05	0.05	0.5	0.5	0.05	0.45	0.45	2	2	3.45	2	2	0.3	1.5	0.3	0.81	0.45	5.08	1	4.4712	0.369	2.97	1.44	0.45];
seg_clock.p_lim_upper=1.1*seg_clock.p_nominal;
seg_clock.p_lim_lower=0.9*seg_clock.p_nominal;
seg_clock.p_estimate = [ 1	2	3	6	7	8	9	10	11	13	14	16	17	18	21	23	25	30	31	32	33	38	39	40	41	42	43	44	45	46	53	54	63	64	67	68	69	70	71 ];
seg_clock.sim_reltol = 1e-2;
seg_clock.sim_abstol = 1e-5;
seg_clock.sim_fun = @seg_clock_ode_mex;
end