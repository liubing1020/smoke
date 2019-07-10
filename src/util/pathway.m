classdef pathway < handle
    properties
        data
        data_times
        data_rel
        meas_fun
        loglh_fun
        sim_jac_pattern
        sim_reltol
        sim_abstol
        sim_fun
        sim_method
        sim_opts
        x_names
        x_init
        x_init_upper
        x_init_lower
        x_init_estimate
        p_nominal
        p_lim_upper
        p_lim_lower
        p_estimate
        species_dep
    end
    
    methods
        function obj=pathway()
        end
    end
end