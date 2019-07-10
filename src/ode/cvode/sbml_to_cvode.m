function sbml_to_cvode(sbml_file, name_prefix, type)
	if(nargin < 3)
		type = 'ode';
	end
	
    if(isempty(sbml_file))
        [sbml_file,sbml_path] = uigetfile('*.xml'); 
        model = TranslateSBML([sbml_path sbml_file]);
    else
        model = TranslateSBML(sbml_file);
    end
    
    
    ns = length(model.species);
    nr = length(model.reaction);
    species = containers.Map;
    
    reaction = cell(nr,1);
    ncp = length(model.compartment);
    compartments = containers.Map;
    
    x0_string = '';
    xnames_string = '';
    for i=1:ns
       species(model.species(i).id) = i; 
       x0_string = [x0_string ' ' num2str(model.species(i).initialConcentration)];
       if(isempty(model.species(i).name))
           xnames_string = [xnames_string ' ''' model.species(i).id ''''];
       else
           xnames_string = [xnames_string ' ''' model.species(i).name ''''];
       end
    end

    for i=1:ncp
        compartments(model.compartment(i).id) = model.compartment(i).size; 
    end

    nparams = 0;
    parameters = containers.Map;
    param_values = '';
    % Global parameters
    params = model.parameter;
    for i=1:length(params)
       parameters(params(i).id) = nparams;
       nparams = nparams + 1;
       param_values = [param_values ' ' num2str(params(i).value)];
    end
    % Local parameters from inside reactions
    for i=1:nr
        params = model.reaction(i).kineticLaw.parameter;
        for j=1:length(params)
            parameters([model.reaction(i).id '|' params(j).id]) = nparams;
            nparams = nparams + 1;
            param_values = [param_values ' ' num2str(model.reaction(i).kineticLaw.parameter(j).value)];
        end
    end
    
    species_reactions = zeros(ns,nr);
    for i=1:nr
        reactants = model.reaction(i).reactant;
        products = model.reaction(i).product;
       [ibr,ibp] = is_bothsides(reactants, products);
       for j = 1:length(reactants)
          if(ibr(j)==0) 
              species_reactions(species(reactants(j).species),i) = -1;
          end
       end
       for j = 1:length(products)
          if(ibp(j)==0) 
              species_reactions(species(products(j).species),i) = 1;
          end
       end       
    end
    
     const_species = find(sum(abs(species_reactions),2)==0);
%     species_partners = zeros(ns,1);
%     species_partners(const_species) = -1;
%     for i=1:ns-1
%         if(species_partners(i)==-1), continue; end
%         for j=i+1:ns
%             if(species_partners(j)==-1), continue; end
%             if sum(abs(species_reactions(i,:)+species_reactions(j,:)))==0
%                 species_partners(j) = i;
%             end
%         end
%     end
% 
%     indep_species = find(species_partners==0);    
%     n_indep = length(indep_species);
    
%     init_str = '';    
%    
%     for i=1:ns
%         if (species_partners(i) == -1)
%             init_str = [init_str 'x[' num2str(i-1) '] = x0[' num2str(i-1) '];\n' ];
%         elseif (species_partners(i) == 0)
%             ode_id = find(indep_species==i);
%             init_str = [init_str 'x[' num2str(i-1) '] = x_in[' num2str(ode_id-1) '];\n'];
%         else
%             partner_ode_id = find(indep_species==species_partners(i));
%             init_str = [init_str 'x[' num2str(i-1) '] = x0[' num2str(i-1) '] + x0[' num2str(species_partners(i)-1) '] - x_in[' num2str(partner_ode_id-1) '];\n'];
%         end
%     end

%     J = sparse(n_indep,n_indep);    
    
    odes = cell(1,ns);
    for i=1:ns
        odes{i} = ['dx[' num2str(i-1) '] = '];
		if any(const_species==i)
			odes{i} = [odes{i} '0'];
		end
	end
	
    nr_total = 0;
    for i=1:nr
		reaction{i} = model.reaction(i).kineticLaw.math;
        reaction{i} = strrep(reaction{i}, '(', ' ( ');
        reaction{i} = strrep(reaction{i}, ')', ' ) ');
        reaction{i} = strrep(reaction{i}, ',', ' , ');
        splreact = regexp(reaction{i}, ' ', 'split');
        rate_species = [];
        for j=1:length(splreact)
            splreact_name = splreact{j};
            if(isKey(species,splreact_name))
                id = species(splreact_name);
                rate_species = [rate_species, id];
                splreact{j} = ['x[' num2str(id-1) ']'];
            elseif(isKey(compartments,splreact{j}))
                vol = compartments(splreact_name);
                splreact{j} = num2str(vol);
            elseif(isKey(parameters,splreact_name))
                id = parameters(splreact_name);
                splreact{j} = ['p[' num2str(id) ']'];
            elseif(isKey(parameters,[model.reaction(i).id '|' splreact_name]))
                id = parameters([model.reaction(i).id '|' splreact_name]);
                splreact{j} = ['p[' num2str(id) ']'];
            end
		end
		
		react_str = strcat(splreact{:});
		
		%%% FOR GILLESPIE %%%%%%%%
		nr_total = nr_total+1;
        
		if(model.reaction(i).reversible)
			splreact = regexp(react_str, '-', 'split');
			if(length(splreact)>2)
				warning('Reversible reaction has more than 2 terms.');
			end
			reactionStr{nr_total} = ['r[' num2str(nr_total-1) '] = ' par_rebalance(splreact{1}) ';\n'];
			nr_total = nr_total+1;
			reactionStr{nr_total} = ['r[' num2str(nr_total-1) '] = ' par_rebalance(splreact{2}) ';\n'];
		else
			reactionStr{nr_total} = ['r[' num2str(nr_total-1) '] = ' react_str ';\n'];
		end
		
		
		%%%%%%%%%%%%%%%%%%%%%%%%%
		
        [ibr,ibp] = is_bothsides(model.reaction(i).reactant, model.reaction(i).product);
        
		%%%% FOR GILLESPIE %%%%%%
		D(nr_total,:) = zeros(1,ns);
		for j=1:length(model.reaction(i).reactant)
			if(ibr(j)==0)
				id = species(model.reaction(i).reactant(j).species);
				D(nr_total,id) = -1;
			end
		end
		for j=1:length(model.reaction(i).product)
			if(ibp(j)==0)
				id =  species(model.reaction(i).product(j).species);
				D(nr_total,id) = 1;
			end
		end
		
		if(model.reaction(i).reversible)
			D(nr_total-1,:) = D(nr_total,:);
			D(nr_total,:) = -D(nr_total,:);
		end
		%%%%%%%%%%%%%%%%%%%%%%%%
		
		%J = setjac(J,species,model.reaction(i).reactant, model.reaction(i).product ,ibr,ibp, rate_species, species_partners);
        for j=1:length(model.reaction(i).reactant)
           if(ibr(j)==0)
            id = species(model.reaction(i).reactant(j).species);
            ode_id = id;
            if(~isempty(ode_id))
                odes{ode_id} = [odes{ode_id} '-(' react_str ') '];
            end
           end
        end
        for j=1:length(model.reaction(i).product)
           if(ibp(j)==0)
            id =  species(model.reaction(i).product(j).species);
            ode_id = id;
            if(~isempty(ode_id))
                odes{ode_id} = [odes{ode_id} '+(' react_str ') '];
            end
           end
        end
    end

    for i=1:ns
       odes{i} = [ odes{i} ';\n']; 
    end
    
%    [i,j,s] = find(J);
%    i_str = ['i = [' num2str(i') '];\n'];
%    j_str = ['j = [' num2str(j') '];\n'];
%    jac_str = [i_str j_str name_prefix '.sim_jac_pattern = sparse(i,j,1);\n'];
    
	if(strcmp(type,'ode'))
		% Pathway file
		pathway_file = [name_prefix '_pathway.m'];
		fh = fopen(pathway_file,'w+t');
		fprintf(fh, ['function ' name_prefix ' = ' name_prefix '_pathway()\n']);
		fprintf(fh, [name_prefix ' = pathway();\n' ]);
		fprintf(fh, [name_prefix '.data = [];\n']);
		fprintf(fh, [name_prefix '.data_times = [];\n']);
		fprintf(fh, [name_prefix '.meas_fun = @(x) x;\n']);
		fprintf(fh, [name_prefix '.x_names = {' xnames_string '};\n']);
		fprintf(fh, [name_prefix '.x_init = [' x0_string '];\n']);
		fprintf(fh, [name_prefix '.x_init_upper = 1000*ones(1,length(' name_prefix '.x_init));\n']);
		fprintf(fh, [name_prefix '.x_init_lower = 0*ones(1,length(' name_prefix '.x_init));\n']);
		fprintf(fh, [name_prefix '.x_init_estimate = [];\n']);
		fprintf(fh, [name_prefix '.p_nominal = [' param_values '];\n']);
		fprintf(fh, [name_prefix '.p_lim_upper = 1000*ones(1,length(' name_prefix '.p_nominal));\n']);
		fprintf(fh, [name_prefix '.p_lim_lower = 0*ones(1,length(' name_prefix '.p_nominal));\n']);
		fprintf(fh, [name_prefix '.p_estimate = 1:length(' name_prefix '.p_nominal);\n']);
		fprintf(fh, [name_prefix '.sim_reltol = 1e-2;\n']);
		fprintf(fh, [name_prefix '.sim_abstol = 1e-5;\n']);
		fprintf(fh, [name_prefix '.sim_fun = @getp;\n']);
		fprintf(fh, [name_prefix '.sim_method = @' name_prefix '_mex_wrapper;\n']);
		fprintf(fh, [name_prefix '.species_dep = zeros(size(' name_prefix '.x_init));\n']);
		fprintf(fh, 'end');
		fclose(fh);
    
		% Wrapper file
		wrapper_file = [name_prefix '_mex_wrapper.m'];
		fh = fopen(wrapper_file,'w+t');
		fprintf(fh,['function [t x] = ' name_prefix '_mex_wrapper(sim_fun,ts,x0,opts)\n']);
		fprintf(fh,'p = sim_fun(1,1);\n');
		fprintf(fh,['try\n [t x] = ' name_prefix '_mex(p, x0, ts);\n']);
		fprintf(fh,'catch\n	t = [];x = [];\n end\n');
		fprintf(fh,'end\n');
		fclose(fh);
		
		
		% ODE function file
		ode_file = [name_prefix '_ode.c'];
		fh = fopen(ode_file,'w+t');
		fprintf(fh,'#ifndef _ODEFUNDEF\n');
		fprintf(fh,'#define _ODEFUNDEF\n');
		fprintf(fh,'#include <iostream>\n');
		fprintf(fh,'#include <nvector/nvector_serial.h>\n');
		fprintf(fh,['#define N_SPECIES ' num2str(ns) '\n']);
		fprintf(fh,['#define N_PARAMS ' num2str(nparams) '\n']);
		fprintf(fh,'#define MAX_CONV_FAIL 1000000\n');
		fprintf(fh,'#define MAX_STEPS 1000000000\n');
		fprintf(fh,'#define MAX_ERRFAILS 15\n');
		fprintf(fh,'#define MIN_STEPSIZE 0.000000000000010000000000000000\n');
		fprintf(fh,'#define MAX_STEPSIZE 100000000000000.00000\n');
		fprintf(fh,'int odefun(double t, N_Vector x_in, N_Vector dx_in, void *f_data){\n');
		fprintf(fh,'double *p = (double*) f_data;\n');
		fprintf(fh,'double *x, *dx;\n');
		fprintf(fh,'x = NV_DATA_S(x_in);\n');
		fprintf(fh,'dx = NV_DATA_S(dx_in);\n');
		fprintf(fh, strcat(odes{:}));
		fprintf(fh,'return 0;};\n');
		fprintf(fh,'#endif\n');
		fclose(fh);
		
		% MEX wrapper for CVODE simulaton
		mex_file = [name_prefix '_mex.cpp'];
		fh = fopen(mex_file,'w+t');
		fprintf(fh,['#include "' ode_file '"\n']);
		fprintf(fh,'#include "mex.h"\n');
		fprintf(fh,'#include "cvode_sim.c"\n');
		fprintf(fh,'using namespace std;\n');
		fprintf(fh,'extern void _main();\n');
		fprintf(fh,'void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){\n');
		fprintf(fh,'double *p_in, *x_in, *x_out, *x_tmp, *t_in, *t_out;\n');
		fprintf(fh,'int nx, nt;\n');
		fprintf(fh,'p_in = mxGetPr(prhs[0]);\n');	
		fprintf(fh,'x_in = mxGetPr(prhs[1]);\n');	
		fprintf(fh,'t_in = mxGetPr(prhs[2]);\n');	
		fprintf(fh,'nx = mxGetN(prhs[1]);\n');
		fprintf(fh,'nt = mxGetN(prhs[2]);\n');
		fprintf(fh,'x_tmp = (double*)malloc(sizeof(double)*nt*nx);\n');
		fprintf(fh,'cvode_sim(t_in,nt,x_in,p_in,odefun,x_tmp);\n');
		fprintf(fh,'plhs[0] = mxCreateDoubleMatrix(nt,1,mxREAL);\n');
		fprintf(fh,'plhs[1] = mxCreateDoubleMatrix(nt,nx,mxREAL);\n');
		fprintf(fh,'t_out = mxGetPr(plhs[0]);\n');
		fprintf(fh,'x_out = mxGetPr(plhs[1]);\n');
		fprintf(fh,'for(int i=0;i<nx;i++){\n');
		fprintf(fh,'for(int j=0;j<nt;j++){\n');
		fprintf(fh,'x_out[i*nt+j] = x_tmp[j*nx+i];\n');
		fprintf(fh,'}\n');
		fprintf(fh,'}\n');
		fprintf(fh,'memcpy(t_out,t_in,nt*sizeof(double));\n');
		fprintf(fh,'}\n');
		fclose(fh);
	end
	
	
	%%% FOR GILLESPIE %%%%%%%%%%%%
	if(strcmp(type,'gillespie'))
		rate_file = [name_prefix '_rates.h'];
		fh = fopen(rate_file,'w+t');
		D_str = num2str(reshape(D',1,numel(D)),'%g,');
		fprintf(fh,['double D[] = {' D_str(1:end-1) '};\n']);
		fprintf(fh,['int nr = ' num2str(nr_total) ';\n']);
		fprintf(fh,'void rhs(double* x, double* p, double* r){\n');
		fprintf(fh,strcat(reactionStr{:}));
		fprintf(fh,'}\n');
	end
	%%%%%%%%%%%%%
end

function sOut = par_rebalance(s)
	sOut = s;
	nOpen = length(find(s=='('));
	nClose = length(find(s==')'));
	while (nOpen < nClose)
		sOut = ['(' sOut];
		nOpen = nOpen + 1;
	end
	while (nOpen > nClose)
		sOut = [sOut ')'];
		nClose = nClose + 1;
	end
end

function [ibr,ibp] = is_bothsides(reactants, products)
    ibr = zeros(length(reactants),1);
    ibp = zeros(length(products),1);
    for i=1:length(reactants)
        for j=1:length(products)
            if(strcmp(reactants(i).species,products(j).species))
                ibr(i) = 1;
                ibp(j) = 1;
            end
        end
    end
end

function J = setjac(J_in,species,reactants, products ,ibr,ibp,rate_species, species_partners)
    J = J_in;
    indep_species = find(species_partners==0);
    for i=1:length(rate_species)                                % For all species affecting rate
        if(species_partners(rate_species(i))==-1), continue;                      % If it's constant skip
        elseif (species_partners(rate_species(i))==0)
            ratesp_ode_id = find(indep_species==rate_species(i));
        else
            ratesp_ode_id = find(indep_species==species_partners(rate_species(i)));
        end
        if(isempty(ratesp_ode_id)), continue; end
        for j=1:length(products)
            if ibp(j)==0
               id = species(products(j).species);
               ode_id = find(indep_species==id,1);
                if(~isempty(ode_id))               
                    J(ode_id,ratesp_ode_id) = 1;
                end
            end
        end
        for j=1:length(reactants)
            if ibr(j)==0
               id = species(reactants(j).species);
               ode_id = find(indep_species==id,1);
               if(~isempty(ode_id))                              
                    J(ode_id,ratesp_ode_id) = 1;
               end
            end
        end
    end
end