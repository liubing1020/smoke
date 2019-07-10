%% Add paths
clear variables
addpath util/
addpath ode/
addpath analysis
addpath modelcheck
addpath ../models/egf-ngf/
addpath ../tools/sres/
diary('egf-ngf-parest-fast-5per-sres-0-10-09-notp.txt');

diary on;

tic;

%% Create the model
cd ode
compile('egfngf','mc')
cd ..
model = egf_ngf_pathway();
option = 'ParamEstimation';

%optMethod = 'sres';
optMethod = 'sres';

% Statistical model checking method
MCparams.method = 'YounesA';
MCparams.rules=[];

% Statistical model checking parameters
MCparams.alpha=0.05;
MCparams.beta=0.05;
MCparams.gamma=0.05;
MCparams.delta=0.05;

model.p_lim_lower = 0.0*model.p_nominal;
model.p_lim_upper = 10.0*model.p_nominal;


for i=1:length(model.p_nominal)
    if(~any(find(model.p_estimate==i)))
    model.p_lim_lower(i)=model.p_nominal(i);
    model.p_lim_upper(i)=model.p_nominal(i);
    end
end



%% Define all the rules here
%MCparams.rules{1} = '([27==1]&F{3}([27==6]))&((F{10}([27==6]&X(G([27==6]))))&(F{}([27==6]&X(G([27==6])))))';
%MCparams.rules{2} = '([25==1]&F{}([25==5]))&F{30}(!([25==5]))';
%MCparams.rules{3} = '([8==1]&F{2}([8==6]))&(F{}([8==6]&F[8==2]))';
%MCparams.rules{4} = '[10==1]& F{}([10==5])';
%MCparams.rules{5} = '([21==1]&F{10}([21==6]))&((F{10}([21==6]&X(G([21==6]))))&(F{}([21==6]&X(G([21==6])))))';
%MCparams.rules{6} = '([4==1]&F{2}([4==6]))&((F{10}([4==6]&X(G([4==6]))))&(F{}([4==6]&X(G([4==6])))))';
%MCparams.rules{7} = '([6==1]&F{10}(![6==5]))&(F{}([6==5]&X(G([6==5]))))';

%%Get rules from Data

%Generate one rule for each species.
%MCparams.rules = dataToIFormula('../../models/egf-ngf/training.csv',0.05);

%% Define probabilities for each rule

a=[];
len=0;
time=0;
[a,len,time] = dataToIFormula('../models/egf-ngf/training.csv',0.1,2);

if numel(MCparams.rules)==0
	MCparams.rules = a;
else
	MCparams.rules = [MCparams.rules,a];
end

MCparams.probs = 0.9*ones(size(MCparams.rules));
MCparams.probGreater = ones(size(MCparams.rules));



% Time between model checking steps (time discretization)
dt = 1;
% Number of discrete time steps
nTime = 61;

switch(option) 
	case 'ParamEstimation'
		%optFcn = @(p) checkParam(p,model,MCparams,dt,nTime,numel(a),len);
		optFcn = @(p) checkParam(p,model,MCparams,dt,nTime,numel(a),len,time);
        
        for k=1:length(model.p_estimate);
			lu(1,k)= model.p_lim_lower(model.p_estimate(k));
			lu(2,k)= model.p_lim_upper(model.p_estimate(k));
        end
        
        %The default parameters in SBML-PET-MODULE
		%popSize = length(model.p_estimate)*7+210;
		%maxGen = 2000;
		%parentNumber = length(model.p_estimate)+30;
        
		
		optParams.popSize = 50;
		optParams.maxGen = 50;
		optParams.minMax = 'max';
		optParams.lu = lu;
		
		%schd = findResource('scheduler', 'configuration', 'local');
		%if(matlabpool('size') == 0) 
			%matlabpool('local',schd.ClusterSize);
		%	matlabpool('local',8);
		%end
		
		[xb, stats] = optimizeParams(optMethod,optFcn,optParams);
		
		%matlabpool close;
        
        disp('Best Param');
        for k=1:length(xb)
            fprintf('%d\t',xb(k));
		end
        
		figure;plot([0 stats.bestStat.nEval(1:length(stats.bestStat.Best))],[0 stats.bestStat.Best]);
		
	case 'SensAnalysis'
		nPCombinations=20000;
        KS = paramSensitivityDiscret(model,MCparams,dt,nTime,nPCombinations,numel(a),len,time);
end

toc;
diary off;