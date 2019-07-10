%% Add paths
clear variables
addpath util/
addpath ode/
addpath modelcheck/
addpath analysis/
addpath ../models/thrombin/
addpath ../tools/sres/

%% Create the model
cd ode
compile('thrombin','mc');
cd ..
model = thrombin_pathway();
option = 'ParamEstimation';
optMethod = 'sres';
%optMethod = 'sres';
%optMethod = 'abcsmc';
%optMethod = 'randsamp';
a=[]; time=0; len=0;
tic;

diary 'thrombin_sres_08-0.5-x.txt';

diary on;
%%The lower and upper limit for parameters

model.p_lim_lower = 0.5*model.p_nominal;
model.p_lim_upper = 2.0*model.p_nominal;


for i=1:length(model.p_nominal)
    if(~any(find(model.p_estimate==i)))
    model.p_lim_lower(i)=model.p_nominal(i);
    model.p_lim_upper(i)=model.p_nominal(i);
    end
end

% Statistical model checking method
MCparams.method = 'YounesA';
MCparams.rules=[];

% Statistical model checking parameters
MCparams.alpha=0.05;
MCparams.beta=0.05;
MCparams.gamma=0.05;
MCparams.delta=0.05;


MCparams.rules{1} = '((F{2}([6,0.004,0.0046]))&(F{5}([6,0.004,0.0046]&F([6,0,0.00006]))))& F{}(([6,0.004,0.0046])&F([6,0,0.00008]&G([6,0,0.00008])))';
MCparams.rules{2} = '(((F{2}([16,0.0007,0.0009]))&(F{5}([16,0.0007,0.0009]&F([16,0,0.00009])))))&F{}(([16,0.0007,0.0009])&F([16,0,0.00009]&G([16,0,0.00009])))';

[a,len,time] = dataToIFormula('../models/thrombin/training.csv',0.1,4);

if numel(MCparams.rules)==0
	MCparams.rules = a;
else
	MCparams.rules = [MCparams.rules,a];
end

%[MCparams.rules,len] = dataToIFormulaPerTime('../../models/thrombin/training.csv',0.1);

MCparams.probs = 0.8*ones(size(MCparams.rules));
MCparams.probGreater = ones(size(MCparams.rules));

% Time between model checking steps (timse discretization)
dt = 50;
% Number of discrete time steps
nTime = 20;

switch(option) 
	case 'ParamEstimation'
				
		optFcn = @(p) checkParam(p,model,MCparams,dt,nTime,numel(a),len,time);
        
        for k=1:length(model.p_estimate);
		lu(1,k)= model.p_lim_lower(model.p_estimate(k));
        lu(2,k)= model.p_lim_upper(model.p_estimate(k));
        end
        
        %The default parameters in SBML-PET-MODULE
		%popSize = length(model.p_estimate)*7+210;
		%maxGen = 2000;
		%parentNumber = length(model.p_estimate)+30;
        

		optParams.popSize = 100;
		optParams.maxGen = 1000;
		optParams.minMax = 'max';
		optParams.lu = lu;
		
		
		[xb, stats] = optimizeParams(optMethod,optFcn,optParams);
		
        
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