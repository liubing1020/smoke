%% Add paths
%clear variables
clear variables
addpath util/
addpath ode/
addpath analysis/
addpath modelcheck/
addpath ../models/segmentationclock/
addpath ../tools/sres/

%% Create the model
cd ode
compile('segclock','mc')
cd ..
model = seg_clock_pathway();
option = 'SensAnalysis';
optMethod = 'ga';
%optMethod = 'sres';
%optMethod = 'abcsmc';
%optMethod = 'randsamp';

tic;
%%The lower and upper limit for parameters

diary 'segclock_paramest_ga_09_p_2.txt';

diary on;

model.p_lim_lower = 0*model.p_nominal;
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
a=[];

% Statistical model checking parameters
MCparams.alpha=0.05;
MCparams.beta=0.05;
MCparams.gamma=0.05;
MCparams.delta=0.05;


%% Define all the rules here (300 minutes)
%MCparams.rules{1} = '([1==10]&F{3}([1==1]))&(F{}([1==1]&F([1==3]&F([1==1]&F([1==3]&F([1==1]&F([1==3])))))))';
%MCparams.rules{2} = '(([3==1])&(F{}([3==6]&F([3==1]&F([3==6]&F([3==1]&F([3==6]&F([3==1]))))))))&(G(![3>=8]))';
%MCparams.rules{3} = '(([4==1])&(F{}([4>=6]&F([4==1]&F([4>=6]&F([4==1]&F([4>=6]&F([4==1]))))))))&(G(![4>=9]))';
%MCparams.rules{4} = '(([11==1])&(F{}([11==8]&F([11==1]&F([11==8]&F([11==1]&F([11==8]&F([11==1]))))))))&(G(![11>=9]))';
%MCparams.rules{5} = '(([13==1])&(F{}([13==8]&F([13==1]&F([13==8]&F([13==1]&F([13==8]&F([13==1]))))))))&(G(![13>=9]))';
%MCparams.rules{6} = '(([15==1])&(F{}([15>=6]&F([15==1]&F([15>=6]&F([15==1]&F([15>=6]&F([15==1]))))))))&(G(![15>=8]))';

%Notch protein
%MCparams.rules{1} = '([1,0.45,0.55]&F{3}([1,0,0.05]))&(F{}([1,0,0.5]&F([1,0.10,0.15]&F([1,0,0.5]&F([1,0.10,0.15])))))';
%nuclear NicD
%MCparams.rules{2} = '(([3,0,0.012])&(F{}([3,0.07,0.08]&F([3,0,0.012]&F([3,0.07,0.08]&F([3,0,0.012]))))))&(G(![3,0.9,-1]))';
%Lunatic fringe mRNA
%MCparams.rules{3} = '(([4,0,0.4])&(F{}([4,2.2,-1]&F([4,0,0.4]&F([4,2.2,-1]&F([4,0,0.4]))))))&(G(![4,3,-1]))';
%Axin2 mRNA
%MCparams.rules{4} = '(([11,0,2])&(F{}([11,13,16]&F([11,0,2]&F([11,13,16]&F([11,0,2]))))))&(G(![11,18,-1]))';
%active ERK
%MCparams.rules{4} = '(([13,0,0.27] &F{3}([13,1.9,2.2]))&(F{}([13,1.9,2.2]&F([13,0,0.27]&F([13,1.9,2.2]&F([13,0,0.27]&(!(F([13,1.5,-1])))))))))&(G(![13,2.4,-1]))';
%Dusp6 mRNA
%MCparams.rules{5} ='(([15,0,1])&(F{}([15,5.5,-1]&F([15,0,1]&F([15,5.5,-1]&F([15,0,1]&(!(F([15,1,-1])))))))))&(G(![15,7,-1]))';

%FINAL SET USE THESE
%Notch protein
%MCparams.rules{1} = '(([1,0.45,0.55]&F{3}([1,0,0.05]))&(F{}([1,0,0.05]&F([1,0.10,0.15]&F([1,0,0.05]&F([1,0.10,0.15]))))))&(!(F{25}([1,0.45,0.55]&F([1,0,0.05]&F([1,0.10,0.15]&F([1,0,0.05]&F([1,0.10,0.15])))))))';
%nuclear NicD
%MCparams.rules{2} = '(([3,0,0.012])&(F{}([3,0.07,0.08]&F([3,0,0.012]&F([3,0.07,0.08]&F([3,0,0.012]))))))&((!(F{21}([3,0.07,0.08]&F([3,0,0.012]&F([3,0.07,0.08])))))&(G(![3,0.9,-1])))';
%Lunatic fringe mRNA
%MCparams.rules{3} = '(([4,0,0.4])&(F{}([4,2.2,-1]&F([4,0,0.4]&F([4,2.2,-1]&F([4,0,0.4]))))))&((!(F{21}(([4,2.2,-1]&F([4,0,0.4]&F([4,2.2,-1]))))))&(G(![4,3,-1])))';
%Axin2 mRNA
%MCparams.rules{4} = '(([11,0,2])&(F{}([11,13,16]&F([11,0,2]&F([11,13,16]&F([11,0,2]))))))&(G(![11,18,-1]))';
%active ERK
%MCparams.rules{4} = '(([13,0,0.27] &F{3}([13,1.9,2.2]))&(F{}([13,1.9,2.2]&F([13,0,0.27]&F([13,1.9,2.2]&F([13,0,0.27]))))))&((!(F{16}(([13,1.9,2.2]&F([13,0,0.27]&F([13,1.9,2.2]))))))&(G(![13,2.4,-1])))';
%Dusp6 mRNA
MCparams.rules{1} ='(([15,0,1])&(F{}([15,5.5,-1]&F([15,0,1]&F([15,5.5,-1]&F([15,0,1]))))))&((!(F{21}([15,5.5,-1]&F([15,0,1]&F([15,5.5,-1])))))&(G(![15,7,-1])))';

len=0;
time=0;
%Axin2 mRNA
%[a,len,time] = dataToIFormula('../../models/segmentationclock/Axin2mRNA.csv',0.1,2);

if numel(MCparams.rules)==0
	MCparams.rules = a;
else
	MCparams.rules = [MCparams.rules,a];
end

%%focus is on just observing oscillations with a certain amplitude (500 minutes)
%Notch protein
%MCparams.rules{1} = '([1,0.45,0.55]&F{3}([1,0,0.05]))&(F{}([1,0,0.05]&F([1,0.10,0.15]&F([1,0,0.05]&F([1,0.10,0.15]&F([1,0,0.05]&F([1,0.10,0.15]&F([1,0,0.05]&F([1,0.10,0.15]&F([1,0,0.05]&F([1,0.10,0.15]&F[1,0,0.05])))))))))))';
%nuclear NicD
%MCparams.rules{2} = '(([3,0,0.012])&(F{}([3,0.07,0.08]&F([3,0,0.012]&F([3,0.07,0.08]&F([3,0,0.012]&F([3,0.07,0.08]&F([3,0,0.012]&F([3,0.07,0.08]&F([3,0,0.012]&F([3,0.07,0.08]&F([3,0,0.012]))))))))))))&(G(![3,0.9,-1]))';
%Lunatic fringe mRNA
%MCparams.rules{3} = '(([4,0,0.4])&(F{}([4,2.2,-1]&F([4,0,0.4]&F([4,2.2,-1]&F([4,0,0.4]&F([4,2.2,-1]&F([4,0,0.4]&F([4,2.2,-1]&F([4,0,0.4]&F([4,2.2,-1]&F([4,0,0.4]))))))))))))&(G(![4,3,-1]))';
%Axin2 mRNA
%MCparams.rules{4} = '(([11,0,2])&(F{}([11,13,16]&F([11,0,2]&F([11,13,16]&F([11,0,2]&F([11,13,16]&F([11,0,2]&F([11,13,16]&F([11,0,2]&F([11,13,16]&F([11,0,2]))))))))))))&(G(![11,18,-1]))';
%active ERK
%MCparams.rules{5} = '(([13,0,0.27] &F{3}([13,1.9,2.2]))&(F{}([13,1.9,2.2]&F([13,0,0.27]&F([13,1.9,2.2]&F([13,0,0.27]&F([13,1.9,2.2]&F([13,0,0.27]&F([13,1.9,2.2]&F([13,0,0.27]&F([13,1.9,2.2]&F([13,0,0.27]))))))))))))&(G(![13,2.4,-1]))';
%Dusp6 mRNA
%MCparams.rules{6} = '(([15,0,1])&(F{}([15,5.5,-1]&F([15,0,1]&F([15,5.5,-1]&F([15,0,1]&F([15,5.5,-1]&F([15,0,1]&F([15,5.5,-1]&F([15,0,1]&F([15,5.5,-1]&F([15,0,1]))))))))))))&(G(![15,7,-1]))';


%Notch protein
%MCparams.rules{1} = '([1==10]&F{3}([1==1]))&(F{}([1==1]&F([1==3]&F([1==1]&F([1==3])))))';
%nuclear NicD
%MCparams.rules{2} = '(([3==1])&(F{}([3>=6]&F([3==1]&F([3>=6]&F([3==1]))))))&(G(![3>=8]))';
%Lunatic fringe mRNA
%MCparams.rules{3} = '(([4==1])&(F{}([4>=6]&F([4==1]&F([4>=6]&F([4==1]))))))&(G(![4>=9]))';
%Axin2 mRNA
%MCparams.rules{4} = '(([11==1])&(F{}([11==8]&F([11==1]&F([11==8]&F([11==1]))))))&(G(![11>=9]))';
%active ERK
%MCparams.rules{5} = '(([13==1] &F{3}([13==8]))&(F{}([13==8]&F([13==1]&F([13==8]&F([13==1]))))))&(G(![13>=9]))';
%Dusp6 mRNA
%MCparams.rules{6} = '(([15==1])&(F{}([15>=6]&F([15==1]&F([15>=6]&F([15==1]))))))&(G(![15>=8]))';


%% Define probabilities for each rule
%MCparams.probs = [0.8 0.8 0.8 0.8 0.8 0.8];
%MCparams.probGreater = [1 1 1 1 1 1];

MCparams.probs = 0.8*ones(size(MCparams.rules));
MCparams.probGreater = ones(size(MCparams.rules));

% Time between model checking steps (timse discretization)
dt = 5;
% Number of discrete time steps
nTime = 40;

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
        
		
		optParams.popSize = 200;
		optParams.maxGen = 300;
		optParams.minMax = 'max';
		optParams.lu = lu;
        
		
		[xb, stats] = optimizeParams(optMethod,optFcn,optParams);
        
       
        disp('Best Param');
        for k=1:length(xb)
            fprintf('%d\t',xb(k));
		end
        
		figure;plot([0 stats.bestStat.nEval(1:length(stats.bestStat.Best))],[0 stats.bestStat.Best]);
	case 'SensAnalysis'
		nPCombinations=500000;
		KS = paramSensitivityDiscret(model,MCparams,dt,nTime,nPCombinations,numel(a),len,time);
end

toc;

diary off;