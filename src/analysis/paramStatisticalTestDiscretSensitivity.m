%performs the statistical test and computes the objective value for the sensitivity analysis 

function objVal = paramStatisticalTestDiscretSensitivity(p,model,MCparams,dt,nTime,nPCombinations,nDataFormulas,nDataSpecies,nTimepointdata)
	nRules = length(MCparams.rules);
	N = size(p,1);
	objVal = zeros(N,1);
    
    for k=1:length(model.p_estimate),
		Par_limits{k}=linspace(model.p_lim_lower(model.p_estimate(k)),model.p_lim_upper(model.p_estimate(k)),nPCombinations+1);
	end
 
batchSize = 10000;
 for a= 1:(N/batchSize)
     
    %schd = findResource('scheduler', 'configuration', 'local');
	%if(matlabpool('size') == 0) 
	%	matlabpool('local',schd.ClusterSize);
    %end
    
    for k=((a-1)*batchSize+1):(a*batchSize)
		RandStream.setGlobalStream(RandStream('mt19937ar','seed',k*sum(clock)));
		posSamples = zeros(1,nRules);
		numSamples = zeros(1,nRules);
		stopFlag = zeros(1,nRules);
		finalResults = -1*ones(1,nRules);
		
		beta_inside = MCparams.beta*ones(1,numel(MCparams.rules));
		counter=0;	
		counter_species = zeros(1,nDataSpecies);
		flag=0;
        
		while(any(stopFlag == 0))
            
            pRand=model.p_nominal;
		    
			for i=1:length(model.p_estimate),
				pRand(model.p_estimate(i)) = Par_limits{i}(p(k,model.p_estimate(i))) + ... 
					(Par_limits{i}(p(k,model.p_estimate(i))+1) - ...
					Par_limits{i}(p(k,model.p_estimate(i))))*rand;
			end
            
			%initial values picked from the initial discretized interval 
            xInit = samplingScheme(model.x_init,'uniform',0.05);
            
			tf = runMC_cvode(pRand,xInit,nTime,MCparams.rules,dt);

            
			% Look for how many postive and total samples are needed for each rule
			posSamples(~stopFlag & tf) = posSamples(~stopFlag & tf)+1;
			numSamples(~stopFlag) = numSamples(~stopFlag)+1;

			% Run the statistical hypothesis testing procedure
			for i=1:nRules
				if(stopFlag(i)==0)
					switch(MCparams.method)
						case 'YounesA'
							% Younes A Algorithm, with no undecided results (takes lesser samples)
							[stopFlag(i),finalResults(i)] = Younes_A( MCparams.alpha, beta_inside(i), MCparams.delta, ... 
								MCparams.probs(i), numSamples(i), posSamples(i),MCparams.probGreater(i));
						case 'YounesB'
							% Younes B Algorithm, with undecided results
							[stopFlag(i),finalResults(i)] = Younes_B( MCparams.alpha, beta_inside(i), MCparams.gamma, ...
								MCparams.delta, MCparams.probs(i), numSamples(i), posSamples(i), MCparams.probGreater(i) );
					end
				end
			end
			
			if(nDataFormulas>0)
			for i=(numel(MCparams.rules)-nDataFormulas)+1:(nDataFormulas/nDataSpecies):numel(MCparams.rules)
				if(counter_species(ceil((i-(numel(MCparams.rules)-nDataFormulas))/(nDataFormulas/nDataSpecies)))~=1)	
					if(ismember(1,finalResults(i:i+(nDataFormulas/nDataSpecies)-1))&& ~ismember(-1,finalResults(i:i+(nDataFormulas/nDataSpecies)-1)))
						for kk=i:i+(nDataFormulas/nDataSpecies)-1
							if(finalResults(kk)==1)
								beta_inside(kk)=beta_inside(kk)/(nDataFormulas/nDataSpecies);
								stopFlag(kk)=0;
								finalResults(kk)=-1;
							end
						end	
						counter_species(ceil((i-(numel(MCparams.rules)-nDataFormulas))/(nDataFormulas/nDataSpecies)))=1;
					end
				end
			end
			end
			
			if(all(finalResults==1)&&counter<1)
				stopFlag = zeros(1,nRules);
				counter = counter+1;
				beta_inside=beta_inside/(numel(MCparams.rules)-nDataFormulas+nDataSpecies);
				flag=1;
			end
			
		end
		
		fprintf('%d \t',max(numSamples));
		
		objVal(k)=sum(finalResults(1:(numel(MCparams.rules)-nDataFormulas)));
		
		if(nDataFormulas>0)
		for i=(numel(MCparams.rules)-nDataFormulas)+1:(nDataFormulas/nDataSpecies):numel(MCparams.rules)
			objVal(k)=objVal(k)+ sum(finalResults(i:i+(nDataFormulas/nDataSpecies)-1))/numel(finalResults(i:i+(nDataFormulas/nDataSpecies)-1));				
		end
		end
		
		objVal(k)=objVal(k)+mean(posSamples./numSamples);
		
		if flag==1 && all(finalResults==1)
		objVal(k)= objVal(k)+1;
		end
		
		%probMeasure = mean(posSamples./numSamples);
		%objVal(k) = (sum(finalResults == 1)) + probMeasure;
        fprintf('population : %d objective value: %f\n',k,objVal(k));
        
    end
    %matlabpool close;
end
    
end