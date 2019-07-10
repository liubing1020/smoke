%sensitivity analysis function, does the LHS sampling

function KS = paramSensitivityDiscret(model,MCparams,dt,nTime,nPCombinations,nDataFormulas,nDataSpecies,nTimepointdata)
	% Sensitivity analysis - discretize the parameter space into a hypercube
	
	% Pick M combinations of random parameters from the combination of
	% parameters whose sensitivity you want to assess, for the remaining
	% prameters fix the nominal discretization
	%ParamLevInit = Levels(discretize(model.p_nominal,pLevels));
    
    for kk=1:length(model.p_estimate)
        for k=1:nPCombinations
            A(k,kk)=k;
        end
	end
	
	%randomize the matrix
	A=randmat(A,1);
	A=shake(A,1);
	
	B = repmat(model.p_nominal,nPCombinations,1);
	
	for i=1:length(model.p_estimate)
		B(:,(model.p_estimate(i)))= A(:,i);
	end
	
	% For each of the chosen set of parameter values compute the corresponding
	% objectve value
	

	
	objValue = paramStatisticalTestDiscretSensitivity(B,model,MCparams,dt,nTime,nPCombinations,nDataFormulas,nDataSpecies,nTimepointdata);
	

	
	% Subjective threshold to decide if the set of parameter are good or bad
	pThreshold = mean(objValue);
	
	% For each parameter value compute the cumulative frequency for the
	% acceptable and unacceptable case.
	pAcceptableTuples = A(objValue>=pThreshold,:);
	pUnAcceptableTuples = A(objValue<pThreshold,:); 
  
	%% Plot the CDF and compute the KS statistic
	for j=1:length(model.p_estimate)
		%figure; hold all;
		%cdfplot(pAcceptableTuples(:,(model.p_estimate(j))));
		%cdfplot(pUnAcceptableTuples(:,(model.p_estimate(j))));
		[~, ~, KS(j)]  = kstest2(pAcceptableTuples(:,j),pUnAcceptableTuples(:,j));
	end 
	disp(KS);
    bar(KS);
end