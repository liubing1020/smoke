function [fval,gval] = checkParam(p,model,MCparams,dt,nTime,nDataFormulas,nDataSpecies,nTimepointdata)
	N = size(p,1);
	gval = -1*ones(N,1);
	
	%assemble the parameter tuples (random parameter sets for the one to 
	% estimate and nominal parameter cubes for the known ones)
	B = repmat(model.p_nominal,N,1);
	
	for i=1:length(model.p_estimate)
		B(:,(model.p_estimate(i)))= p(:,i);
    end
	   
	fval = paramStatisticalTest(B,model,MCparams,dt,nTime,nDataFormulas,nDataSpecies,nTimepointdata);	

end