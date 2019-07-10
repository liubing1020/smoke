function [xb, stats] = optimizeParams(method,optFcn,optParams)
	switch(method)
		case 'sres'
			parentNumber = ceil(optParams.popSize/7);
			pf = 0.45;
			varphi = 1;
			[xb, stats.stats, stats.Gm] = sres(optFcn, optParams.minMax, ...
				optParams.lu, optParams.popSize, optParams.maxGen, ...
				parentNumber, pf, varphi);
			stats.bestStat.nEval = (1:optParams.maxGen)*optParams.popSize;
			stats.bestStat.Best = cummax(stats.stats(:,1));
			
		case 'ga'
			if strcmp(optParams.minMax,'max')
				optFcnGA = @(x) -optFcn(x);
			else
				optFcnGA = @(x) optFcn(x);
			end
			options = gaoptimset('PopulationSize',optParams.popSize, ...
                'Vectorized','on','UseParallel','always', ...
                'Generations', optParams.maxGen , 'Display', 'off', ...
                'PopInitRange',optParams.lu,'OutputFcns',{@savegastate});
			[xb,stats.fval,stats.exitflag,stats.output,stats.population,stats.scores] = ga(optFcnGA,size(optParams.lu,2),[],[],[],[],...
				optParams.lu(1,:),optParams.lu(2,:),[],options);
			[stats.state] = savegastate([],[],[]);
			stats.bestStat.nEval = (1:optParams.maxGen)*optParams.popSize;
			stats.bestStat.Best = -(stats.state.Best);
		case 'abcsmc'
			[xb,stats] = abcsmc(optFcn, optParams.lu, 100);
		case 'randsamp'
			[xb,stats] = randsamp(optFcn, optParams.lu, 10000);
	end
end