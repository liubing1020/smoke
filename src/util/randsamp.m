function [xb,stats] = randsamp(optFcn, lu, N)
	np = size(lu,2);
	ps = rand(N,np).*(repmat(lu(2,:),N,1)-repmat(lu(1,:),N,1))+repmat(lu(1,:),N,1);
	
	fval = optFcn(ps);
	[~,maxIdx] = max(fval);
	xb = ps(maxIdx,:);

	stats.bestStat.nEval = 1:N;
	stats.bestStat.Best = cummax(fval);
end