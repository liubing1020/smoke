function d = discretize(x,levels)
	d = zeros(size(x));
	for i=1:length(x)
		d(i) = min(find(x(i) < [levels{i},Inf],1)-1,length(levels{i})-1);
	end
end