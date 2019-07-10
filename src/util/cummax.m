function cmv  = cummax(v)
	cmv = zeros(1,length(v));
	for i=1:length(v)
		cmv(i) = max([0; v(1:i-1)]);
	end
end