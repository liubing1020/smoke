function cvode_compile(modelName,cvodePath)
	if nargin < 2
		cvodePath = './';
	end
	simPath = [cvodePath './'];

	compileString = ['mex '...
			'-g' ...									% optimize
			' -output ' modelName ...					% output name
			' -I' simPath ...							% cvode_sim location
			' -I' cvodePath '/cv_src/include ' ...		% CVODE includes
			cvodePath '/lib/CVODE.lib ' ...				% CVODE lib
			modelName '.cpp'];							% model mex file

	eval(compileString);
end