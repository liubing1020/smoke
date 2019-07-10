function compile(modelName,type)
	if nargin < 1
		modelDefine = 'MODEL_THROMBIN';
	else
		switch(modelName)
			case 'segclock'
				modelDefine = 'MODEL_SEGCLOCK';
			case 'thrombin'
				modelDefine = 'MODEL_THROMBIN';
			case 'egfngf'
				modelDefine = 'MODEL_EGFNGF';
			case 'repress'
				modelDefine = 'MODEL_REPRESS';
			case 'tlr'
				modelDefine = 'MODEL_TLR';
		end
	end
	
	if nargin < 2
		target = 'runMC_cvode';
	else
		switch(type)
			case 'mc'
				target = 'runMC_cvode';
			case 'graph'
				target = 'Gen_Graph';
		end
	end

	if ispc
		libpath = 'cvode/lib/CVODE.lib';
	elseif ismac
		libpath = 'cvode/lib/CVODE_mac.a';
	else
		libpath = 'cvode/lib/CVODE_unix.a';
	end
	
	fprintf('Compiling %s for %s pathway...\n',target,modelDefine);


	a = ['mex COMPFLAGS="$COMPFLAGS" -g ' ...
			'-output ' target ' '...
			'-D',modelDefine,' ' ...
			'-Icvode ' ...
			'-Icvode/cv_src/include ' ...
			'-I../modelcheck ' ...
			'-I../util ' ...
			libpath ' '  ...
			target '.cpp' ];

	eval(a);
end