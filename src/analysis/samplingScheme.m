function xInitOut = samplingScheme(xInit,option,range)
    switch(option)
        case 'uniform'
			xInitOut = max(0,xInit + 2*range*(rand(size(xInit))-0.5).*xInit);
        %make changes
        case 'guassian'
            xInitOut = (1-range).*xInit + (2*range*rand()).*xInit;
        case 'lognormal'
            xInitOut = (1-range).*xInit + (2*range*rand()).*xInit;
    end
end