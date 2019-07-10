function [ stopFlag,finalResults ] = Younes_A( alpha, beta, delta, probs, numSamples, posSamples, probGreater )
%Younes_A
	
	acceptance_number = floor((log(beta/(1-alpha)) + ...
						numSamples*(log((1-(probs+delta))/(1-(probs-delta))))) / ...
						((log((probs-delta)/(probs+delta))- log((1-(probs-delta))/(1-(probs+delta))))));
    rejection_number = ceil((log((1-beta)/alpha) + ...
						numSamples*(log((1-(probs+delta))/(1-(probs-delta)))))/ ...
						((log((probs-delta)/(probs+delta))- log((1-(probs-delta))/(1-(probs+delta))))));

	if (posSamples>acceptance_number)
		if(probGreater==1)
			finalResults=1;
		elseif(probGreater==0)
			finalResults=0;
        end
		stopFlag=1;
        
	elseif (posSamples<rejection_number)
		if(probGreater==1)
			finalResults=0;
		elseif(probGreater==0)
			finalResults=1;
		end
		stopFlag = 1;
        
    %added an additional condition based on p-value when the number of samples exceeds 100
    elseif (numSamples>=100000)
        a=binocdf(posSamples,numSamples,probs);
        b = 1-a;
        if(probGreater==1)
            if(a>b)
                finalResults=1;
            else
                finalResults=0;
            end
            
		elseif(probGreater==0)
            if(a>b)
                finalResults=0;
            else
                finalResults=1;
            end
        end  
        stopFlag=1;
        
	else
		stopFlag=0;
		finalResults=-1;
	end
	
end

