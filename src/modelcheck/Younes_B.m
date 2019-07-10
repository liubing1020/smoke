function [ stopFlag,finalResults ] = Younes_B( alpha, beta, gamma, delta, probs, numSamples, posSamples, probGreater )
%Younes with undecided- algorithm
	acceptance_number_1 = floor((log(gamma/(1-alpha)) + ...
					numSamples*(log((1-(probs))/(1-(probs-delta))))) / ...
					((log((probs-delta)/(probs))- log((1-(probs-delta))/(1-(probs))))));
	rejection_number_1 = ceil((log((1-gamma)/alpha) + ... 
					numSamples*(log((1-(probs))/(1-(probs-delta))))) / ...
					((log((probs-delta)/(probs))- log((1-(probs-delta))/(1-(probs))))));
	acceptance_number_2 = floor((log(beta/(1-gamma)) + ... 
					numSamples*(log((1-(probs+delta))/(1-(probs))))) / ...
					((log((probs)/(probs+delta))- log((1-(probs))/(1-(probs+delta))))));
	rejection_number_2 = ceil((log((1-beta)/gamma) + ... 
					numSamples*(log((1-(probs+delta))/(1-(probs))))) / ... 
					((log((probs)/(probs+delta))- log((1-(probs))/(1-(probs+delta))))));
 
	if (posSamples>acceptance_number_1 && posSamples>acceptance_number_2)
		if(probGreater==1)
			finalResults=1;
		else
			finalResults=0;
		end
		stopFlag=1;
	elseif (posSamples<rejection_number_2 && posSamples<rejection_number_1)
		if(probGreater==1)
			finalResults=0;
		else
			finalResults=1;
		end
		stopFlag = 1;
	elseif (((posSamples>acceptance_number_1) && (posSamples<rejection_number_2)) || ...
			((posSamples>acceptance_number_2) && (posSamples<rejection_number_1)) )
		stopFlag=1;
		finalResults=2;
        
     elseif (numSamples>=100)
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

