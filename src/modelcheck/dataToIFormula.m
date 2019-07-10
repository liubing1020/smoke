function [Expr,len,time] = dataToIFormula(dataFile,noise,PointsToClub)

%[numData,textData,rawData]=xlsread(dataFile);

data = importdata(dataFile);
speciesNames = data.colheaders(2:end);
timePoints = data.data(:,1);
measurements = data.data(:,2:end);


tt=1;
for i=1:(length(speciesNames))
	openbrac=0;
    j=1;
	while j<(length(timePoints)+1)
		if((j+PointsToClub)<=(length(timePoints)+1))
			for k=1:PointsToClub
				   if k==1
						   openbrac=1;
						   if(openbrac~=0)
							   A='(';
						   end
						   if(PointsToClub>1)
						   A=strcat(A,'I{',num2str(timePoints(j)+1),'}[',(speciesNames{i}(2:end)),',',num2str((1-noise)*measurements(j,i)),',',num2str((1+noise)*measurements(j,i)),']');
						   elseif(PointsToClub==1)           
						   A=strcat(A,'I{',num2str(timePoints(j)+1),'}[',(speciesNames{i}(2:end)),',',num2str((1-noise)*measurements(j,i)),',',num2str((1+noise)*measurements(j,i)),'])');
						   end
				   elseif k>1
					   A=strcat(A,'&');
					   if(openbrac==0)
							A=strcat('(',A);   
						   openbrac=~openbrac;
					   end

					   A=strcat(A,'I{',num2str(timePoints(j)+1),'}[',(speciesNames{i}(2:end)),',',num2str((1-noise)*measurements(j,i)),',',num2str((1+noise)*measurements(j,i)),']');

					   if(openbrac==1)
							A=strcat(A,')');
							openbrac=~openbrac;
					   end
				   end
				   j=j+1;
			end
		 Expr{tt}= A;
		 tt=tt+1;
		elseif((j+PointsToClub)>=(length(timePoints)+1))
		  for k=1:((length(timePoints)+1)-j)
				   if k==1
						   openbrac=1;
						   if(openbrac~=0)
							   A='(';
						   end
						   
						   if(((length(timePoints)+1)-j)>1)
						   A=strcat(A,'I{',num2str(timePoints(j)+1),'}[',(speciesNames{i}(2:end)),',',num2str((1-noise)*measurements(j,i)),',',num2str((1+noise)*measurements(j,i)),']');
				           elseif(((length(timePoints)+1)-j)==1)           
						   A=strcat(A,'I{',num2str(timePoints(j)+1),'}[',(speciesNames{i}(2:end)),',',num2str((1-noise)*measurements(j,i)),',',num2str((1+noise)*measurements(j,i)),'])');
						   end
				   elseif k>1
					   A=strcat(A,'&');
					   if(openbrac==0)
							A=strcat('(',A);   
						   openbrac=~openbrac;
					   end

					   A=strcat(A,'I{',num2str(timePoints(j)+1),'}[',(speciesNames{i}(2:end)),',',num2str((1-noise)*measurements(j,i)),',',num2str((1+noise)*measurements(j,i)),']');

					   if(openbrac==1)
							A=strcat(A,')');
							openbrac=~openbrac;
					   end
				   end
				  j=j+1;
			end
		 Expr{tt}= A;
		 tt=tt+1;
		end
	 A=[];
	end

end
time = length(timePoints);
len = length(speciesNames);
end
 



