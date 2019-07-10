function [Expr,len] = dataToIFormulaPerTime(dataFile,noise)

%read data
%[numData,textData,rawData]=xlsread(dataFile);
data = importdata(dataFile);

speciesNames = data.colheaders(2:end);
timePoints = data.data(:,1);
measurements = data.data(:,2:end);
%make a seperate I formula for each experimental data point
index=1;
for i=1:length(speciesNames)
    for j=1:length(timePoints)
        A=strcat('I{',num2str(timePoints(j)+1),'}[',(speciesNames{i}(2:end)),',',...
			num2str((1-noise)*measurements(j,i)),',',num2str((1+noise)*measurements(j,i)),']');
        Expr{index}= A;
        index=index+1;
    end

end

len = length(speciesNames);
end




