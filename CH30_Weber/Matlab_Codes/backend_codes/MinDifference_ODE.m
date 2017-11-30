% This function determines the best parameter set by minimizing the
% normalized sum of squared residuals 
function [min_difference,BestParameters] = MinDifference_ODE(Input,Nm,Output_Times,x0,Data_Set,Parameter_Guess)
% Define 'difference' function to be fed into optimization 
F=@(x)Difference_ODE(10.^x,Input,Nm,Output_Times,x0,Data_Set);
% Redefine current parameters in logspace
CurParameters = log10(Parameter_Guess);
% Options for fminsearch
options = optimset('display','final','MaxIter',500);
% Perform optimization using fminsearch before computing chi-squared function
[Parameters,difference]=fminsearch(F,CurParameters,options);    
[min_difference] = Difference_ODE(10.^Parameters,Input,Nm,Output_Times,x0,Data_Set);
CurParameters = Parameters;
BestParameters = 10.^CurParameters; % Puts parameters back in linear space
end



