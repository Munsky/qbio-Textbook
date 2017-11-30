%This function will compute normalized-sum-of-squares difference of an experimental data set from the model.
function [difference] = Difference_ODE(Parameters,Input,Nm,Output_Times,x0,Data_Set)

[Y]=MeanTrajectory_ODE(Parameters,Nm,Input,Output_Times,x0);  %Get mRNA mean trajectory from ODE solver
S=std(Data_Set)/sqrt(size(Data_Set,1));  % Get standard error of data at each time point
SSE=((Y-mean(Data_Set))./S).^2;   % Calculate the sum of err sqrd with the mean trajectory of the data set
difference=sum(SSE(2:end)); % Remove first point (nan value) and sum the squared error

%% Add effect of log-normal prior on parameters
LogPar = log10(abs(Parameters));
sigma_log10 = 2;
mu_log10 = zeros(size(Parameters));
ErrLogPar = (LogPar-mu_log10).^2/(2*sigma_log10^2);  % The prior probability of current parameter set -- up to a constant.
difference = difference + sum(ErrLogPar);            % log(P(D|lam)*P(lam)) = log(P(D|lam)) + log(P(lam)) 

end