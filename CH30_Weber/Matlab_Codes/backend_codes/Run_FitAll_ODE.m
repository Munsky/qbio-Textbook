% Run_FitAll_ODE: Optimization with multiple sweeps to fit parameters for all models
% Set initial difference as infinity for all four models
BestResults_Difference(1,1:4)=inf;

% This section will iterate through the optimization code for each
% model for a total of 10 iterations, using best fit parameters as a
% starting point 
[NewResults_Parameters,NewResults_Difference] = ...
    TestAllFun(Input,Output_Times_data,x0,Data_Set,All_Guesses);
All_Guesses = NewResults_Parameters;
for sweep=2:30
    disp(' ');
    disp(['Running sweep ',num2str(sweep),' out of 30.'])
    OldResults_Difference=NewResults_Difference;
    [NewResults_Parameters,NewResults_Difference] = ...
        TestAllFun(Input,Output_Times_data,x0,Data_Set,All_Guesses);  % Try again      
    for N=4:-1:1
            % If the new results are better than any previous results then
            % store the new results
        if  NewResults_Difference(N)<BestResults_Difference(N)
            BestResults_Difference(N) = NewResults_Difference(N);
            BestResults_Parameters(N,:) = NewResults_Parameters(N,:);
            All_Guesses(N,:) = NewResults_Parameters(N,:);
        elseif NewResults_Difference(N)==OldResults_Difference(N) %Local minimum
            All_Guesses(N,:) = 2*BestResults_Parameters(N,:).*rand(1,8); % Randomize initial guess
        elseif NewResults_Difference(N)<OldResults_Difference(N)
            All_Guesses(N,:) = NewResults_Parameters(N,:);
            % Otherwise, if the best result has not changed, randomize the
            % initial parameter guesses to get out of the local minimum
        elseif NewResults_Difference(N)>OldResults_Difference(N)
            All_Guesses(N,:) = abs(sqrt(2).*randn(1,8)); % Randomize initial guess
        end
    end
end

% One last time to put them all together
disp(' ');
disp('Now to collect the final results')
[BestResults_Parameters,BestResults_Difference] = ...
    TestAllFun(Input,Output_Times_data,x0,Data_Set,BestResults_Parameters);
