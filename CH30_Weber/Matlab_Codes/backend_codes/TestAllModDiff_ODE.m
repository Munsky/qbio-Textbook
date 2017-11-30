% This function will test all 4 different models (each with a different
% kinase-dependent rate) and perform the optimization to determine best
% parameter values for each model
function[Model_Results_Parameters, Model_Results_Differences] = TestAllModDiff_ODE(Input,Output_Times,x0,Data_Set,All_Guesses)
for i=1:size(All_Guesses,1)
    % Save model differences and parameters in a matrix
    [Model_Results_Differences(i),Model_Results_Parameters(i,:)] = MinDifference_ODE(Input,i,Output_Times,x0,Data_Set,All_Guesses(i,:))
    % Plot and label the fitted model and actual data for each model
    subplot(2,2,i)
       plot(Output_Times,MeanTrajectory_ODE( Model_Results_Parameters(i,:),...
        i, Input,Output_Times,x0),'b-',Output_Times,mean(Data_Set),'r.')
    title(sprintf('Fitted Model Number : %i',i)) % Adding title
    xlabel('Time', 'FontWeight', 'bold') % Adding label to x axis
    ylabel('mRNA count','FontWeight', 'bold') % Adding label to y axis
    legend('Fitted Model','Actual Data', 'Location', 'SouthEast') % Adding legend
    drawnow % Draws plot immediately
    % Save the parameter results for each model and update each sweep
    eval(sprintf('save generated_results/BestParams_%i_ODE.mat Model_Results_Parameters',size(All_Guesses,1)))
end 
end

