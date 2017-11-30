%% Script to Run All Functions and Generate All Plots for Chapter 30: Gene Regulation Model Tutorial

% Lisa Weber: PhD Graduate Student, Colorado State University, Fort Collins, CO
% This code runs all functions and produces all plots in the gene regulation
% tutorial chapter.  These codes are not necessarily optimized for
% efficiency, but do adequately perform the tasks outlined.

addpath(genpath('Matlab_Codes'));

%% Functions (and m-files) to be called throughout code
% Functions for (Ordinary Differential Equation (ODE) Analysis:
Type = 'ODE';  
xdotfun = str2func(['Model_',Type]);  %Model for ODEs
odefun = str2func(['MeanTrajectory_',Type]); %Solves ODEs
ODEDifffun = str2func(['Difference_',Type]); %Computes chi-squared function (ODEs)
MinODEDifffun = str2func(['MinDifference_',Type]); %Minimizes chi-squared function (ODEs)
TestAllFun = str2func(['TestAllModDiff_',Type]);  %Performs optimization to fit all four models
RunModFitAll = ['Run_FitAll_',Type,'.m']; %Run file for fitting all four models
RunMetHast = ['Run_MH_',Type,'.m']; %Run file for ODE Metropolis-Hastings Search

% Functions for Stochastic Simulation Algorithm (SSA) Analysis:
Type = 'SSA';  
SSA_Mod = str2func(['Model_',Type]); %Model for SSA (S, W, w0, W1)
SSAfun = str2func(['SSA']); %Performs SSA for single trajectory
SSAHistfun = str2func(['Histogram_',Type]); %Performs SSA for many trajectories

% Functions for Finite-State Projection (FSP) Analysis:
Type = 'FSP'; 
TestAllFunFSP = str2func(['TestAllModDiff_',Type]);  %Performs optimization to fit all four models
RunMetHastFSP = ['Run_MH_',Type,'.m']; %Run file for FSP Metropolis-Hastings Search

% Additional functions used to perform the tasks outlined in this chapter:
% Analysis / Algorithm Codes: SSA.m, MetHaste.m, get_FSP_OBJ.m, FSP_Generate_Data.m
% Plotting Codes from File Exchange on http://www.mathworks.com: 
% error_ellipse_fill.m, shadedplot.m, suplabel.m

%% ALL TASKS: Model Input, Data Set, Model Number, and Initial Conditions

% Experiments: Input signal functions to evaluate kinase-dependent rate
Input1 = @(t)(1-cos(2*pi/30*t))*(t>5)*(t<70);  % Input signal 1: Sinusoidal
Input2 = @(t)(2)*(t>5)*(t<70);                 % Input signal 2: Step
Input3 = @(t)((2/65)*t-(2/13))*(t>5)*(t<70);   % Input signal 3: Ramp

Nm = 2; % True Model Number: defines which transition rate is input-dependent

% Actual parameter values which were used to generate data
% Parameters = [k12 k23 k21 k32 kr2 kr3 b g]; b=beta and g = gamma
Parameters_True = [0.2 -0.076 0.05 0.05 16 5 0.1 1]; 

x0 = [1;0;0;0]; % Initial condition x0 = [S1 S2 S3 R]'

% Vectors with output times specified 
Output_Times = linspace(0,100,100); 
Output_Times_data = linspace(0,100,10); % data has 10 time points

% Import simulated experimental data file and data histogram file
Data_Set = csvread('sim_data_bimodal_fsp.csv'); % Simulated Data
save('sim_data_bimodal.mat','Data_Set'); % Save file as 'Data_Set'
load('sim_data_bimodal.mat'); % Load Simulated Data
Data_Set_Hists = csvread('sim_data_bimodal_fsp_Data_Set_Hists.csv'); % Simulated Data Histograms
save('sim_data_bimodal_Hists.mat','Data_Set_Hists'); % Save file as 'Data_Set_Hists'
load('sim_data_bimodal_Hists.mat'); % Load Simulated Data Histograms

%% Task 2: ODE Chi-Squared Function (Difference) Calculation
% Computes the chi-squared function using the true parameter set that was
% used to generate the simulated experimental data. 
% These parameters are often not available, but are shown here to get an
% idea of the magnitude of chi-squared function values.

% Calculate chi-squared function based on actual parameter values
[difference] = ODEDifffun(abs(Parameters_True),Input1,Nm,Output_Times_data,x0,Data_Set);

%% Task 2: Test all 4 models - with different kinase-dependent rates - for All Three Inputs
% Runs fitting for all four models by performing optimization to minimize
% normalized chi-squared function. Runs sweep 30 times. Each sweep
% randomizes the initial parameter guess.

for input_num = 1:3  % Loop through all three input signals
    figure()
    % Set initial parameter guesses as 'ones' for all four models
    for N=4:-1:1
        Par = ones(1,8);
        All_Guesses(N,1:8) = Par;
    end
    
    % Specify which input is being used
    eval(['Input = Input',num2str(input_num)]);
    
    % Perform the optimization to fit each model and determine best parameter values
    run(RunModFitAll) 
    
    % Load results from fitting and rename file for each specific input
    load BestParams_4_ODE
    eval(['Model_Results_Parameters_I',num2str(input_num),' = Model_Results_Parameters;'])
end

%% Task 3: Run the Metropolis-Hastings Search for All Three Inputs with ODE Objective Function
% Runs the MH search of the parameter space using the chi-squared function 
% from the ODE analysis as the objective function

% N_Chain = length of the MH chain to be run for current parameter analysis
% N_Burn = # MH steps to run before starting to record results
% N_Thin = # MH steps to skip between each result output
% Npart_ODE = # MH runs with chain length of N_Chain
% Total chain length of MH run is N_Chain*N_Thin*Npart_ODE
N_Chain = 500; N_Burn = 0; N_Thin = 10; Npart_ODE=200;

% Output times for computing chi-squared function (same time points as data)
Output_Times = Output_Times_data; 

for Nm = 1:4  % Loop through all model numbers
    for input_num = 1:3 % Loop through all input signals
        
        % Load best fit parameters from minimizing (chi-squared function) residuals
        % and use as starting point for MH search
        Init_Param = eval(['Model_Results_Parameters_I',num2str(input_num),'(Nm,:);']);
        eval(['Input = Input',num2str(input_num),';']);

        % Run the MH search for the model
        run(RunMetHast)

        % Load resulting data file and save function evaluations and parameters
        MH_Results = csvread(STR_ODE_MH);
        eval(['Rslts_ODE_I',num2str(input_num),'(Nm).Par_Chain = MH_Results(:,1:8);']);
        eval(['Rslts_ODE_I',num2str(input_num),'(Nm).Fun_Chain = MH_Results(:,9);']);
    end
end

%% Task 3: Calculation of means and covariances from ODE MH Search for All Three Inputs
% Calculate means and covariances for all parameters from MH search of
% parameter space with an ODE objective function (chi-squared function). 
% These are calculated for the MH searches for all three experiments (each 
% with a different input signal).

for Nm = 1:4  % Loop through all four models
    % Load and save parameters from MH result files in workspace for mean and covariance calculations
    for i = 1:length(Par)
        for input_num = 1:3   % Loop through all three input signals
            STR_ODE_MH = ['MetHast_',num2str(Nm),'_ODE_Results_I',num2str(input_num),'.csv'];  
            MH_ODE_Results = csvread(STR_ODE_MH);
            eval(['Rslts_ODE_I',num2str(input_num),'(Nm).Par_Chain = MH_ODE_Results(:,1:8);']);
            eval(['Rslts_ODE_I',num2str(input_num),'(Nm).Fun_Chain = MH_ODE_Results(:,9);']);
            eval(['Rslts_ODE_I',num2str(input_num),'(Nm).Par',num2str(i),' = ','Rslts_ODE_I',num2str(input_num),'(Nm).Par_Chain(:,i);']);
            eval(['I',num2str(input_num),'_ODE','_M',num2str(Nm),'_Par',num2str(i),' = ','Rslts_ODE_I',num2str(input_num),'(Nm).Par',num2str(i),';']);
        end
    end

    % Calculate means and covariances for all parameters from MH search results for
    % each input. Excludes burn-in, which is estimated as 1/4 of total chain length.
    for i = 1:length(Par)
        for j = 1:length(Par)
            for input_num = 1:3   % Loop through all three input signals
                A = eval(['I',num2str(input_num),'_ODE','_M',num2str(Nm),'_Par',num2str(i),';']); 
                B = eval(['I',num2str(input_num),'_ODE','_M',num2str(Nm),'_Par',num2str(j),';']);
                MN=mean([A(floor(end/4):end),B(floor(end/4):end)]); COV=cov(A(floor(end/4):end),B(floor(end/4):end));
                eval(['RSLTS_ODE_M',num2str(Nm),'_I',num2str(input_num),'(i,j).MN = MN;']); 
                eval(['RSLTS_ODE_M',num2str(Nm),'_I',num2str(input_num),'(i,j).COV = COV;']);
            end
        end
    end

    % Generate diagonal matrix with all covariance results for each input. 
    % Excludes burn-in, which is estimated as 1/4 of total chain length.
    for input_num = 1:3  % Loop through all three input signals
        eval(['cov_all_M',num2str(Nm),'_I',num2str(input_num),...
            '_act = blkdiag(cov(Rslts_ODE_I',num2str(input_num),'(Nm).Par_Chain(floor(end/4):end,:)),0);']);
        eval(['cov_all_M',num2str(Nm),'_I',num2str(input_num),...
            ' = blkdiag(cov(Rslts_ODE_I',num2str(input_num),'(Nm).Par_Chain(floor(end/4):end,:)),0);']);
        eval(['cov_all_M',num2str(Nm),'_I',num2str(input_num),'(end,1:2) = [1 -1]*-0.05;']);
        eval(['cov_all_M',num2str(Nm),'_I',num2str(input_num),' = max(cov_all_M',num2str(Nm),'_I',num2str(input_num),',-0.05);']);
        eval(['cov_all_M',num2str(Nm),'_I',num2str(input_num),' = min(cov_all_M',num2str(Nm),'_I',num2str(input_num),',0.05);']);
    end
end

%% Task 3: FYI: Plots of All Parameter Combinations from ODE MH Search for All Three Inputs
% Excludes burn-in, which is estimated to be 1/4 of total chain length.

for Nm = 1:4
    figure()
    for ip1 = 1:7
        for ip2 = ip1+1:8
            subplot(7,7,(ip1-1)*7+ip2-1)
            plot(Rslts_ODE_I1(Nm).Par_Chain(floor(end/4):end,ip2),Rslts_ODE_I1(Nm).Par_Chain(floor(end/4):end,ip1));hold on
            plot(Rslts_ODE_I2(Nm).Par_Chain(floor(end/4):end,ip2),Rslts_ODE_I2(Nm).Par_Chain(floor(end/4):end,ip1));
            plot(Rslts_ODE_I3(Nm).Par_Chain(floor(end/4):end,ip2),Rslts_ODE_I3(Nm).Par_Chain(floor(end/4):end,ip1));
            plot(log10(abs(Parameters_True(ip2))),log10(abs(Parameters_True(ip1))),'o',...
               'markerfacecolor',[0.4,0.4,0.4],'MarkerEdgeColor',[0.3,0.3,0.3]);
            xlabel(sprintf('Par %i',ip2))
            ylabel(sprintf('Par %i',ip1))
        end
    end
    suptitle('ODE Metropolis-Hastings Results: Trajectory of Searches') %labels the set of scatter plots

    figure()
    for i=1:length(Par)-1
        for j=i+1:length(Par)
            subplot(7,7,(i-1)*7+j-1)
            scatter(Rslts_ODE_I1(Nm).Par_Chain(floor(end/4):end,j),Rslts_ODE_I1(Nm).Par_Chain(floor(end/4):end,i),10); hold on 
            scatter(Rslts_ODE_I2(Nm).Par_Chain(floor(end/4):end,j),Rslts_ODE_I2(Nm).Par_Chain(floor(end/4):end,i),10);
            scatter(Rslts_ODE_I3(Nm).Par_Chain(floor(end/4):end,j),Rslts_ODE_I3(Nm).Par_Chain(floor(end/4):end,i),10);
            COL = 'b'; conf = 0.90; %For ellipse: color and confidence interval
            [~,XX_1,YY_1] = error_ellipse_fill(RSLTS_ODE_I1(j,i).COV,RSLTS_ODE_I1(j,i).MN,'style',COL,'conf',conf);
            COL = 'g'; conf = 0.90; %For ellipse: color and confidence interval
            [~,XX_2,YY_2] = error_ellipse_fill(RSLTS_ODE_I2(j,i).COV,RSLTS_ODE_I2(j,i).MN,'style',COL,'conf',conf);
            COL = 'r'; conf = 0.90; %For ellipse: color and confidence interval
            [~,XX_3,YY_3] = error_ellipse_fill(RSLTS_ODE_I3(j,i).COV,RSLTS_ODE_I3(j,i).MN,'style',COL,'conf',conf);
            plot(log10(abs(Parameters_True(ip2))),log10(abs(Parameters_True(ip1))),'o',...
                    'markerfacecolor',[0.4,0.4,0.4],'MarkerEdgeColor',[0.3,0.3,0.3]);
            xlabel(sprintf('Par %i',j))
            ylabel(sprintf('Par %i',i))
        end
    end
    suptitle('ODE Metropolis-Hastings Results: Scatter') %labels the set of scatter plots
end

%% Task 3: FYI: Plot of ODE Objective Function Evaluations from MH Search for All Three Inputs
% Plot of full chain, including burn-in period.

for Nm = 1:4
    figure()
    ode1 = plot(1:N_Chain*Npart_ODE,Rslts_ODE_I1(Nm).Fun_Chain,'b'); hold on
    ode2 = plot(1:N_Chain*Npart_ODE,Rslts_ODE_I2(Nm).Fun_Chain,'g');
    ode3 = plot(1:N_Chain*Npart_ODE,Rslts_ODE_I3(Nm).Fun_Chain,'r');
    xlabel('MH Chain Length', 'FontWeight', 'bold')
    ylabel('Objective Function Value', 'FontWeight', 'bold')
    title('Objective Function Results from Metropolis-Hastings','FontWeight','bold')
    legend([ode1,ode2,ode3],{'Input 1','Input 2','Input3'},'FontSize',11,'Location','SouthEast') 
end

%% Task 3: Find two parameter sets from MH ODE Objective Function Evaluations (using Input 1)
% Acceptable range of objective function values include all values between
% abs(max likelihood) and abs(max likelihood - alpha). 
% Alpha is arbitrarily chosen after visual inspection of objective function evaluations.

Nm = 2; % For Model 2
[mx] = max(Rslts_ODE_I1(Nm).Fun_Chain(floor(end/4):end)); % Max likelihood
alpha = 4; range = abs(mx-alpha); % Range of acceptable obj fun values

% Location of acceptable obj fun values (in the full chain)
Igood = find(range > (abs(Rslts_ODE_I1(Nm).Fun_Chain)));
OF_good = Rslts_ODE_I1(Nm).Fun_Chain(Igood); % Acceptable obj fun values

% Location of min and max from acceptable obj fun values (in the full chain)
I_min = Igood(OF_good==min(OF_good));
I_max = Igood(OF_good==max(OF_good));

% Two acceptable parameter sets farthest apart in chain
Par_set_1 = Rslts_ODE_I1(Nm).Par_Chain(I_min(1),:);
Par_set_2 = Rslts_ODE_I1(Nm).Par_Chain(I_max(1),:); 

%% Task 4: Run 1000 SSA trajectories using Input 1
% Run many SSA trajectories using ODE best fit parameters, Par_set_1 (from ODE MH), and Par_set_2 (from ODE MH)

% For All SSA Runs
Nm = 2; Input = Input1; % Model 2 and Input 1
Output_Times = linspace(0,100,100); % Reset output times
Num_Runs = 1000;  % Number of SSA trajectories
Parameters = Model_Results_Parameters_I1(Nm,:); % Parameters from ODE best fit for specific model

% Multiple runs of many SSA trajectories with Best Fit Parameters (from ODE MH)
for i=1:4
    Parameters = Model_Results_Parameters_I1(Nm,:); % Parameters from ODE best fit for specific model
    ModFun = @(t,x)SSA_Mod(t,x,Parameters,Input,Nm); % Function with stoichoimetries and propensity functions
    eval(['m_RNA_Array_T',num2str(i),' = SSAHistfun(Num_Runs,ModFun,x0,Output_Times);']) % Run many SSA trajectories
    eval(['R_max_T',num2str(i),' = max(m_RNA_Array_T',num2str(i),');']); 
    eval(['R_max_T',num2str(i),'= max(R_max_T',num2str(i),');']) % Max RNA number from SSA
end

% % Run SSA with Parameter Set 1 (from ODE MH)
% Parameters = 10.^Par_set_1; % Parameter set 1 from ODE MH
% ModFun = @(t,x)SSA_Mod(t,x,Parameters,Input1,Nm); % Function with stoichoimetries and propensity functions
% m_RNA_Array1 = SSAHistfun(Num_Runs,ModFun,x0,Output_Times); % Run SSA
% R_max1 = max(m_RNA_Array1); R_max1 = max(R_max1); % Max RNA number from SSA
% 
% % Run SSA with Parameter Set 2 (from ODE MH)
% Parameters = 10.^Par_set_2; % Parameter set 1 from ODE MH
% ModFun = @(t,x)SSA_Mod(t,x,Parameters,Input1,Nm); % Function with stoichoimetries and propensity functions
% m_RNA_Array2 = SSAHistfun(Num_Runs,ModFun,x0,Output_Times); % Run SSA
% R_max2 = max(m_RNA_Array2); R_max2 = max(R_max2); % Max RNA number from SSA

%% Task 5: Run FSP MH search with Best Fit Parameters for All Three Inputs
% Generate A matrix and run FSP to obtain mRNA distributions at every time
% point for all three experiments (each with a different input signal)

Nf = 10; % Number of time points for FSP MH search
Output_Times = linspace(0,100,Nf); % Output time array
N = 80;  % First guess of the max number of RNA
N_states = 3; % Number of states

% Npart_FSP = # MH runs with chain length of N_Chain
% As described above in ODE MH section, total chain length is N_Chain*N_Thin*Npart_FSP
N_Chain = 500; N_Burn = 0; N_Thin = 10;  
Npart_FSP=200;

% GENERATE A MATRIX
% Generate components of A matrix for states 1,2,3
I1 = repmat([1 0 0]',N+1,1); I2 = repmat([0 1 0]',N+1,1); I3 = repmat([0 0 1]',N+1,1);
Mats_A.A12 = -spdiags(I1,0,N_states*(N+1),N_states*(N+1))+spdiags(I1(1:end-1),-1,N_states*(N+1),N_states*(N+1));
Mats_A.A23 = -spdiags(I2,0,N_states*(N+1),N_states*(N+1))+spdiags(I2(1:end-1),-1,N_states*(N+1),N_states*(N+1));
Mats_A.A21 = -spdiags(I2,0,N_states*(N+1),N_states*(N+1))+spdiags(I2,1,N_states*(N+1),N_states*(N+1));
Mats_A.A32 = -spdiags(I3,0,N_states*(N+1),N_states*(N+1))+spdiags(I3,1,N_states*(N+1),N_states*(N+1));

% Generate components of A matrix for RNA degradation
Ig = zeros(N_states*(N+1),1);Ig(1:3:end,1) = [0:N]';Ig(2:3:end,1) = [0:N]';Ig(3:3:end,1) = [0:N]';
Mats_A.Ag = -spdiags(Ig,0,N_states*(N+1),N_states*(N+1))+spdiags(Ig,3,N_states*(N+1),N_states*(N+1));

% Generate components of A matrix for RNA production
Mats_A.Ar2 = -spdiags(I2,0,N_states*(N+1),N_states*(N+1))+spdiags(I2,-3,N_states*(N+1),N_states*(N+1));
Mats_A.Ar3 = -spdiags(I3,0,N_states*(N+1),N_states*(N+1))+spdiags(I3,-3,N_states*(N+1),N_states*(N+1));

for Nm = 1:4  % Loop through all four models 
    for input_num = 1:3 % Loop through all three input signals
        
        % Load parameters from ODE best fit for all three input signals
        Init_Param = eval(['Model_Results_Parameters_I',num2str(input_num),'(Nm,:);']);
        eval(['Input = Input',num2str(input_num),';']);
        
        % Run the MH search with an FSP objective function
        run(RunMetHastFSP) 

        % Store parameter matrix and obj function chain from MH search
        MH_FSP_Results = csvread(STR_FSP_MH);
        eval(['Rslts_FSP_I',num2str(input_num),'(Nm).Par_Chain = MH_FSP_Results(:,1:8);']);
        eval(['Rslts_FSP_I',num2str(input_num),'(Nm).Fun_Chain = MH_FSP_Results(:,9);']);
    end
end

%% Task 5: Calculation of means and covariances from FSP MH Search for All Three Inputs
% Calculate means and covariances for all parameters from MH search of
% parameter space with an FSP objective function. These are calculated for
% the MH searches for all three experiments (each with a different input
% signal).

for Nm = 1:4   % Loop through all four models
% Load and save parameters from MH result files in workspace for mean and covariance calculations
    for i = 1:length(Par)
        for input_num = 1:3   % Loop through all three input signals
            STR_FSP_MH = ['MetHast_',num2str(Nm),'_FSP_Results_I',num2str(input_num),'.csv'];  
            MH_FSP_Results = csvread(STR_FSP_MH);
            eval(['Rslts_FSP_I',num2str(input_num),'(Nm).Par_Chain = MH_FSP_Results(:,1:8);']);
            eval(['Rslts_FSP_I',num2str(input_num),'(Nm).Fun_Chain = MH_FSP_Results(:,9);']);
            eval(['Rslts_FSP_I',num2str(input_num),'(Nm).Par',num2str(i),' = ','Rslts_FSP_I',num2str(input_num),'(Nm).Par_Chain(:,i);']);
            eval(['I',num2str(input_num),'_FSP','_M',num2str(Nm),'_Par',num2str(i),' = ','Rslts_FSP_I',num2str(input_num),'(Nm).Par',num2str(i),';']);
        end
    end

% Calculate means and covariances for all parameters from MH search results for
% each input. Excludes burn-in, which is estimated as 1/4 of total chain length.
    for i = 1:length(Par)
        for j = 1:length(Par)
            for input_num = 1:3   % Loop through all three input signals
                A = eval(['I',num2str(input_num),'_FSP','_M',num2str(Nm),'_Par',num2str(i),';']); 
                B = eval(['I',num2str(input_num),'_FSP','_M',num2str(Nm),'_Par',num2str(j),';']);
                MN=mean([A(floor(end/4):end),B(floor(end/4):end)]); COV=cov(A(floor(end/4):end),B(floor(end/4):end));
                eval(['RSLTS_FSP_M',num2str(Nm),'_I',num2str(input_num),'(i,j).MN = MN;']); 
                eval(['RSLTS_FSP_M',num2str(Nm),'_I',num2str(input_num),'(i,j).COV = COV;']);
            end
        end
    end

% Generate diagonal matrix with all covariance results for each input. 
% Excludes burn-in, which is estimated as 1/4 of total chain length.
    for input_num = 1:3  % Loop through all three input signals
        eval(['cov_all_fsp_M',num2str(Nm),'_I',num2str(input_num),...
            '_act = blkdiag(cov(Rslts_FSP_I',num2str(input_num),'(Nm).Par_Chain(floor(end/4):end,:)),0);']);
        eval(['cov_all_fsp_M',num2str(Nm),'_I',num2str(input_num),...
            ' = blkdiag(cov(Rslts_FSP_I',num2str(input_num),'(Nm).Par_Chain(floor(end/4):end,:)),0);']);
        eval(['cov_all_fsp_M',num2str(Nm),'_I',num2str(input_num),'(end,1:2) = [1 -1]*-0.05;']);
        eval(['cov_all_fsp_M',num2str(Nm),'_I',num2str(input_num),' = max(cov_all_fsp_M',num2str(Nm),'_I',num2str(input_num),',-0.05);']);
        eval(['cov_all_fsp_M',num2str(Nm),'_I',num2str(input_num),' = min(cov_all_fsp_M',num2str(Nm),'_I',num2str(input_num),',0.05);']);
    end
end

%% Task 5: FYI: Plots of All Parameter Combinations from FSP MH Search for All Three Inputs
% Excludes burn-in, which is estimated to be 1/4 of total chain length.

for Nm = 1:4  % Loop through all four models
    figure()
    for ip1 = 1:7
        for ip2 = ip1+1:8
            subplot(7,7,(ip1-1)*7+ip2-1)
            plot(Rslts_FSP_I1(Nm).Par_Chain(floor(end/4):end,ip2),Rslts_FSP_I1(Nm).Par_Chain(floor(end/4):end,ip1));hold on
            plot(Rslts_FSP_I2(Nm).Par_Chain(floor(end/4):end,ip2),Rslts_FSP_I2(Nm).Par_Chain(floor(end/4):end,ip1));
            plot(Rslts_FSP_I3(Nm).Par_Chain(floor(end/4):end,ip2),Rslts_FSP_I3(Nm).Par_ChainMH_FSP_Results_I3(floor(end/4):end,ip1));
            plot(log10(abs(Parameters_True(ip2))),log10(abs(Parameters_True(ip1))),'o',...
                'markerfacecolor',[0.4,0.4,0.4],'MarkerEdgeColor',[0.3,0.3,0.3]);
            xlabel(sprintf('Par %i',ip2))
            ylabel(sprintf('Par %i',ip1))
        end
    end
    suptitle('FSP Metropolis-Hastings Results: Trajectory of Searches') % Label the set of scatter plots

    figure()
    for i = 1:length(Parameters)-1
        for j = i+1:length(Parameters)
            subplot(7,7,(i-1)*7+j-1)
            scatter(Rslts_FSP_I1(Nm).Par_Chain(floor(end/4):end,j),Rslts_FSP_I1(Nm).Par_Chain(floor(end/4):end,i),10); hold on 
            scatter(Rslts_FSP_I2(Nm).Par_Chain(floor(end/4):end,j),Rslts_FSP_I2(Nm).Par_Chain(floor(end/4):end,i),10);
            scatter(Rslts_FSP_I3(Nm).Par_Chain(floor(end/4):end,j),Rslts_FSP_I3(Nm).Par_Chain(floor(end/4):end,i),10);
            plot(log10(abs(Parameters_True(j))),log10(abs(Parameters_True(i))),'o',...
                'markerfacecolor',[0.4,0.4,0.4],'MarkerEdgeColor',[0.3,0.3,0.3],'markersize',10);
            COL = 'b'; conf = 0.90; % For Input 1 ellipses: color and confidence interval
            [~,XX_1,YY_1] = error_ellipse_fill(RSLTS_FSP_I1(j,i).COV,RSLTS_FSP_I1(j,i).MN,'style',COL,'conf',conf);
            COL = 'g'; conf = 0.90; % For Input 2 ellipses: color and confidence interval
            [~,XX_2,YY_2] = error_ellipse_fill(RSLTS_FSP_I2(j,i).COV,RSLTS_FSP_I2(j,i).MN,'style',COL,'conf',conf);
            COL = 'r'; conf = 0.90; % For Input 3 ellipses: color and confidence interval
            [~,XX_3,YY_3] = error_ellipse_fill(RSLTS_FSP_I3(j,i).COV,RSLTS_FSP_I3(j,i).MN,'style',COL,'conf',conf);
            xlabel(sprintf('Par %i',j))
            ylabel(sprintf('Par %i',i))
        end
    end
    suptitle('FSP Metropolis-Hastings Results: Scatter') % Label the set of scatter plots
end

%% Task 5: FYI: Plot of FSP Objective Function Evaluations from MH Search for All Three Inputs
% Plot of full chain, including burn-in period.

for Nm = 1:4
    figure()
    f1=plot(1:N_Chain*Npart_FSP,Rslts_FSP_I1(Nm).Fun_Chain); hold on
    f2=plot(1:N_Chain*Npart_FSP,Rslts_FSP_I2(Nm).Fun_Chain);
    f3=plot(1:N_Chain*Npart_FSP,Rslts_FSP_I3(Nm).Fun_Chain);
    xlabel('MH Chain Length', 'FontWeight', 'bold')
    ylabel('Objective Function Value', 'FontWeight', 'bold')
    title('Objective Function Results from FSP Metropolis-Hastings','FontWeight','bold')
    legend([f1,f2,f3],{'Input 1','Input 2','Input3'},'FontSize',11,'Location','SouthEast') 
end

%% Task 5/6: Identify/Label ODE and FSP Best Fit Params For All Three Inputs
% Best fit parameters from ODE fitting, ODE MH search, and FSP MH search

for Nm = 1:4  % Loop through all four models
    for input_num = 1:3  % Loop through all three inputs
    % Best parameters from ODE fitting using all 3 input signals
        eval(['BestParams_ODE_M',num2str(Nm),'_I',num2str(input_num),' = log10(Model_Results_Parameters_I',num2str(input_num),'(Nm,:));']);
    % Best parameters from ODE Met-Haste using all 3 input signals
        eval(['[~,b',num2str(input_num),'] = max(Rslts_ODE_I',num2str(input_num),'(Nm).Fun_Chain);']); 
        eval(['BestParams_ODE_MH_M',num2str(Nm),'_I',num2str(input_num),' = Rslts_ODE_I',num2str(input_num),'(Nm).Par_Chain(b',num2str(input_num),',:);']);
    % Best parameters from FSP Met-Haste using all 3 input signals
        eval(['[~,b',num2str(input_num),'] = max(Rslts_FSP_I',num2str(input_num),'(Nm).Fun_Chain);']); 
        eval(['BestParams_FSP_M',num2str(Nm),'_I',num2str(input_num),' = Rslts_FSP_I',num2str(input_num),'(Nm).Par_Chain(b',num2str(input_num),',:);']);
    end
end

%% Task 5: Run FSP with Par Set 1 and Par Set 2 (from ODE MH), Best Fit FSP Params, and True Params For All Three Inputs 
% Run FSP for three experiments (each with a different input signal) to
% obtain full mRNA distributions at the same time points as the data.  This
% is done using both parameter sets obtained from the ODE MH search of
% parameter space, along with the best fit FSP params and true parameters for comparison.

for Nm = 1:4  % Loop through all four models
    for input_num = 1:3 % Loop through all three input signals
        eval(['Input = Input',num2str(input_num),';']);
    % FSP objective function to determine full mRNA distributions
        OBJ_FSP_MH_P = @(X)get_FSP_OBJ(10.^X,Input,Nm,Output_Times_data,Data_Set_Hists,Mats_A,N,N_states,'P_RNA');
    % Run FSP with Parameter Set 1 (from MH-ODE)
        eval(['P_RNA_1_M',num2str(Nm),'_I',num2str(input_num),' = OBJ_FSP_MH_P(Par_set_1);']);
    % Run FSP with Parameter Set 2 (from MH-ODE)
        eval(['P_RNA_2_M',num2str(Nm),'_I',num2str(input_num),' = OBJ_FSP_MH_P(Par_set_2);']);
    % Run FSP with True Parameters
        eval(['P_true_M',num2str(Nm),'_I',num2str(input_num),' = OBJ_FSP_MH_P(log10(abs(Parameters_True)));']);
    % Run FSP with FSP Best Fit Parameters
        eval(['P_FSP_M',num2str(Nm),'_I',num2str(input_num),' = OBJ_FSP_MH_P(BestParams_FSP_M2_I1);']);
    end
end

%% Task 6: Predictions - ODE PLOTS - TWO NEW SINUSOIDAL INPUT SIGNALS - Best Fit Params For All Three Inputs
% Plot mean mRNA (ODE solution) for two new sinusoidal input signals: 
% One with a fast frequency (Input_P1) and one with a slow frequency (Input_P2) 
% using parameters from ODE best fit and FSP best fit for all three inputs

Input_P1 = @(t)(1-cos(2*pi/10*t))*(t>5)*(t<70);   % Prediction Input signal 1 (fast)
Input_P2 = @(t)(1-cos(2*pi/60*t))*(t>5)*(t<70);   % Prediction Input signal 2 (slow)
Output_Times = linspace(0,100,100);

for Nm = 1:4 % Loop through all four models
    for input_num = 1:3  % Loop through best fit parameters from all three inputs
    % ODE Solution for Prediction Input 1
        [MRNA_P1_M2_True] = odefun(abs(Parameters_True),2,Input_P1,Output_Times,x0); %Solve ODEs for mean mRNA
        eval(['[MRNA_P1_ODE_M',num2str(Nm),'_I',num2str(input_num),'] = odefun(Model_Results_Parameters_I',num2str(input_num),'(Nm,:),Nm,Input_P1,Output_Times,x0);']);
        eval(['[MRNA_P1_FSP_M',num2str(Nm),'_I',num2str(input_num),'] = odefun(10.^BestParams_FSP_M',num2str(Nm),'_I',num2str(input_num),',Nm,Input_P1,Output_Times,x0);']);
    % ODE Solution for Prediction Input 2
        [MRNA_P2_M2_True] = odefun(abs(Parameters_True),2,Input_P2,Output_Times,x0); %Solve ODEs for mean mRNA
        eval(['[MRNA_P2_ODE_M',num2str(Nm),'_I',num2str(input_num),'] = odefun(Model_Results_Parameters_I',num2str(input_num),'(Nm,:),Nm,Input_P2,Output_Times,x0);']);
        eval(['[MRNA_P2_FSP_M',num2str(Nm),'_I',num2str(input_num),'] = odefun(10.^BestParams_FSP_M',num2str(Nm),'_I',num2str(input_num),',Nm,Input_P2,Output_Times,x0);']);
    end
end

% Make adjustments to ODE plots here
COL_TRUE = 'k'; COL1 = [0.3,0.3,0.3]; COL2 = [0.5,0.5,0.5]; COL3 = [0.75,0.75,0.75];
STYLE_TRUE = '-'; STYLE1 = '-'; STYLE2 = ':'; 
LOC = 'SouthEast'; XLIM = 60; LIN = 2; 

% Plot mean mRNA (ODE solution) of a new input (fast) signal (Input_P1) using parameters 
% from ODE best fit and FSP best fit for all three inputs and Model 2
figure()
subplot(1,2,1)
plot(Output_Times,MRNA_P1_M2_True,'Color',COL_TRUE,'linewidth',LIN,'LineStyle',STYLE_TRUE); hold on % Plot ODE soln for True Parameters
plot(Output_Times,MRNA_P1_ODE_M2_I1,'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1); % Plot ODE soln for ODE Best Parameters: Input 1
plot(Output_Times,MRNA_P1_ODE_M2_I2,'Color',COL2,'linewidth',LIN,'LineStyle',STYLE1); % Plot ODE soln for ODE Best Parameters: Input 2
plot(Output_Times,MRNA_P1_ODE_M2_I3,'Color',COL3,'linewidth',LIN,'LineStyle',STYLE1); % Plot ODE soln for ODE Best Parameters: Input 3
plot(Output_Times,MRNA_P1_FSP_M2_I1,'Color',COL1,'linewidth',LIN,'LineStyle',STYLE2); % Plot ODE soln for FSP Best Parameters: Input 1
plot(Output_Times,MRNA_P1_FSP_M2_I2,'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2); % Plot ODE soln for FSP Best Parameters: Input 2
plot(Output_Times,MRNA_P1_FSP_M2_I3,'Color',COL3,'linewidth',LIN,'LineStyle',STYLE2); % Plot ODE soln for FSP Best Parameters: Input 3
title('Input Signal with Fast Frequency') % Adding title
xlabel('Time', 'FontWeight', 'bold') % Adding label to x axis
ylabel('Average mRNA Count','FontWeight', 'bold') % Adding label to y axis
legend({'M2:True Parameters','M2:I1:ODE','M2:I2:ODE','M2:I3:ODE','M2:I1:FSP','M2:I2:FSP','M2:I3:FSP'},...
    'FontSize',12,'Location',LOC); grid off;

% Plot mean mRNA (ODE solution) of a new input (slow) signal (Input_P2) using parameters 
% from ODE best fit and FSP best fit for all three inputs and Model 2
subplot(1,2,2)
plot(Output_Times,MRNA_P2_M2_True,'Color',COL_TRUE,'linewidth',LIN,'LineStyle',STYLE_TRUE); hold on % Plot ODE soln for True Parameters
plot(Output_Times,MRNA_P2_ODE_M2_I1,'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1); % Plot ODE soln for ODE Best Parameters: Input 1
plot(Output_Times,MRNA_P2_ODE_M2_I2,'Color',COL2,'linewidth',LIN,'LineStyle',STYLE1); % Plot ODE soln for ODE Best Parameters: Input 2
plot(Output_Times,MRNA_P2_ODE_M2_I3,'Color',COL3,'linewidth',LIN,'LineStyle',STYLE1); % Plot ODE soln for ODE Best Parameters: Input 3
plot(Output_Times,MRNA_P2_FSP_M2_I1,'Color',COL1,'linewidth',LIN,'LineStyle',STYLE2); % Plot ODE soln for FSP Best Parameters: Input 1
plot(Output_Times,MRNA_P2_FSP_M2_I2,'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2); % Plot ODE soln for FSP Best Parameters: Input 2
plot(Output_Times,MRNA_P2_FSP_M2_I3,'Color',COL3,'linewidth',LIN,'LineStyle',STYLE2); % Plot ODE soln for FSP Best Parameters: Input 3
title('Input Signal with Slow Frequency') % Adding title
xlabel('Time', 'FontWeight', 'bold') % Adding label to x axis
ylabel('Average mRNA Count','FontWeight', 'bold') % Adding label to y axis
legend({'M2:True Parameters','M2:I1:ODE','M2:I2:ODE','M2:I3:ODE','M2:I1:FSP','M2:I2:FSP','M2:I3:FSP'},...
    'FontSize',12,'Location',LOC); grid off;
suptitle('Experimental Prediction of Mean mRNA Level: ODE Solution')

%% Task 6: Predictions - FSP PLOTS - FOUR NEW SINUSOIDAL INPUT SIGNALS - Best Fit Params For Input 1
% Plot full distributions (FSP solution) for four new sinusoidal input signals 
% using parameters from FSP best fit and true parameters for all three comparison inputs: 
% Input_P1: fast frequency, Input_P2: slow frequency 
% Input_P3 and Input_P4: mid-level frequencies  

Output_Times_expt = linspace(0,100,11);
Input_P1 = @(t)(1-cos(2*pi/10*t))*(t>5)*(t<70);   % Prediction Input signal 1 (fast)
Input_P2 = @(t)(1-cos(2*pi/60*t))*(t>5)*(t<70);   % Prediction Input signal 2 (slow)
Input_P3 = @(t)(1-cos(2*pi/20*t))*(t>5)*(t<70);   % Prediction Input signal 3 (mid/fast)
Input_P4 = @(t)(1-cos(2*pi/50*t))*(t>5)*(t<70);   % Prediction Input signal 4 (mid/slow)

for Nm = 1:4  % Loop through all four models
    for input_num = 1:3 % Loop through all three input signals
    % FSP for Prediction Input 1 
        OBJ_FSP_MH_P = @(X)get_FSP_OBJ(10.^X,Input_P1,Nm,Output_Times_expt,Data_Set_Hists,Mats_A,N,N_states,'P_RNA');
        eval(['P_true_P1_M',num2str(Nm),' = OBJ_FSP_MH_P(log10(abs(Parameters_True)));']);
        eval(['P_RNA_P1_M',num2str(Nm),'_I',num2str(input_num),' = OBJ_FSP_MH_P(BestParams_FSP_M',num2str(Nm),'_I',num2str(input_num),');']);
    % FSP for Prediction Input 2
        OBJ_FSP_MH_P = @(X)get_FSP_OBJ(10.^X,Input_P2,Nm,Output_Times_expt,Data_Set_Hists,Mats_A,N,N_states,'P_RNA');
        eval(['P_true_P2_M',num2str(Nm),' = OBJ_FSP_MH_P(log10(abs(Parameters_True)));']);
        eval(['P_RNA_P2_M',num2str(Nm),'_I',num2str(input_num),' = OBJ_FSP_MH_P(BestParams_FSP_M',num2str(Nm),'_I',num2str(input_num),');']);
    % FSP for Prediction Input 3
        OBJ_FSP_MH_P = @(X)get_FSP_OBJ(10.^X,Input_P3,Nm,Output_Times_expt,Data_Set_Hists,Mats_A,N,N_states,'P_RNA');
        eval(['P_true_P3_M',num2str(Nm),' = OBJ_FSP_MH_P(log10(abs(Parameters_True)));']);
        eval(['P_RNA_P3_M',num2str(Nm),'_I',num2str(input_num),' = OBJ_FSP_MH_P(BestParams_FSP_M',num2str(Nm),'_I',num2str(input_num),');']);
    % FSP for Prediction Input 4
        OBJ_FSP_MH_P = @(X)get_FSP_OBJ(10.^X,Input_P4,Nm,Output_Times_expt,Data_Set_Hists,Mats_A,N,N_states,'P_RNA');
        eval(['P_true_P4_M',num2str(Nm),' = OBJ_FSP_MH_P(log10(abs(Parameters_True)));']);
        eval(['P_RNA_P4_M',num2str(Nm),'_I',num2str(input_num),' = OBJ_FSP_MH_P(BestParams_FSP_M',num2str(Nm),'_I',num2str(input_num),');']);
    end
end

% Make adjustments to FSP plots here
COL_TRUE = 'k'; COL1 = [0.2,0.2,0.2]; COL2 = [0.4,0.4,0.4]; COL3 = [0.6,0.6,0.6]; COL4 = [0.8,0.8,0.8]; 
STYLE_TRUE = '-'; STYLE1 = '-'; STYLE2 = '-'; STYLE3 = '-'; STYLE4 = '-';
LOC = 'NorthEast'; XLIM = 60; LIN = 1.5; 

% PLOT FSP DISTRIBUTION RESULTS: Prediction Input Signal 1 and 2
figure()
    for i = 1:10
    % Plot of FSP distributions for all four models compared to true
    % parameters for Prediction Input Signal 1
        subplot(2,10,i)
        plot([0:XLIM], P_true_P1_M2(i+1,1:(XLIM+1)),'Color',COL_TRUE,'linewidth',LIN,'LineStyle',STYLE_TRUE);hold on;
        plot([0:XLIM], P_RNA_P1_M1_I1(i+1,1:(XLIM+1)),'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1);
        plot([0:XLIM], P_RNA_P1_M2_I1(i+1,1:(XLIM+1)),'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2);
        plot([0:XLIM], P_RNA_P1_M3_I1(i+1,1:(XLIM+1)),'Color',COL3,'linewidth',LIN,'LineStyle',STYLE3);
        plot([0:XLIM], P_RNA_P1_M4_I1(i+1,1:(XLIM+1)),'Color',COL4,'linewidth',LIN,'LineStyle',STYLE4);
        set(gca,'xlim',[0,XLIM],'ylim',[0,0.15])
        title(sprintf('t = %i min',round(Output_Times_expt(i+1))),'FontWeight','bold','FontSize',12)
        legend({'M2:TP','M1:I1:FSP','M2:I1:FSP','M3:I1:FSP','M4:I1:FSP'},'FontSize',8,'Location',LOC)
    % Plot of FSP distributions for all four models compared to true
    % parameters for Prediction Input Signal 2
        subplot(2,10,(i+10))
        plot([0:XLIM], P_true_P2_M2(i+1,1:(XLIM+1)),'Color',COL_TRUE,'linewidth',LIN,'LineStyle',STYLE_TRUE);hold on;
        plot([0:XLIM], P_RNA_P2_M1_I1(i+1,1:(XLIM+1)),'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1);
        plot([0:XLIM], P_RNA_P2_M2_I1(i+1,1:(XLIM+1)),'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2);
        plot([0:XLIM], P_RNA_P2_M3_I1(i+1,1:(XLIM+1)),'Color',COL3,'linewidth',LIN,'LineStyle',STYLE3);
        plot([0:XLIM], P_RNA_P2_M4_I1(i+1,1:(XLIM+1)),'Color',COL4,'linewidth',LIN,'LineStyle',STYLE4);
        set(gca,'xlim',[0,XLIM],'ylim',[0,0.15])
        legend({'M2:TP','M1:I1:FSP','M2:I1:FSP','M3:I1:FSP','M4:I1:FSP'},'FontSize',8,'Location',LOC)
    end    
subplot(2,10,1) 
ylabel('Probability','FontWeight','bold','FontSize',14) %Adding label to y axis
subplot(2,10,11) 
ylabel('Probability','FontWeight','bold','FontSize',14) %Adding label to y axis
[ax1,h1]=suplabel('mRNA Level','x',[0.09 0.09 0.84 0.84]); set(h1,'FontWeight','bold','FontSize',14) 
suptitle({'Distribution Predictions for New Sinusoidal Input Signal at Different Times',...
    '(Top Row) Fast Frequency (Bottom Row) Slow Frequency'}) %Label the set of scatter plots

% PLOT FSP DISTRIBUTION RESULTS: Prediction Input Signal 3 and 4
figure()
    for i = 1:10
    % Plot of FSP distributions for all four models compared to true
    % parameters for Prediction Input Signal 3
        subplot(2,10,i)
        plot([0:XLIM], P_true_P3_M2(i+1,1:(XLIM+1)),'Color',COL_TRUE,'linewidth',LIN,'LineStyle',STYLE_TRUE);hold on;
        plot([0:XLIM], P_RNA_P3_M1_I1(i+1,1:(XLIM+1)),'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1);
        plot([0:XLIM], P_RNA_P3_M2_I1(i+1,1:(XLIM+1)),'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2);
        plot([0:XLIM], P_RNA_P3_M3_I1(i+1,1:(XLIM+1)),'Color',COL3,'linewidth',LIN,'LineStyle',STYLE3);
        plot([0:XLIM], P_RNA_P3_M4_I1(i+1,1:(XLIM+1)),'Color',COL4,'linewidth',LIN,'LineStyle',STYLE4);
        set(gca,'xlim',[0,XLIM],'ylim',[0,0.15])
        title(sprintf('t = %i min',round(Output_Times_expt(i+1))),'FontWeight','bold','FontSize',12)
        legend({'M2:TP','M1:I1:FSP','M2:I1:FSP','M3:I1:FSP','M4:I1:FSP'},'FontSize',8,'Location',LOC)
    % Plot of FSP distributions for all four models compared to true
    % parameters for Prediction Input Signal 4
        subplot(2,10,(i+10))
        plot([0:XLIM], P_true_P4_M2(i+1,1:(XLIM+1)),'Color',COL_TRUE,'linewidth',LIN,'LineStyle',STYLE_TRUE);hold on;
        plot([0:XLIM], P_RNA_P4_M1_I1(i+1,1:(XLIM+1)),'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1);
        plot([0:XLIM], P_RNA_P4_M2_I1(i+1,1:(XLIM+1)),'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2);
        plot([0:XLIM], P_RNA_P4_M3_I1(i+1,1:(XLIM+1)),'Color',COL3,'linewidth',LIN,'LineStyle',STYLE3);
        plot([0:XLIM], P_RNA_P4_M4_I1(i+1,1:(XLIM+1)),'Color',COL4,'linewidth',LIN,'LineStyle',STYLE4);
        set(gca,'xlim',[0,XLIM],'ylim',[0,0.15])
        legend({'M2:TP','M1:I1:FSP','M2:I1:FSP','M3:I1:FSP','M4:I1:FSP'},'FontSize',8,'Location',LOC)
    end    
subplot(2,10,1) 
ylabel('Probability','FontWeight','bold','FontSize',14) %Adding label to y axis
subplot(2,10,11) 
ylabel('Probability','FontWeight','bold','FontSize',14) %Adding label to y axis
[ax1,h1]=suplabel('mRNA Level','x',[0.09 0.09 0.84 0.84]); set(h1,'FontWeight','bold','FontSize',14) 
suptitle({'Distribution Predictions for New Sinusoidal Input Signal at Different Times',...
    '(Top Row) Mid-Fast Frequency (Bottom Row) Mid-Slow Frequency'}) %Label the set of scatter plots

%% Task 6: Predictions - FSP PLOTS - MANY DIFFERENT INPUT SIGNALS - Best Fit Params For All Three Inputs
% Plot full distributions (FSP solution) for all four models and for many different sinusoidal input signals, 
% each with a different frequency, using parameters from FSP best fit and true parameters for all three inputs

Output_Times_expt = linspace(0,100,11);
for pd = [5:5:65]  % Loop through multiple different frequencies for comparison
    for Nm = 2:4 % Loop through all model numbers
        freq = 1/pd;  % Frequency for specific experiment input signal
        Input_Exp = @(t)(1-cos(2*pi*freq*t))*(t>5)*(t<70);   % Sinusoidal input signal
    % Run FSP with specific experiment input signal and true parameters
        OBJ_FSP_MH_P = @(X)get_FSP_OBJ(10.^X,Input_Exp,Nm,Output_Times_expt,Data_Set_Hists,Mats_A,N,N_states,'P_RNA');
        eval(['P_EXPT_true_M',num2str(Nm),' = OBJ_FSP_MH_P(log10(abs(Parameters_True)));']);
    % FSP for Predictions
        for input_num = 1:3  % Loop through FSP best fit parameters from all three original input signals 
            OBJ_FSP_MH_P = @(X)get_FSP_OBJ(10.^X,Input_Exp,Nm,Output_Times_expt,Data_Set_Hists,Mats_A,N,N_states,'P_RNA');
            eval(['P_RNA_EXPT_M',num2str(Nm),'_I',num2str(input_num),' = OBJ_FSP_MH_P(BestParams_FSP_M',num2str(Nm),'_I',num2str(input_num),');']);
        end
    end
end
    
% Make adjustments to FSP plots here
    COL_TRUE = 'k'; COL1 = [0.2,0.2,0.2]; COL2 = [0.4,0.4,0.4]; COL3 = [0.6,0.6,0.6]; COL4 = [0.8,0.8,0.8]; 
    STYLE_TRUE = '-'; STYLE1 = '-'; STYLE2 = '-'; STYLE3 = '-'; STYLE4 = '-';
    LOC = 'NorthEast'; XLIM = 60; LIN = 1.5; 
    
% PLOT FSP DISTRIBUTION RESULTS FOR ALL EXPERIMENTS (INPUT SIGNALS)
    figure()
    for i = 1:10
        subplot(2,5,i)
        plot([0:XLIM], P_EXPT_true_M2(i+1,1:(XLIM+1)),'Color',COL_TRUE,'linewidth',LIN,'LineStyle',STYLE_TRUE);hold on;
        plot([0:XLIM], P_RNA_EXPT_M1_I1(i+1,1:(XLIM+1)),'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1);
        plot([0:XLIM], P_RNA_EXPT_M2_I1(i+1,1:(XLIM+1)),'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2);
        plot([0:XLIM], P_RNA_EXPT_M3_I1(i+1,1:(XLIM+1)),'Color',COL3,'linewidth',LIN,'LineStyle',STYLE3);
        plot([0:XLIM], P_RNA_EXPT_M4_I1(i+1,1:(XLIM+1)),'Color',COL4,'linewidth',LIN,'LineStyle',STYLE4);        
        set(gca,'xlim',[0,XLIM],'ylim',[0,0.15])
        legend({'M2:TP','M1:I1:FSP','M2:I1:FSP','M3:I1:FSP','M4:I1:FSP'},'FontSize',9,'Location',LOC)
        title(sprintf('t = %i min',round(Output_Times_expt(i+1))))
    end
    subplot(2,5,1)
    ylabel('Probability','FontWeight','bold','FontSize',14) % Adding label to y axis
    subplot(2,5,6) 
    ylabel('Probability','FontWeight','bold','FontSize',14) % Adding label to y axis
    [ax1,h1]=suplabel('mRNA Level','x',[0.09 0.09 0.84 0.84]); set(h1,'FontWeight','bold','FontSize',14) 
    suptitle(sprintf('Frequencey = %i',pd)) % Label the set of scatter plots
end
% total_diff(pd/5) = sum(std(EXP_Var))

%% Task 6: Predictions for Figure 11 - FSP PLOTS - INPUT 2/3/4 - Best Fit Params For Input 1 
% Predict full distributions (FSP solution) for three experiments: step, ramp, and fast sinusoidal
% using parameters from FSP best fit for input signal 1 (true input)

Output_Times_expt = linspace(0,100,11);
Input_Exp_1 = Input2;   % Prediction Input signal 2 (step)
Input_Exp_2 = Input3;   % Prediction Input signal 3 (ramp)
Input_Exp_3 = @(t)(1-cos(2*pi/10*t))*(t>5)*(t<70);   % Prediction Input signal 4 (fast) 
    
for Nm = 1:4 % Loop through all four models
% Predict FSP distributions for experiments 1-3 using true parameters
    OBJ_FSP_MH_P_1 = @(X)get_FSP_OBJ(10.^X,Input_Exp_1,Nm,Output_Times_expt,Data_Set_Hists,Mats_A,N,N_states,'P_RNA');
    OBJ_FSP_MH_P_2 = @(X)get_FSP_OBJ(10.^X,Input_Exp_2,Nm,Output_Times_expt,Data_Set_Hists,Mats_A,N,N_states,'P_RNA');
    OBJ_FSP_MH_P_3 = @(X)get_FSP_OBJ(10.^X,Input_Exp_3,Nm,Output_Times_expt,Data_Set_Hists,Mats_A,N,N_states,'P_RNA');
    eval(['P_EXPT_true_M',num2str(Nm),'_1 = OBJ_FSP_MH_P_1(log10(abs(Parameters_True)));']);
    eval(['P_EXPT_true_M',num2str(Nm),'_2 = OBJ_FSP_MH_P_2(log10(abs(Parameters_True)));']);
    eval(['P_EXPT_true_M',num2str(Nm),'_3 = OBJ_FSP_MH_P_3(log10(abs(Parameters_True)));']);
% Predict FSP distributions for experiments 1-3 using FSP best fit parameters
    for input_num = 1 % Using FSP best fit parameters from original experiment with input signal 1
        OBJ_FSP_MH_P_1 = @(X)get_FSP_OBJ(10.^X,Input_Exp_1,Nm,Output_Times_expt,Data_Set_Hists,Mats_A,N,N_states,'P_RNA');
        eval(['P_RNA_EXPT_M',num2str(Nm),'_I',num2str(input_num),'_1 = OBJ_FSP_MH_P_1(BestParams_FSP_M',num2str(Nm),'_I',num2str(input_num),');']);
        OBJ_FSP_MH_P_2 = @(X)get_FSP_OBJ(10.^X,Input_Exp_2,Nm,Output_Times_expt,Data_Set_Hists,Mats_A,N,N_states,'P_RNA');
        eval(['P_RNA_EXPT_M',num2str(Nm),'_I',num2str(input_num),'_2 = OBJ_FSP_MH_P_2(BestParams_FSP_M',num2str(Nm),'_I',num2str(input_num),');']);
        OBJ_FSP_MH_P_3 = @(X)get_FSP_OBJ(10.^X,Input_Exp_3,Nm,Output_Times_expt,Data_Set_Hists,Mats_A,N,N_states,'P_RNA');
        eval(['P_RNA_EXPT_M',num2str(Nm),'_I',num2str(input_num),'_3 = OBJ_FSP_MH_P_3(BestParams_FSP_M',num2str(Nm),'_I',num2str(input_num),');']);
    end
end

% Make adjustments to FSP plots here
    COL_TRUE = 'k'; COL1 = [0.2,0.2,0.2]; COL2 = [0.4,0.4,0.4]; COL3 = [0.6,0.6,0.6]; COL4 = [0.8,0.8,0.8]; 
    STYLE_TRUE = '-'; STYLE1 = '-'; STYLE2 = '-'; STYLE3 = '-'; STYLE4 = '-';
    LOC = 'NorthEast'; XLIM = 60; LIN = 1.5; 

% PLOT FSP RESULTS USING INPUT 2
    figure()
    for i = 2:11
        subplot(2,5,(i-1))
        plot([0:XLIM], P_EXPT_true_M2_1(i,1:(XLIM+1)),'Color',COL_TRUE,'linewidth',LIN,'LineStyle',STYLE_TRUE);hold on;
        plot([0:XLIM], P_RNA_EXPT_M1_I1_1(i,1:(XLIM+1)),'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1);
        plot([0:XLIM], P_RNA_EXPT_M2_I1_1(i,1:(XLIM+1)),'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2);
        plot([0:XLIM], P_RNA_EXPT_M3_I1_1(i,1:(XLIM+1)),'Color',COL3,'linewidth',LIN,'LineStyle',STYLE3);
        plot([0:XLIM], P_RNA_EXPT_M4_I1_1(i,1:(XLIM+1)),'Color',COL4,'linewidth',LIN,'LineStyle',STYLE4);
        set(gca,'xlim',[0,XLIM],'ylim',[0,0.15])
        title(sprintf('t = %i min',round(Output_Times_expt(i))))
        legend({'M2:TP','M1:I1:FSP','M2:I1:FSP','M3:I1:FSP','M4:I1:FSP'},'FontSize',9,'Location',LOC)
    end
    subplot(2,5,1)
    ylabel('Probability','FontWeight','bold','FontSize',14) %Adding label to y axis
    subplot(2,5,6) 
    ylabel('Probability','FontWeight','bold','FontSize',14) %Adding label to y axis
    [ax1,h1]=suplabel('mRNA Level','x',[0.09 0.09 0.84 0.84]); set(h1,'FontWeight','bold','FontSize',14) 
    suptitle('Step Input Signal') % Label the set of scatter plots
    
% PLOT FSP RESULTS USING INPUT 3
    figure()
    for i = 2:11
        subplot(2,5,(i-1))
        plot([0:XLIM], P_EXPT_true_M2_2(i,1:(XLIM+1)),'Color',COL_TRUE,'linewidth',LIN,'LineStyle',STYLE_TRUE);hold on;
        plot([0:XLIM], P_RNA_EXPT_M1_I1_2(i,1:(XLIM+1)),'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1);
        plot([0:XLIM], P_RNA_EXPT_M2_I1_2(i,1:(XLIM+1)),'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2);
        plot([0:XLIM], P_RNA_EXPT_M3_I1_2(i,1:(XLIM+1)),'Color',COL3,'linewidth',LIN,'LineStyle',STYLE3);
        plot([0:XLIM], P_RNA_EXPT_M4_I1_2(i,1:(XLIM+1)),'Color',COL4,'linewidth',LIN,'LineStyle',STYLE4);
        set(gca,'xlim',[0,XLIM],'ylim',[0,0.15])
        title(sprintf('t = %i min',round(Output_Times_expt(i))))
        legend({'M2:TP','M1:I1:FSP','M2:I1:FSP','M3:I1:FSP','M4:I1:FSP'},'FontSize',9,'Location',LOC)
    end
    subplot(2,5,1)
    ylabel('Probability','FontWeight','bold','FontSize',14) % Adding label to y axis
    subplot(2,5,6) 
    ylabel('Probability','FontWeight','bold','FontSize',14) % Adding label to y axis
    [ax1,h1]=suplabel('mRNA Level','x',[0.09 0.09 0.84 0.84]); set(h1,'FontWeight','bold','FontSize',14) 
    suptitle('Ramp Input Signal') % Label the set of scatter plots
    
% PLOT FSP RESULTS USING INPUT 4
    figure()
    for i = 2:11
        subplot(2,5,(i-1))
        plot([0:XLIM], P_EXPT_true_M2_3(i,1:(XLIM+1)),'Color',COL_TRUE,'linewidth',LIN,'LineStyle',STYLE_TRUE);hold on;
        plot([0:XLIM], P_RNA_EXPT_M1_I1_3(i,1:(XLIM+1)),'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1);
        plot([0:XLIM], P_RNA_EXPT_M2_I1_3(i,1:(XLIM+1)),'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2);
        plot([0:XLIM], P_RNA_EXPT_M3_I1_3(i,1:(XLIM+1)),'Color',COL3,'linewidth',LIN,'LineStyle',STYLE3);
        plot([0:XLIM], P_RNA_EXPT_M4_I1_3(i,1:(XLIM+1)),'Color',COL4,'linewidth',LIN,'LineStyle',STYLE4);
        set(gca,'xlim',[0,XLIM],'ylim',[0,0.15])
        title(sprintf('t = %i min',round(Output_Times_expt(i))))
        legend({'M2:TP','M1:I1:FSP','M2:I1:FSP','M3:I1:FSP','M4:I1:FSP'},'FontSize',9,'Location',LOC)
    end
    subplot(2,5,1)
    ylabel('Probability','FontWeight','bold','FontSize',14) % Adding label to y axis
    subplot(2,5,6) 
    ylabel('Probability','FontWeight','bold','FontSize',14) % Adding label to y axis
    [ax1,h1]=suplabel('mRNA Level','x',[0.09 0.09 0.84 0.84]); set(h1,'FontWeight','bold','FontSize',14) 
    suptitle('Sinusoidal (Fast) Input Signal') % Label the set of scatter plots

%% FIGURES FOR CHAPTER: Figure 2: ODE SOLUTION - TWO ARBITRARY PARAMETER SETS FOR INPUT 1 AND MODEL 2
% Plot input signal 1 and ODE solution for average mRNA count for two arbitrary parameter sets
figure()
time = [0 100]; xval = Output_Times_data; % For plotting input signal
Nm = 2;  % Model Number

% Plot Input Signal 1: Sinusoidal
subplot(2,1,1)
yval = [0 3];
for i=1:length(xval)
    line([xval(i),xval(i)],[yval(1),yval(2)],'LineStyle',':',...
        'LineWidth',1,'Color',[0.3,0.3,0.3]);hold on
end
fI1 = ezplot(Input1,time);hold on
set(gca,'ylim',[0,3]); set(fI1,'Color','k','LineWidth',2)
xlabel('Time (min)','FontWeight','bold','FontSize',12); 
ylabel('Y_{1}(t)','FontWeight','bold','FontSize',12)
title('Input Signal 1: Sinusoidal','FontWeight','bold','FontSize',12)

% Plot ODE solution for average mRNA count for two arbitary parameter sets
subplot(2,1,2)
yval = [0 15];
for i=1:length(xval)
    line([xval(i),xval(i)],[yval(1),yval(2)],'LineStyle',':',...
        'LineWidth',1,'Color',[0.3,0.3,0.3]);hold on
end
Output_Times = linspace(0,100,100);  % Vector with output times specified
Parameters_A = [3, -0.4, 2.5, 1, 9, 2, 1, 0.8]; % Arbitrary set of parameters (A)
Parameters_B = [1, -0.6, 1, 0.2, 5, 2, 1, 0.2];  % Arbitrary set of parameters (B)
[M_A] = odefun(abs(Parameters_A),Nm,Input1,Output_Times,x0); % Solve ODEs for mean mRNA
[M_B] = odefun(abs(Parameters_B),Nm,Input1,Output_Times,x0); % Solve ODEs for mean mRNA
m1=plot(Output_Times,M_A,'Color',[0.2,0.2,0.2],'linewidth',2); hold on % Plot ODE soln for 1st parameter set
m2=plot(Output_Times,M_B,'Color',[0.8,0.8,0.8],'linewidth',2) % Plot ODE soln for 2nd parameter set on same plot
legend([m1,m2],{'Parameter Set A', 'Parameter Set B'},'FontSize',12, 'Location', 'SouthEast');  
xlabel('Time (min)', 'FontWeight','bold','FontSize',12) 
ylabel('Average mRNA Count','FontWeight','bold','FontSize',12)
title('Mean mRNA Level: ODE Solution','FontSize',12)

%% FIGURES FOR CHAPTER: Figure 3: ODE FITTING RESULTS FOR MODEL 2
% Includes two sets of error bars: large (light) error bars are std dev 
% of the data, small (dark) error bars are standard error of the mean
% Black astericks indicate the mean of the data

Nm = 2; Input = Input1; % Model and input can be changed.
Output_Times = linspace(0,100,100);  % Vector with output times specified

figure()
% Plot ODE fitting using Best Fit Parameters: Input 1
MVT_test = odefun(Model_Results_Parameters_I1(Nm,:),Nm,Input,Output_Times,x0);  
mvt = plot(Output_Times,MVT_test,'Color',[0.7,0.7,0.7],'linewidth',2); hold on
% Plot error bars for standard deviations and standard error
err = std(Data_Set)./sqrt(length(Data_Set)); % Standard error of the mean of the data
e1 = errorbar(Output_Times_data,mean(Data_Set),std(Data_Set),'Color',[0.7,0.7,0.7],'linewidth',0.8,'linestyle','none');
e2 = errorbar(Output_Times_data,mean(Data_Set),err,'Color',[0.3,0.3,0.3],'linewidth',1.5,'linestyle','none');
mvt2 = scatter(Output_Times_data,mean(Data_Set),'k','*'); % Plot mean of data
set(gca,'xlim',[Output_Times(1),Output_Times(end)]);
legend([e1,e2],{'Mean(data) \pm Std Deviation (light)','Mean(data) \pm Std Error (dark)'},'FontSize',12,'Location','SouthEast');
xlabel('Time (min)','FontWeight','bold','FontSize',14)
ylabel('Average mRNA Count','FontWeight','bold','FontSize',14)
title('Parameter Optimization (ODE Fitting)','FontWeight','bold','FontSize',14)

%% FIGURES FOR CHAPTER: Figure 4: ODE MET-HASTE RESULTS FOR INPUT 1 AND MODEL 2
% Excluding burn-in period (estimated to be 1/4 of total chain length)

Nm = 2;  % Model Number
figure()
subplot(1,2,1)  % Plot ODE MH scatter plot with 90% confidence ellipse
% Make a scatter plot of Par5 vs. Par8
mh1=scatter(I1_ODE_M2_Par5(floor(end/4):end),I1_ODE_M2_Par8(floor(end/4):end),10,...
    'MarkerEdgeColor',[0.7,0.7,0.7]); hold on 
% Call routine to plot ellipses for the 90% confidense intervals.
COL = [0.2,0.2,0.2]; conf = 0.90; % For ellipses: color and confidence interval  
[~,XX,YY] = error_ellipse_fill(RSLTS_ODE_M2_I1(5,8).COV,RSLTS_ODE_M2_I1(5,8).MN,...
    'style',COL,'conf',conf); 
% Plot parameter set 1 from ODE MH
plot(Par_set_1(5),Par_set_1(8),'wx','markersize',12,'linewidth',8); hold on   
mh2=plot(Par_set_1(5),Par_set_1(8),'kx','markersize',12,'linewidth',4); hold on
% Plot parameter set 2 from ODE MH
plot(Par_set_2(5),Par_set_2(8),'w+','markersize',12,'linewidth',8); hold on   
mh3=plot(Par_set_2(5),Par_set_2(8),'k+','markersize',12,'linewidth',4); hold on
% Plot of true parameter values for parameters 5 and 8
mh4=plot(log10(abs(Parameters_True(5))),log10(abs(Parameters_True(8))),'o',...   
    'markerfacecolor',[0.4,0.4,0.4],'MarkerEdgeColor',[0.3,0.3,0.3],'markersize',12);
set(gca,'position',[0.11, 0.16905, 0.33466, 0.68843])
legend([mh1,mh2,mh3,mh4],{'ODE','Parameter Set 1','Parameter Set 2','True Parameters'},...
    'FontSize',11,'Location','SouthEast')
xlabel('\it\bf log_{10}k_{r2}','FontWeight','bold','FontSize',16); 
ylabel('\bf log_{10}\gamma','FontWeight','bold','FontSize',16)
title('Parameter Uncertainty','FontWeight','bold','FontSize',12)


% Plot normalized covariances for ODE MH search for all parameter combinations 
subplot(1,2,2)
colormap gray
set(gca,'position',[0.53393, 0.16905, 0.4, 0.68843]);
pcolor(abs(cov_all_M2_I1))
set(gca,'xtick',[1.5:1:8.5],'xticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',14)
set(gca,'ytick',[1.5:1:8.5],'yticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',14)
title('Covariance Comparison: ODE','FontWeight','bold','FontSize',12)
colorbar
for i=1:8
    text(i+.45,i+.45,'+','FontSize',14)
    for j=i+0:8
        if abs(cov_all_M2_I1(i,j))>0.02
            coll = 'k';
        else
            coll = 'w';
        end
        if cov_all_M2_I1(i,j)>0
            text(i+.45,j+.45,'+','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'+','FontSize',16,'Color',coll)
        else
            text(i+.45,j+.45,'-','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'-','FontSize',16,'Color',coll)
        end       
    end
end
suptitle('ODE Metropolis-Hastings Results') % Label the set of scatter plots

%% FIGURES FOR CHAPTER: Figure 5: SSA RESULTS USING MULTIPLE PARAMETERS SETS, INPUT 1, AND MODEL 2

Nm = 2; Input = Input1; % Model and Input Signal
Output_Times = linspace(0,100,100);  % Vector with output times

figure()
% Plot single SSA trajectory using ODE best fit parameters
subplot(1,3,1) 
plot(Output_Times,m_RNA_Array_T1_I1(820,:),'Color',[0,0,0],'linewidth',1.5);
set(gca,'xlim',[Output_Times(1),Output_Times(end)],'ylim',[0 22])
title('Single SSA Trajectory', 'FontWeight','bold','FontSize',13) 
xlabel('Time (min)','FontWeight','bold','FontSize',13) 
ylabel('Average mRNA Count','FontWeight','bold','FontSize',13) 

% Plot mean mRNA from multiple runs of many SSA trajectories using ODE best
% fit parameters compared to ODE solution for mean mRNA
subplot(1,3,2)
% Plot mean mRNA from ODE solution using ODE best fit parameters
MVT_test =  odefun(Model_Results_Parameters_I1(Nm,:),Nm,Input,Output_Times,x0);
h1=plot(Output_Times,MVT_test,'Color',[0.1,0.1,0.1],'linewidth',3); hold on
% Plot the mean of multiple runs of 1000 SSA trajectories using ODE best fit parameters
h2=plot(Output_Times,mean(m_RNA_Array_T1_I1),':','Color',[0.1,0.1,0.1],'linewidth',2); 
h3=plot(Output_Times,mean(m_RNA_Array_T2_I1),':','Color',[0.3,0.3,0.3],'linewidth',2);
h4=plot(Output_Times,mean(m_RNA_Array_T3_I1),':','Color',[0.5,0.5,0.5],'linewidth',2); 
h5=plot(Output_Times,mean(m_RNA_Array_T4_I1),':','Color',[0.7,0.7,0.7],'linewidth',2);  
% Plot the mean of the data at each time point
h6=plot(Output_Times_data,mean(Data_Set),'o','Color',[0.1,0.1,0.1],'linewidth',1.5);
title('1000 SSA Trajectories','FontWeight','bold','FontSize',13) 
xlabel('Time (min)','FontWeight','bold','FontSize',13) 
ylabel('Average mRNA Count','FontWeight','bold','FontSize',13) %
set(gca,'xlim',[Output_Times(1),Output_Times(end)],'ylim',[0 14])
legend([h1,h2,h3,h4,h5,h6],{'ODE: Best Fit','SSA 1','SSA 2',...
    'SSA 3','SSA 4','Data'},'FontSize',11,'Location','SouthEast') 

% Plot data distribution at a specific time point along with FSP distributions
% at the same time point, along with FSP distribution at the same time point
% using the true parameter set
subplot(1,3,3)
% Data Distribution at a specific time point
[Hb1] = histogram(Data_Set(1:end,5),'Normalization','probability','FaceColor',[0.6,0.6,0.6]); hold on 
bins1 = Hb1.BinEdges+(0.5*Hb1.BinWidth); hvalb1 = Hb1.Values;
% FSP Distributions at the same time point as data distribution
OBJ_FSP_MH_P = @(X)get_FSP_OBJ(10.^X,Input,Nm,Output_Times,zeros(100,50),Mats_A,N,N_states,'P_RNA');
P_RNA_1 = OBJ_FSP_MH_P(Par_set_1);
fsp1 = plot([0:30],P_RNA_1(44,1:31),'Color',[0.1,0.1,0.1],'linewidth',3);
P_RNA_2 = OBJ_FSP_MH_P(Par_set_2);
fsp2 = plot([0:30],P_RNA_2(44,1:31),'Color',[0.4,0.4,0.4],'linewidth',3);
P_RNA_t = OBJ_FSP_MH_P(log10(abs(Parameters_True)));
fsp_t = plot([0:30],P_RNA_t(44,1:31),'k:','linewidth',3);
set(gca,'ylim',[0 0.15],'xlim',[0 30]);
legend([Hb1,fsp1,fsp2,fsp_t],{'Data','FSP: Par Set 1','FSP: Par Set 2','FSP: True Params'},...
    'FontSize',10,'Location','NorthEast') 
xlabel('mRNA Level','FontWeight','bold','FontSize',13) 
ylabel('Probability','FontWeight','bold','FontSize',13)
title('Distributions (t=44 min)','FontWeight','bold','FontSize',13) 
suptitle('SSA Results') % Label the set of scatter plots

%% FIGURES FOR CHAPTER: Figure 6: DATA DISTRIBUTIONS COMPARED TO FSP DISTRIBUTIONS
% Plot data distributions compared to full fsp distributions using FSP best fit
% parameters from model 2 and input signal 1 (at 5 time points)

figure()
subplot(1,5,1) 
histogram(Data_Set(1:end,2),'Normalization','probability','FaceColor',[0.5,0.5,0.5]); hold on
plot([0:40], P_FSP_M2_I1(2,1:41),'k','linewidth',3);
set(gca,'xlim',[0,40],'ylim',[0,0.15],'FontSize',12) 
title('t = 11 min','FontWeight','bold','FontSize',13) 
ylabel('Probability','FontWeight','bold','FontSize',14) 
legend({'Data','FSP'},'FontSize',13,'Location','North') 

subplot(1,5,2) 
histogram(Data_Set(1:end,4),'Normalization','probability','FaceColor',[0.5,0.5,0.5]); hold on
plot([0:40], P_FSP_M2_I1(4,1:41),'k','linewidth',3);
set(gca,'xlim',[0,40],'ylim',[0,0.15],'FontSize',12)
title('t = 33 min','FontWeight','bold','FontSize',13)
legend({'Data','FSP'},'FontSize',13,'Location','North')

subplot(1,5,3)  
histogram(Data_Set(1:end,6),'Normalization','probability','FaceColor',[0.5,0.5,0.5]); hold on
plot([0:40], P_FSP_M2_I1(6,1:41),'k','linewidth',3);
set(gca,'xlim',[0,40],'ylim',[0,0.15],'FontSize',12)
title('t = 56 min','FontWeight','bold','FontSize',13)
legend({'Data','FSP'},'FontSize',13,'Location','North') 

subplot(1,5,4)  
histogram(Data_Set(1:end,8),'Normalization','probability','FaceColor',[0.5,0.5,0.5]);hold on
plot([0:40], P_FSP_M2_I1(8,1:41),'k','linewidth',3);
set(gca,'xlim',[0,40],'ylim',[0,0.15],'FontSize',12)
title('t = 78 min','FontWeight','bold','FontSize',13) 
legend({'Data','FSP'},'FontSize',13,'Location','North') 

subplot(1,5,5)  
histogram(Data_Set(1:end,10),'Normalization','probability','FaceColor',[0.5,0.5,0.5]); hold on
plot([0:40], P_FSP_M2_I1(10,1:41),'k','linewidth',3);
set(gca,'xlim',[0,40],'ylim',[0,0.15],'FontSize',12)
title('t = 100 min','FontWeight','bold','FontSize',13) 
legend({'Data','FSP'},'FontSize',13,'Location','North') 

[ax1,h1]=suplabel('mRNA Level','x',[0.11 0.12 0.82 0.82]);  % One x label for all scatter plots
set(h1,'FontWeight','bold','FontSize',14)
% suptitle('Data Distributions at Different Times') % Label the set of scatter plots

%% FIGURES FOR CHAPTER: Figure 7: FSP MET-HASTE RESULTS FOR INPUT 1 AND MODEL 2
% Excluding burn-in period (estimated to be 1/4 of total chain length)

% Plot ODE and FSP MH scatter plots with 90% confidence ellipses
figure()
subplot(1,3,1)
% Make a scatter plot of Par5 vs. Par8 for ODE and/or FSP MH searches
s2=scatter(I1_FSP_M2_Par5(floor(end/4):end),I1_FSP_M2_Par8(floor(end/4):end),10,'MarkerEdgeColor',[0.5,0.5,0.5]);hold on 
% Plot parameter set 1 from ODE MH
plot(Par_set_1(5),Par_set_1(8),'wx','markersize',12,'linewidth',8); hold on
s3=plot(Par_set_1(5),Par_set_1(8),'kx','markersize',12,'linewidth',4); hold on
% Plot parameter set 2 from ODE MH
plot(Par_set_2(5),Par_set_2(8),'w+','markersize',12,'linewidth',8); hold on
s4=plot(Par_set_2(5),Par_set_2(8),'k+','markersize',12,'linewidth',4); hold on
% Plot of true parameter values for parameters 5 and 8
s5=plot(log10(abs(Parameters_True(5))),log10(abs(Parameters_True(8))),'o',...
    'markerfacecolor',[0.1,0.1,0.1],'MarkerEdgeColor',[0.1,0.1,0.1],'markersize',12);
COL = 'k'; conf = 0.90; % For FSP ellipses: color and confidence interval
[~,XX_FSP,YY_FSP] = error_ellipse_fill(RSLTS_FSP_M2_I1(5,8).COV,RSLTS_FSP_M2_I1(5,8).MN,'style',COL,'conf',conf);  % Call routine to plots ellipses for the 90% confidense intervals.
xlabel('\it\bf  log_{10}k_{r2}','FontWeight','bold','FontSize',16); 
ylabel('\bf log_{10}\gamma','FontWeight','bold','FontSize',16)
title('Parameter Uncertainty','FontWeight','bold','FontSize',12)
legend([s2,s3,s4,s5],{'FSP','Parameter Set 1','Parameter Set 2','True Parameters'},...
    'FontSize',11,'Location','SouthEast')

% Plot normalized covariances for ODE MH search for all parameter combinations 
subplot(1,3,2)
colormap gray
pcolor(abs(cov_all_M2_I1))
set(gca,'xtick',[1.5:1:8.5],'xticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',14)
set(gca,'ytick',[1.5:1:8.5],'yticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',14)
title('Covariance Comparison: ODE','FontWeight','bold','FontSize',12)
for i=1:8
    text(i+.45,i+.45,'+','FontSize',14)
    for j=i+0:8
        if abs(cov_all_M2_I1(i,j))>0.02
            coll = 'k';
        else
            coll = 'w';
        end
        if cov_all_M2_I1(i,j)>0
            text(i+.45,j+.45,'+','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'+','FontSize',16,'Color',coll)
        else
            text(i+.45,j+.45,'-','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'-','FontSize',16,'Color',coll)
        end       
    end
end

% Plot normalized covariances for FSP MH search for all parameter combinations 
subplot(1,3,3)
colormap gray
pcolor(abs(cov_all_fsp_M2_I1))
set(gca,'xtick',[1.5:1:8.5],'xticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',14)
set(gca,'ytick',[1.5:1:8.5],'yticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',14)
title('Covariance Comparison: FSP','FontWeight','bold','FontSize',12)
colorbar
for i=1:8
    text(i+.45,i+.45,'+','FontSize',14)
    for j=i+0:8
        if abs(cov_all_fsp_M2_I1(i,j))>0.02
            coll = 'k';
        else
            coll = 'w';
        end
        if cov_all_fsp_M2_I1(i,j)>0
            text(i+.45,j+.45,'+','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'+','FontSize',16,'Color',coll)
        else
            text(i+.45,j+.45,'-','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'-','FontSize',16,'Color',coll)
        end       
    end
end
suptitle('Comparison between ODE and FSP Metropolis-Hastings Results') % Label the set of scatter plots

%% FIGURES FOR CHAPTER: Figure 8: BAR PLOT COMPARING PARAMETER RESULTS FROM ODE AND FSP MET-HASTE SEARCHES
% Excluding burn-in period (estimated to be 1/4 of total chain length)

Nm = 2; %Specify Model Number
% Bar plot of true parameter values
TRUEPAR = abs([Parameters_True]);
% Bar plot of mean FSP parameter values (from MH search)
FSPBAR = 10.^mean(Rslts_FSP_I1(Nm).Par_Chain(floor(end/4):end,:));
% Bar plot of mean ODE parameter values (from MH search)
ODEBAR = 10.^mean(Rslts_ODE_I1(Nm).Par_Chain(floor(end/4):end,:));
% Standard deviation of FSP parameter values (from MH search)
std_FSP = std(10.^Rslts_FSP_I1(Nm).Par_Chain(floor(end/4):end,:));
% Standard deviation of ODE parameter values (from MH search)
std_ODE = std(10.^Rslts_ODE_I1(Nm).Par_Chain(floor(end/4):end,:));
BARS = [TRUEPAR;ODEBAR;FSPBAR]';
ERRORS = [zeros(size(TRUEPAR(1,:)));std_ODE;std_FSP]';

fq = figure;
b=bar([TRUEPAR;ODEBAR;FSPBAR]','BarWidth',1); hold on
b(1).FaceColor = [0.2,0.2,0.2]; b(2).FaceColor = [0.6,0.6,0.6]; b(3).FaceColor = [0.9,0.9,0.9];
legend({'True Parameter Values','ODE','FSP'},'Location','NorthEast','FontSize',14) 
set(gca,'yscale','log','xtick',[1:1:8.5],'xticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',14)
ylabel('Parameter Mean Value', 'FontWeight', 'bold','FontSize',14)
title('Comparison of Mean Parameter Values From Metropolis-Hastings Search','FontWeight','bold','FontSize',14)

% Alignment of error bars with corresponding individual bar on plot
numgroups = size(BARS,1); 
numbars = size(BARS,2); 
groupwidth = min(0.8, numbars/(numbars+1.5));
for i = 1:numbars
      x = (1:numgroups)- groupwidth/2 + (2*i-1) * groupwidth / (2*numbars);
      errorbar(x,BARS(:,i),min((BARS(:,i)-1e-6),ERRORS(:,i)),ERRORS(:,i), 'k','linestyle','none','linewidth',0.8);
end
set(gca,'ylim',[1e-3,1e5])

%% FIGURES FOR CHAPTER: Figure 9: ODE MET-HASTE RESULTS - ALL THREE INPUT SIGNALS  (Sinusoidal, Step, Ramp)

% Plot MH Results with ODE Objective Function (chi-squared function) - All three input signals
figure()
subplot(2,4,[1 5])
% Make scatter plots of Par5 vs. Par8
mh1_I1=scatter(I1_ODE_M2_Par5(floor(end/4):end),I1_ODE_M2_Par8(floor(end/4):end),10,'Marker','none'); hold on % Input 1: Sinusoidal
mh1_I2=scatter(I2_ODE_M2_Par5(floor(end/4):end),I2_ODE_M2_Par8(floor(end/4):end),10,'Marker','none'); hold on % Input 2: Step
mh1_I3=scatter(I3_ODE_M2_Par5(floor(end/4):end),I3_ODE_M2_Par8(floor(end/4):end),10,'Marker','none'); hold on % Input 3: Ramp
% Call routine to plot ellipses for the 90% confidense intervals for all three inputs 
COL = [0.2,0.2,0.2]; conf = 0.90; % For Input 1 ellipses: color and confidence interval
[h1,XX,YY] = error_ellipse_fill(RSLTS_ODE_M2_I1(5,8).COV,RSLTS_ODE_M2_I1(5,8).MN,'style',COL,'conf',conf); 
COL = [0.5,0.5,0.5]; conf = 0.90; % For Input 2 ellipses: color and confidence interval
[h2,XX_I2,YY_I2] = error_ellipse_fill(RSLTS_ODE_M2_I2(5,8).COV,RSLTS_ODE_M2_I2(5,8).MN,'style',COL,'conf',conf); 
COL = [0.8,0.8,0.8]; conf = 0.90; % For Input 3 ellipses: color and confidence interval
[h3,XX_I3,YY_I3] = error_ellipse_fill(RSLTS_ODE_M2_I3(5,8).COV,RSLTS_ODE_M2_I3(5,8).MN,'style',COL,'conf',conf); 
% Plot of true parameter values for parameters 5 and 8
mh2=plot(log10(abs(Parameters_True(5))),log10(abs(Parameters_True(8))),'o',...
    'markerfacecolor','k','MarkerEdgeColor','k','markersize',12);
set(gca,'xlim',[-1.0 4.5],'ylim',[-5 2.0])
xlabel('\it\bf  log_{10}k_{r2}','FontWeight','bold','FontSize',14); 
ylabel('\bf log_{10}\gamma','FontWeight','bold','FontSize',14)
title('Parameter Uncertainty','FontWeight','bold','FontSize',12)
legend([h1,h2,h3,mh2],{'ODE:Input 1','ODE:Input 2','ODE:Input 3','True Parameters'},...
    'FontSize',10,'Location','SouthEast')

% For plotting input signals
time = [0 100]; xval = Output_Times_data; yval = [0 3];

% Plot Input Signal 1: Sinusoidal
subplot(2,4,2)
for i=1:length(xval)
    line([xval(i),xval(i)],[yval(1),yval(2)],'LineStyle',':','Color',[0.5,0.5,0.5]);hold on
end
fI1 = ezplot(Input1,time);hold on
set(gca,'ylim',[0,3]); set(fI1,'Color','k','LineWidth',2)
xlabel('Time (min)','FontWeight','bold','FontSize',12); 
ylabel('Y_{1}(t)','FontWeight','bold','FontSize',12)
title('Input Signal 1: Sinusoidal','FontWeight','bold','FontSize',12)

% Plot Input Signal 2: Step
subplot(2,4,3)
for i=1:length(xval)
    line([xval(i),xval(i)],[yval(1),yval(2)],'LineStyle',':','Color',[0.5,0.5,0.5]);hold on
end
fI2 = ezplot(Input2,time);hold on
set(gca,'ylim',[0,3]); set(fI2,'Color','k','LineWidth',2)
xlabel('Time (min)','FontWeight','bold','FontSize',12); 
ylabel('Y_{2}(t)','FontWeight','bold','FontSize',12)
title('Input Signal 2: Step','FontWeight','bold','FontSize',12)

% Plot Input Signal 3: Ramp
subplot(2,4,4)
for i=1:length(xval)
    line([xval(i),xval(i)],[yval(1),yval(2)],'LineStyle',':','Color',[0.5,0.5,0.5]);hold on
end
fI3 = ezplot(Input3,time);hold on
set(gca,'ylim',[0,3]); set(fI3,'Color','k','LineWidth',2)
xlabel('Time (min)','FontWeight','bold','FontSize',12); 
ylabel('Y_{3}(t)','FontWeight','bold','FontSize',12)
title('Input Signal 3: Ramp','FontWeight','bold','FontSize',12)

% Plot normalized covariances for ODE MH search for all parameter combinations: Input 1
subplot(2,4,6)
colormap gray
pcolor(abs(cov_all_M2_I1))
set(gca,'xtick',[1.5:1:8.5],'xticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',12)
set(gca,'ytick',[1.5:1:8.5],'yticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',12)
title('Covariances: Input 1','FontWeight','bold','FontSize',12)
for i=1:8
    text(i+.45,i+.45,'+','FontSize',14)
    for j=i+0:8
        if abs(cov_all_M2_I1(i,j))>0.02
            coll = 'k';
        else
            coll = 'w';
        end
        if cov_all_M2_I1(i,j)>0
            text(i+.45,j+.45,'+','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'+','FontSize',16,'Color',coll)
        else
            text(i+.45,j+.45,'-','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'-','FontSize',16,'Color',coll)
        end       
    end
end

% Plot normalized covariances for ODE MH search for all parameter combinations: Input 2
subplot(2,4,7)
colormap gray
pcolor(abs(cov_all_M2_I2))
set(gca,'xtick',[1.5:1:8.5],'xticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',12)
set(gca,'ytick',[1.5:1:8.5],'yticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',12)
title('Covariances: Input 2','FontWeight','bold','FontSize',12)
for i=1:8
    text(i+.45,i+.45,'+','FontSize',14)
    for j=i+0:8
        if abs(cov_all_M2_I2(i,j))>0.02
            coll = 'k';
        else
            coll = 'w';
        end
        if cov_all_M2_I2(i,j)>0
            text(i+.45,j+.45,'+','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'+','FontSize',16,'Color',coll)
        else
            text(i+.45,j+.45,'-','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'-','FontSize',16,'Color',coll)
        end       
    end
end

% Plot normalized covariances for ODE MH search for all parameter combinations: Input 3
subplot(2,4,8)
colormap gray
pcolor(abs(cov_all_M2_I3))
set(gca,'xtick',[1.5:1:8.5],'xticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',12)
set(gca,'ytick',[1.5:1:8.5],'yticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',12)
title('Covariances: Input 3','FontWeight','bold','FontSize',12)
colorbar
for i=1:8
    text(i+.45,i+.45,'+','FontSize',14)
    for j=i+0:8
        if abs(cov_all_M2_I3(i,j))>0.02
            coll = 'k';
        else
            coll = 'w';
        end
        if cov_all_M2_I3(i,j)>0
            text(i+.45,j+.45,'+','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'+','FontSize',16,'Color',coll)
        else
            text(i+.45,j+.45,'-','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'-','FontSize',16,'Color',coll)
        end       
    end
end
suptitle('ODE Metropolis-Hastings Results: 3 Different Inputs') % Label the set of scatter plots

%% FIGURES FOR CHAPTER: Figure 10: FSP MET-HASTE RESULTS - ALL THREE INPUT SIGNALS (Sinusoidal, Step, Ramp)

% Plot MH Results with FSP Objective Function - All three input signals
figure()   
% Make scatter plots of Par5 vs. Par8
subplot(2,4,[1 5])
s2_I1=scatter(I1_FSP_M2_Par5(floor(end/4):end),I1_FSP_M2_Par8(floor(end/4):end),10,'Marker','none');hold on % Input 1: Sinusoidal
s2_I2=scatter(I2_FSP_M2_Par5(floor(end/4):end),I2_FSP_M2_Par8(floor(end/4):end),10,'Marker','none');hold on % Input 2: Step
s2_I3=scatter(I3_FSP_M2_Par5(floor(end/4):end),I3_FSP_M2_Par8(floor(end/4):end),10,'Marker','none');hold on % Input 3: Ramp
% Call routine to plot ellipses for the 90% confidense intervals for all three inputs 
COL = [0.2,0.2,0.2]; conf = 0.90; % For Input 1 ellipses: color and confidence interval
[h1_fsp,XX_FSP,YY_FSP] = error_ellipse_fill(RSLTS_FSP_M2_I1(5,8).COV,RSLTS_FSP_M2_I1(5,8).MN,'style',COL,'conf',conf);
COL = [0.5,0.5,0.5]; conf = 0.90; % For Input 2 ellipses: color and confidence interval
[h2_fsp,XX_FSP_I2,YY_FSP_I2] = error_ellipse_fill(RSLTS_FSP_M2_I2(5,8).COV,RSLTS_FSP_M2_I2(5,8).MN,'style',COL,'conf',conf); 
COL = [0.8,0.8,0.8]; conf = 0.90; % For Input 3 ellipses: color and confidence interval
[h3_fsp,XX_FSP_I3,YY_FSP_I3] = error_ellipse_fill(RSLTS_FSP_M2_I3(5,8).COV,RSLTS_FSP_M2_I3(5,8).MN,'style',COL,'conf',conf); 
% Plot of true parameter values for parameters 5 and 8
s3=plot(log10(abs(Parameters_True(5))),log10(abs(Parameters_True(8))),'o',...
    'markerfacecolor','k','MarkerEdgeColor','k','markersize',12);
set(gca,'xlim',[0.5 3.0],'ylim',[-0.5 1.7])
xlabel('\it\bf log_{10}k_{r2}','FontWeight','bold','FontSize',14); 
ylabel('\bf log_{10}\gamma','FontWeight','bold','FontSize',14)
title('Parameter Uncertainty','FontWeight','bold','FontSize',12)
legend([h1_fsp,h2_fsp,h3_fsp,s3],{'FSP:Input1','FSP:Input2','FSP:Input3','True Parameters'},...
    'FontSize',10,'Location','SouthEast')

% For plotting input signals
time = [0 100]; xval = Output_Times_data; yval = [0 3];

% Plot Input Signal 1: Sinusoidal
subplot(2,4,2)
for i=1:length(xval)
    line([xval(i),xval(i)],[yval(1),yval(2)],'LineStyle',':','Color',[0.5,0.5,0.5]);hold on
end
fI1 = ezplot(Input1,time);hold on
set(gca,'ylim',[0,3]); set(fI1,'Color','k','LineWidth',2)
xlabel('Time (min)','FontWeight','bold','FontSize',12); 
ylabel('Y_{1}(t)','FontWeight','bold','FontSize',12)
title('Input Signal 1: Sinusoidal','FontWeight','bold','FontSize',12)

% Plot Input Signal 2: Step
subplot(2,4,3)
for i=1:length(xval)
    line([xval(i),xval(i)],[yval(1),yval(2)],'LineStyle',':','Color',[0.5,0.5,0.5]);hold on
end
fI2 = ezplot(Input2,time);hold on
set(gca,'ylim',[0,3]); set(fI2,'Color','k','LineWidth',2)
xlabel('Time (min)','FontWeight','bold','FontSize',12); 
ylabel('Y_{2}(t)','FontWeight','bold','FontSize',12)
title('Input Signal 2: Step','FontWeight','bold','FontSize',12)

% Plot Input Signal 3: Ramp
subplot(2,4,4)
for i=1:length(xval)
    line([xval(i),xval(i)],[yval(1),yval(2)],'LineStyle',':','Color',[0.5,0.5,0.5]);hold on
end
fI3 = ezplot(Input3,time);hold on
set(gca,'ylim',[0,3]); set(fI3,'Color','k','LineWidth',2)
xlabel('Time (min)','FontWeight','bold','FontSize',12); 
ylabel('Y_{3}(t)','FontWeight','bold','FontSize',12)
title('Input Signal 3: Ramp','FontWeight','bold','FontSize',12)

% Plot normalized covariances for FSP MH search for all parameter combinations: Input 1
subplot(2,4,6)
colormap gray
pcolor(abs(cov_all_fsp_M2_I1))
set(gca,'xtick',[1.5:1:8.5],'xticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',12)
set(gca,'ytick',[1.5:1:8.5],'yticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',12)
title('FSP Covariances: Input 1','FontWeight','bold','FontSize',12)
for i=1:8
    text(i+.45,i+.45,'+','FontSize',14)
    for j=i+0:8
        if abs(cov_all_fsp_M2_I1(i,j))>0.02
            coll = 'k';
        else
            coll = 'w';
        end
        if cov_all_fsp_M2_I1(i,j)>0
            text(i+.45,j+.45,'+','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'+','FontSize',16,'Color',coll)
        else
            text(i+.45,j+.45,'-','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'-','FontSize',16,'Color',coll)
        end       
    end
end

% Plot normalized covariances for FSP MH search for all parameter combinations: Input 2
subplot(2,4,7)
colormap gray
pcolor(abs(cov_all_fsp_M2_I2))
set(gca,'xtick',[1.5:1:8.5],'xticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',12)
set(gca,'ytick',[1.5:1:8.5],'yticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',12)
title('FSP Covariances: Input 2','FontWeight','bold','FontSize',12)
for i=1:8
    text(i+.45,i+.45,'+','FontSize',14)
    for j=i+0:8
        if abs(cov_all_fsp_M2_I2(i,j))>0.02
            coll = 'k';
        else
            coll = 'w';
        end
        if cov_all_fsp_M2_I2(i,j)>0
            text(i+.45,j+.45,'+','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'+','FontSize',16,'Color',coll)
        else
            text(i+.45,j+.45,'-','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'-','FontSize',16,'Color',coll)
        end       
    end
end

% Plot normalized covariances for FSP MH search for all parameter combinations: Input 3
subplot(2,4,8)
colormap gray
pcolor(abs(cov_all_fsp_M2_I3))
set(gca,'xtick',[1.5:1:8.5],'xticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',12)
set(gca,'ytick',[1.5:1:8.5],'yticklabel',{'\it\bf k_{12}','\it\bf k_{23}',...
    '\it\bf k_{21}','\it\bf k_{32}','\it\bf k_{r2}','\it\bf k_{r3}','\bf\beta','\bf\gamma'},...
    'FontSize',12)
title('FSP Covariances: Input 3','FontWeight','bold','FontSize',12)
colorbar
for i=1:8
    text(i+.45,i+.45,'+','FontSize',14)
    for j=i+0:8
        if abs(cov_all_fsp_M2_I3(i,j))>0.02
            coll = 'k';
        else
            coll = 'w';
        end
        if cov_all_fsp_M2_I3(i,j)>0
            text(i+.45,j+.45,'+','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'+','FontSize',16,'Color',coll)
        else
            text(i+.45,j+.45,'-','FontSize',16,'Color',coll)
            text(j+.45,i+.45,'-','FontSize',16,'Color',coll)
        end       
    end
end
suptitle('FSP Metropolis-Hastings Results: 3 Different Inputs') % Label the set of scatter plots

%% FIGURES FOR CHAPTER: Figure 11: Predictions - Multiple Time / Input Comparisons
% Plot predictions for multiple different experiments (different input
% signals) at different times using the best parameters from the FSP
% Metropolis-Hastings search of parameter space using input 1 

% For plotting input signals
time = [0,100]; xval = Output_Times_expt; yval = [0 3]; 

figure()
% Plot Input Signal 2: Step
subplot(3,3,1)
for i=1:length(xval)
    line([xval(i),xval(i)],[yval(1),yval(2)],'LineStyle',':','Color',[0.5,0.5,0.5]);hold on
end
fI2 = ezplot(Input_Exp_1,time);hold on
set(gca,'ylim',[0,3]); set(fI2,'Color','k','LineWidth',2)
xlabel('Time (min)','FontWeight','bold','FontSize',11); 
ylabel('Y_{2}(t)','FontWeight','bold','FontSize',11)
title('Input Signal 2: Step','FontWeight','bold','FontSize',12)

% Plot Input Signal 3: Ramp
subplot(3,3,2)
for i=1:length(xval)
    line([xval(i),xval(i)],[yval(1),yval(2)],'LineStyle',':','Color',[0.5,0.5,0.5]);hold on
end
fI3 = ezplot(Input_Exp_2,time);hold on
set(gca,'ylim',[0,3]); set(fI3,'Color','k','LineWidth',2) 
xlabel('Time (min)','FontWeight','bold','FontSize',11); 
ylabel('Y_{3}(t)','FontWeight','bold','FontSize',11)
title('Input Signal 3: Ramp','FontWeight','bold','FontSize',12)

% Plot Input Signal 4: Sinusoidal (Fast Frequency)
subplot(3,3,3)
for i=1:length(xval)
    line([xval(i),xval(i)],[yval(1),yval(2)],'LineStyle',':','Color',[0.5,0.5,0.5]);hold on
end
fI4 = ezplot(Input_Exp_3,time);hold on
set(gca,'ylim',[0,3]); set(fI4,'Color','k','LineWidth',2) 
xlabel('Time (min)','FontWeight','bold','FontSize',11); 
ylabel('Y_{4}(t)','FontWeight','bold','FontSize',11)
title('Input Signal: Sinusoidal (Fast)','FontWeight','bold','FontSize',12)

% Make adjustments to FSP plots here
COL_TRUE = 'k'; COL1 = [0.2,0.2,0.2]; COL2 = [0.4,0.4,0.4]; COL3 = [0.6,0.6,0.6]; COL4 = [0.8,0.8,0.8]; 
STYLE_TRUE = '-'; STYLE1 = '-'; STYLE2 = '-'; STYLE3 = '-'; STYLE4 = '-';
LOC = 'NorthEast'; XLIM = 40; LIN = 1.5; LEG_SZ = 11;

% Compare predictions for models 1-4 using FSP best fit parameters from MH search (with input signal 1)
% Predictions using input signal 2 (step)
subplot(3,3,4)
plot([0:XLIM], P_RNA_EXPT_M1_I1_1(4,1:(XLIM+1)),'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1); hold on
plot([0:XLIM], P_RNA_EXPT_M2_I1_1(4,1:(XLIM+1)),'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2);
plot([0:XLIM], P_RNA_EXPT_M3_I1_1(4,1:(XLIM+1)),'Color',COL3,'linewidth',LIN,'LineStyle',STYLE3);
plot([0:XLIM], P_RNA_EXPT_M4_I1_1(4,1:(XLIM+1)),'Color',COL4,'linewidth',LIN,'LineStyle',STYLE4);
set(gca,'xlim',[0,XLIM],'ylim',[0,0.15])
title(sprintf('t = %i min',round(Output_Times_expt(4))),'FontWeight','bold','FontSize',12)
legend({'Model 1','Model 2','Model 3','Model 4'},'FontSize',LEG_SZ,'Location',LOC)
 
% Predictions using input signal 3 (ramp)
subplot(3,3,5)
plot([0:XLIM], P_RNA_EXPT_M1_I1_2(4,1:(XLIM+1)),'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1); hold on
plot([0:XLIM], P_RNA_EXPT_M2_I1_2(4,1:(XLIM+1)),'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2);
plot([0:XLIM], P_RNA_EXPT_M3_I1_2(4,1:(XLIM+1)),'Color',COL3,'linewidth',LIN,'LineStyle',STYLE3);
plot([0:XLIM], P_RNA_EXPT_M4_I1_2(4,1:(XLIM+1)),'Color',COL4,'linewidth',LIN,'LineStyle',STYLE4);
set(gca,'xlim',[0,XLIM],'ylim',[0,0.15])
title(sprintf('t = %i min',round(Output_Times_expt(4))),'FontWeight','bold','FontSize',12)
legend({'Model 1','Model 2','Model 3','Model 4'},'FontSize',LEG_SZ,'Location',LOC)

% Predictions using sinusoidal input signal (with fast frequency)
subplot(3,3,6)
plot([0:XLIM], P_RNA_EXPT_M1_I1_3(4,1:(XLIM+1)),'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1); hold on
plot([0:XLIM], P_RNA_EXPT_M2_I1_3(4,1:(XLIM+1)),'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2);
plot([0:XLIM], P_RNA_EXPT_M3_I1_3(4,1:(XLIM+1)),'Color',COL3,'linewidth',LIN,'LineStyle',STYLE3);
plot([0:XLIM], P_RNA_EXPT_M4_I1_3(4,1:(XLIM+1)),'Color',COL4,'linewidth',LIN,'LineStyle',STYLE4);
set(gca,'xlim',[0,XLIM],'ylim',[0,0.15])
title(sprintf('t = %i min',round(Output_Times_expt(4))),'FontWeight','bold','FontSize',12)
legend({'Model 1','Model 2','Model 3','Model 4'},'FontSize',LEG_SZ,'Location',LOC)

% Predictions using input signal 2 (step)
subplot(3,3,7)
plot([0:XLIM], P_RNA_EXPT_M1_I1_1(11,1:(XLIM+1)),'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1); hold on
plot([0:XLIM], P_RNA_EXPT_M2_I1_1(11,1:(XLIM+1)),'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2);
plot([0:XLIM], P_RNA_EXPT_M3_I1_1(11,1:(XLIM+1)),'Color',COL3,'linewidth',LIN,'LineStyle',STYLE3);
plot([0:XLIM], P_RNA_EXPT_M4_I1_1(11,1:(XLIM+1)),'Color',COL4,'linewidth',LIN,'LineStyle',STYLE4);
set(gca,'xlim',[0,XLIM],'ylim',[0,0.15])
title(sprintf('t = %i min',round(Output_Times_expt(11))),'FontWeight','bold','FontSize',12)
legend({'Model 1','Model 2','Model 3','Model 4'},'FontSize',LEG_SZ,'Location',LOC)

% Predictions using input signal 3 (ramp)
subplot(3,3,8)
plot([0:XLIM], P_RNA_EXPT_M1_I1_2(11,1:(XLIM+1)),'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1); hold on
plot([0:XLIM], P_RNA_EXPT_M2_I1_2(11,1:(XLIM+1)),'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2);
plot([0:XLIM], P_RNA_EXPT_M3_I1_2(11,1:(XLIM+1)),'Color',COL3,'linewidth',LIN,'LineStyle',STYLE3);
plot([0:XLIM], P_RNA_EXPT_M4_I1_2(11,1:(XLIM+1)),'Color',COL4,'linewidth',LIN,'LineStyle',STYLE4);
set(gca,'xlim',[0,XLIM],'ylim',[0,0.15])
title(sprintf('t = %i min',round(Output_Times_expt(11))),'FontWeight','bold','FontSize',12)
legend({'Model 1','Model 2','Model 3','Model 4'},'FontSize',LEG_SZ,'Location',LOC)

% Predictions using sinusoidal input signal (with fast frequency)
subplot(3,3,9)
plot([0:XLIM], P_RNA_EXPT_M1_I1_3(11,1:(XLIM+1)),'Color',COL1,'linewidth',LIN,'LineStyle',STYLE1); hold on
plot([0:XLIM], P_RNA_EXPT_M2_I1_3(11,1:(XLIM+1)),'Color',COL2,'linewidth',LIN,'LineStyle',STYLE2);
plot([0:XLIM], P_RNA_EXPT_M3_I1_3(11,1:(XLIM+1)),'Color',COL3,'linewidth',LIN,'LineStyle',STYLE3);
plot([0:XLIM], P_RNA_EXPT_M4_I1_3(11,1:(XLIM+1)),'Color',COL4,'linewidth',LIN,'LineStyle',STYLE4);
set(gca,'xlim',[0,XLIM],'ylim',[0,0.15])
title(sprintf('t = %i min',round(Output_Times_expt(11))),'FontWeight','bold','FontSize',12)
legend({'Model 1','Model 2','Model 3','Model 4'},'FontSize',LEG_SZ,'Location',LOC)

% Predictions using sinusoidal input signal (with slow frequency)
subplot(3,3,4) 
ylabel('Probability','FontWeight','bold','FontSize',14) %Adding label to y axis
subplot(3,3,7) 
ylabel('Probability','FontWeight','bold','FontSize',14) %Adding label to y axis
[ax1,h1]=suplabel('mRNA Level','x',[0.09 0.09 0.84 0.84]); set(h1,'FontWeight','bold','FontSize',14) 
suptitle({'Distribution Predictions Using FSP Best Fit Parameters (from Input 1)'}) %Label the set of scatter plots

