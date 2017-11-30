% This function defines the stoichiometry of the system and determines the
% propensity functions
function [W0,W1,P,S] = Model_SSA(t,x,Parameters,Input,Nm)

% Runs the parameter function below, which adjusts parameters to include
% kinase-dependent rate based on model number
[New_Params] = Dispay_Params(t,Parameters,Input,Nm);

% Defines parameters in the New_Params vector: [k12 k23 k21 k32 kr2 kr3 b g];
k12 = New_Params(1); k23 = New_Params(2);
k21 = New_Params(3); k32 = New_Params(4);
kr2 = New_Params(5); kr3 = New_Params(6); 
b = New_Params(7); g = New_Params(8);
     
% Stoichiometry Matrix - Please note this is transposed
S = [0 0 0 0;  % Initial fast rate w/ no state change
    -1 1 0 0;  % Rxn 1: S1 -> S2
    1 -1 0 0;  % Rxn 2: S2 -> S1
    0 -1 1 0;  % Rxn 3: S2 -> S3
    0 1 -1 0;  % Rxn 4: S3 -> S2
    0 0 0 1;  % Rxn 5: S2 -> mRNA
    0 0 0 1;  % Rxn 6: S3 -> mRNA
    0 0 0 -1]'; % Rxn 7: mRNA -> phi 

% Propensity functions, where Propensity  = W0 + W1*X
W0 = [100 0 0 0 0 0 0 0]'; % Sets initial fast rate of 100/min
W1 = [0 0 0 0;...  % The initial fast rate doesn't change state of system
      k12 0 0 0;... % Rxn 1: S1 -> S2
      0 k21 0 0;... % Rxn 2: S2 -> S1
      0 k23 0 0;... % Rxn 3: S2 -> S3
      0 0 k32 0;... % Rxn 4: S3 -> S2
      0 kr2 0 0;... % Rxn 5: S2 -> mRNA
      0 0 kr3 0;... % Rxn 6: S3 -> mRNA
      0 0 0 g];  % Rxn 7: mRNA -> phi
P = W0 + W1*x; % Determines propensity functions for each reaction
end


% This function will adjust the parameters to take into consideration the
% kinase-dependent rate based on the model number, which are then put into a
% new vector called New_Params that's used in the propensity calcs
function [New_Params] = Dispay_Params(t,Parameters,Input,Nm)

% Define the variables in the parameters vector
     k12 = Parameters(1); k23 = Parameters(2);
     k21 = Parameters(3); k32 = Parameters(4);
     kr2 = Parameters(5); kr3 = Parameters(6);
     b = Parameters(7); g = Parameters(8);
     
% Modify whichever k is kinase-dependent based on model number
     switch Nm
        case 1
            k12=max(0,-k12+b*Input(t));
        case 2
            k23=max(0,-k23+b*Input(t));
        case 3
            k21=max(0,k21-b*Input(t));
        case 4
            k32=max(0,k32-b*Input(t));
     end
     
% Defines the parameters with adjustments based on model number and puts them 
% into a new vector called New_Params, which is used to calculate propensities
New_Params = [k12 k23 k21 k32 kr2 kr3 b g];
end