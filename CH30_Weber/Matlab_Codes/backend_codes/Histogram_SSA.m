%This function runs many SSA trajectories and collects the results for
%the mRNA and not the entire X_Array
function [m_RNA_Array] = Histogram_SSA(Num_Runs, ModFun, x0, Output_Times)

%'Num_Runs' is number of SSA trajectories
[~,~,~,S] = ModFun(Output_Times(1),x0); %Get stoichiometry from model function
N_Species = size(S,1);  % Number of species
N_Rxns = size(S,2);  % Number of reactions
Nt = size(Output_Times,2); %Number of time points
X_Array = zeros(Num_Runs,N_Species,Nt); %Preallocate for simulation results to populate 
                                 
%Runs the SSA for 'Num_Runs' times, so that X_Array is now [Num_Runs,N_Species,Output_Times] 
%instead of [N_Species,Output_Times] like in single SSA trajectory run
for i=1:Num_Runs
    [X_Array(i,:,:)] = SSA(ModFun,x0,Output_Times);
end

%Record the mRNA results from the simulation and squeezes down from
%3D to 2D since we are only recording results for 1 species
m_RNA_Array(:,:) = squeeze(X_Array(:,4,:));
end
