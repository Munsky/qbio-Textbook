%Run_MH_FSP: : Runs Metropolis-Hastings with the FSP Objective Function in Logspace

disp(' ')
disp('Running the FSP MH search on the gene expression models.')
disp('This could take as much as a few hours.')
disp('You may wish to use shorter MH chains or fewer skips during the debug stage.')
disp('Then once it is working, you can leave it to run while you do something else.')

N_Data = size(Data_Set_Hists,1); %Number of cells measured in data
% Init_Param = Model_Results_Parameters(Nm,:); %Initial parameters from testing all models
rsl_FSP = log10(Init_Param); %Put initial parameters in log space 
proppdf = [];
delta = 0.02;
proprnd = @(x) x+delta*randn(size(x)).*(rand(size(x))>0.5);   % Proposal random sampler

%Objective function for MH
OBJ_FSP_MH = @(X)-get_FSP_OBJ(10.^X,Input,Nm,Output_Times,Data_Set_Hists,Mats_A,N,N_states);
%Objective function to calculate full probability distributions
OBJ_FSP_MH_P = @(X)get_FSP_OBJ(10.^X,Input,Nm,Output_Times,Data_Set_Hists,Mats_A,N,N_states,'P_RNA');
%Objective function for optimization before starting full MH
OBJ_FSP_FMIN = @(X)get_FSP_OBJ(10.^X,Input,Nm,Output_Times,Data_Set_Hists,Mats_A,N,N_states);
%Options for optimization before starting full MH
options =optimset('display','iter','MaxIter',500);
%Chain length for optimization to get good guess before starting full MH
N_Chain_search = 500;  N_Burn_search = 0; N_Thin_search = 1;   

figure();
for ipart = 1:10
    [Par_Chain_FSP,accept_FSP,Fun_Chain_FSP] = MetHaste(rsl_FSP,N_Chain_search,...
        'logpdf',OBJ_FSP_MH,'proppdf',proppdf,'proprnd',proprnd,...
        'N_Burn',N_Burn_search,'N_Thin',N_Thin_search,'symmetric','True');
    plot(Fun_Chain_FSP); hold on
    drawnow
    [~,J_FSP] = max(Fun_Chain_FSP);
    rsl_FSP = fminsearch(OBJ_FSP_FMIN,Par_Chain_FSP(J_FSP,:),options);
    if accept_FSP>0.2
        delta = delta*2;
        proprnd = @(x) x+delta*randn(size(x)).*(rand(size(x))>0.5);   %Proposal random sampler
    elseif accept_FSP<0.01
        delta = delta*0.8;
    end   
end

%% Perform the MH Search and Generate the Results 
STR_FSP_MH = ['generated_results/MetHast_',num2str(Nm),'_FSP_Results_I',num2str(input_num),'.csv'];  
fid = fopen(STR_FSP_MH,'w');
figure();
for ipart = 1:Npart_FSP
    waitbar(ipart/Npart_FSP)
    [Par_Chain_FSP,accept_FSP,Fun_Chain_FSP] = MetHaste(rsl_FSP,N_Chain,...
        'logpdf',OBJ_FSP_MH,'proppdf',proppdf,'proprnd',proprnd,...
        'N_Burn',0,'N_Thin',N_Thin,'symmetric','True');
    rsl_FSP=Par_Chain_FSP(end,:);
    plot(Fun_Chain_FSP); hold on
    drawnow
   
    for k=1:N_Chain
        for l = 1:size(Par_Chain_FSP,2)
            fprintf(fid,num2str(Par_Chain_FSP(k,l)));fprintf(fid,',');
        end
        fprintf(fid,num2str(Fun_Chain_FSP(k)));fprintf(fid,'\n ');
    end  
    
    mns_FSP(ipart,:) = mean(Par_Chain_FSP);   
end