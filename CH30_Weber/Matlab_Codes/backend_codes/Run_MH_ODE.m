%Run_MH_ODE: Runs Metropolis-Hastings with the ODE Objective Function in Logspace

disp(' ')
disp('Running the MH search on the gene expression models.')
disp('This could take as much as a few hours.')
disp('You may wish to use shorter MH chains or fewer skips during the debug stage.')
disp('Then once it is working, you can leave it to run while you do something else.')

N_Data = size(Data_Set,1); %Number of cells measured in data  
Init_Param = eval(['Model_Results_Parameters_I',num2str(input_num),'(Nm,:)']);
rsl = log10(Init_Param); %Put initial parameters in logspace
proppdf = [];
delta = 0.02;
proprnd = @(x) x+delta*randn(size(x)).*(rand(size(x))>0.5);   %Proposal random sampler

%Objective function for MH
OBJ_MH = @(x)-Difference_ODE(10.^x,Input,Nm,Output_Times,x0,Data_Set); 
%Objective function for optimization before starting full MH
OBJ_FMIN = @(x)Difference_ODE(10.^x,Input,Nm,Output_Times,x0,Data_Set); 
%Options for optimization before starting full MH
options =optimset('display','iter','MaxIter',500); 
%Chain length for optimization to get good guess before starting full MH
N_Chain_search = 500;  N_Burn_search = 0; N_Thin_search = 1; 

% Optimization to get good parameter guesses for the initial MH
figure();
for ipart = 1:10
    [Par_Chain,accept,Fun_Chain] = MetHaste(rsl,N_Chain_search,...
        'logpdf',OBJ_MH,'proppdf',proppdf,'proprnd',proprnd,...
        'N_Burn',N_Burn_search,'N_Thin',N_Thin_search,'symmetric','True');
    plot(Fun_Chain); hold on
    drawnow
    [~,J] = max(Fun_Chain);
    rsl = fminsearch(OBJ_FMIN,Par_Chain(J,:),options);
    if accept>0.2
        delta = delta*2;
        proprnd = @(x) x+delta*randn(size(x)).*(rand(size(x))>0.5);   %Proposal random sampler
    elseif accept<0.01
        delta = delta*0.8;
    end   
end

%% Perform the MH Search and Generate the Results 
STR_ODE_MH = ['generated_results/MetHast_',num2str(Nm),'_ODE_Results_I',num2str(input_num),'.csv'];  
fid = fopen(STR_ODE_MH,'w');
figure();
for ipart = 1:Npart_ODE
    waitbar(ipart/Npart_ODE)
    [Par_Chain,accept,Fun_Chain] = MetHaste(rsl,N_Chain,...
        'logpdf',OBJ_MH,'proppdf',proppdf,'proprnd',proprnd,...
        'N_Burn',0,'N_Thin',N_Thin,'symmetric','True');
    rsl=Par_Chain(end,:);
    plot(Fun_Chain); hold on
    drawnow
   
    for k=1:N_Chain
        for l = 1:size(Par_Chain,2)
            fprintf(fid,num2str(Par_Chain(k,l)));fprintf(fid,',');
        end
        fprintf(fid,num2str(Fun_Chain(k)));fprintf(fid,'\n ');
    end
    
    mns(ipart,:) = mean(Par_Chain);   
end


