%% Generate SSA Data Using FSP
OBJ_FSP_MH_P = @(X)get_FSP_OBJ(10.^X,Input,Nm,Output_Times,Data_Set_Hists,Mats_A,N,N_states,'P_RNA');
P_true = OBJ_FSP_MH_P(log10(abs(Parameters_True)));
Nt = size(P_true,1); %Total number of time points
Ncells = 100; %Total number of cells
Data_Set = zeros(Ncells,Nt); %Preallocate for data
Data_Set_Hists = zeros(Nt,50); %Preallocate for data histograms
for it = 1:Nt
    for ic = 1:Ncells
        P_cumul = cumsum(P_true(it,:));
        r = rand;
        j = find(P_cumul>r,1,'first');
        Data_Set(ic,it) = j-1;
        Data_Set_Hists(it,j) = Data_Set_Hists(it,j)+1;
    end
end
fid = ['sim_data_bimodal_fsp.csv'];  
dlmwrite(fid,Data_Set);
fid = ['sim_data_bimodal_fsp_Data_Set_Hists.csv'];  
dlmwrite(fid,Data_Set_Hists);

