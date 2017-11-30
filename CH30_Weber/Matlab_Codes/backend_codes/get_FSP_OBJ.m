% This function outputs either the FSP objective function for the
% Metropolis-Hastings search of parameter space or outputs full probability
% distributions from the FSP analysis

function [the_output] = get_FSP_OBJ(Parameters,Input,Nm,Output_Times,Data_Set_Hists,Mats_A,N,N_states,whichoutput)
% Modifications to parameters and A matrix based on Model Number
switch Nm
     case 1
        %Define parameters
        k12_0 = Parameters(1); k23 = Parameters(2); k21 = Parameters(3); k32 = Parameters(4);
        kr2 = Parameters(5); kr3 = Parameters(6); b = Parameters(7); g = Parameters(8);
        k12=@(t)max(0,-k12_0+b*Input(t)); %Time-dependent transition rate
        
        % Constant terms of A matrix (not time dependent)
        Acon = k23*Mats_A.A23 + k32*Mats_A.A32 + k21*Mats_A.A21 + kr2*Mats_A.Ar2 + kr3*Mats_A.Ar3 + g*Mats_A.Ag;
        Jacfun = @(t,x)(Acon+k12(t)*Mats_A.A12);
        odefun_fsp = @(t,x)(Acon+k12(t)*Mats_A.A12)*x;
        J_Pattern = isfinite(abs(Mats_A.A12+Acon));
    case 2
        %Define parameters
        k12 = Parameters(1); k23_0 = Parameters(2); k21 = Parameters(3); k32 = Parameters(4);
        kr2 = Parameters(5); kr3 = Parameters(6); b = Parameters(7); g = Parameters(8);
        k23=@(t)max(0,-k23_0+b*Input(t)); %Time-dependent transition rate
        
        % Constant terms of A matrix (not time dependent)
        Acon = k12*Mats_A.A12 + k32*Mats_A.A32 + k21*Mats_A.A21 + kr2*Mats_A.Ar2 + kr3*Mats_A.Ar3 + g*Mats_A.Ag;
        Jacfun = @(t,x)(Acon+k23(t)*Mats_A.A23);
        odefun_fsp = @(t,x)(Acon+k23(t)*Mats_A.A23)*x;
        J_Pattern = isfinite(abs(Mats_A.A23+Acon));
     case 3
        %Define parameters
        k12 = Parameters(1); k23 = Parameters(2); k21_0 = Parameters(3); k32 = Parameters(4);
        kr2 = Parameters(5); kr3 = Parameters(6); b = Parameters(7); g = Parameters(8);
        k21=@(t)max(0,k21_0-b*Input(t)); %Time-dependent transition rate
        
        % Constant terms of A matrix (not time dependent)
        Acon = k12*Mats_A.A12 + k23*Mats_A.A23 + k32*Mats_A.A32 + kr2*Mats_A.Ar2 + kr3*Mats_A.Ar3 + g*Mats_A.Ag;
        Jacfun = @(t,x)(Acon+k21(t)*Mats_A.A21);
        odefun_fsp = @(t,x)(Acon+k21(t)*Mats_A.A21)*x;
        J_Pattern = isfinite(abs(Mats_A.A21+Acon));
     case 4
        %Define parameters
        k12 = Parameters(1); k23 = Parameters(2); k21 = Parameters(3); k32_0 = Parameters(4);
        kr2 = Parameters(5); kr3 = Parameters(6); b = Parameters(7); g = Parameters(8);
        k32=@(t)max(0,k32_0-b*Input(t)); %Time-dependent transition rate
        
        % Constant terms of A matrix (not time dependent)
        Acon = k12*Mats_A.A12 + k23*Mats_A.A23 + k21*Mats_A.A21 + kr2*Mats_A.Ar2 + kr3*Mats_A.Ar3 + g*Mats_A.Ag;
        Jacfun = @(t,x)(Acon+k32(t)*Mats_A.A32);
        odefun_fsp = @(t,x)(Acon+k32(t)*Mats_A.A32)*x;
        J_Pattern = isfinite(abs(Mats_A.A32+Acon));
end

P0 = zeros(N_states*(N+1),1); P0(1)=1;
options = odeset('JPattern',J_Pattern,'Jacobian',Jacfun,'MaxStep',1);
[~,XOUT] = ode15s(odefun_fsp,Output_Times,P0,options);

% Full probability distribution
P_RNA = zeros(length(Output_Times),N+1);
for j=1:N_states
    P_RNA(:,:) = P_RNA(:,:) + XOUT(:,j:N_states:end);
end

% Data_Set_Hist has same form as P_RNA
if ~exist('whichoutput','var')    
    Nd = length(Data_Set_Hists(1,:));  % Max response in any cell.
    TMP = Data_Set_Hists.*log(P_RNA(:,1:Nd));
    OBJ = -sum(sum(TMP(Data_Set_Hists~=0)));
    % Account for parameter priors (lognormal prior with mean 10^0 and relative deviation +/- 10^2.
    LogPar = log10(Parameters);
    ErrLogPar = LogPar.^2/8;
    OBJ = OBJ + sum(ErrLogPar);
    the_output=OBJ;
elseif strcmp(whichoutput,'OBJ')
     Nd = length(Data_Set_Hists(1,:));  % Max response in any cell.
    TMP = Data_Set_Hists.*log(P_RNA(:,1:Nd));
    OBJ = -sum(sum(TMP(Data_Set_Hists~=0)));
    % Account for parameter priors (lognormal prior with mean 10^0 and relative deviation +/- 10^2.
    LogPar = log10(Parameters);
    ErrLogPar = LogPar.^2/8;
    OBJ = OBJ + sum(ErrLogPar);
    the_output=OBJ;
elseif strcmp(whichoutput,'P_RNA') 
    the_output=P_RNA;  % Full mRNA probability distributions.
end    
end

    
