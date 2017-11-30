% This function runs a single SSA trajectory and records the entire X_Array,
% which includes all species: x = [S1 S2 S3 m]'
function [X_Array] = SSA(ModFun, x0, Output_Times)
Nt = size(Output_Times,2); % Number of time points (from Output_Times)
t = Output_Times(1);  % Initial time - 1st Output_Time
tstop = Output_Times(Nt); % Stop time - Last Output_time
N_species = size(x0,1); % Number of Species
X_Array = zeros(N_species,Nt);  % Preallocate trajectory for efficiency
x = x0;  % Initial Condition for x = [S1 S2 S3 m]'
count = 1; % Initialize counter

while t<tstop  
    [~,~,P,S] = ModFun(t,x); % Get stoichiometry and linear propensity functions from model function
    w0 = sum(P); % Compute the sum of the propensity functions
    t = t+1/w0*log(1/rand); % Update time of next reaction (where tau = 1/w0*log(1/rand))
    if t<=tstop                      
        while t>Output_Times(count)     
            X_Array(:,count) = x;
            count = count+1; 
        end       
        r2w0=rand*w0; % Generate second random number and multiply by propensity sum        
        i=1; % Initialize reaction counter
        while sum(P(1:i))<r2w0 % Increment counter until sum of propensities exceeds r2w0
            i=i+1; % Update the counter
        end
        x = x+S(:,i); % Update the state 
    end
end

% Record x results into matrix
X_Array(:,count:end) = repmat(x,1,Nt-(count-1));
end



