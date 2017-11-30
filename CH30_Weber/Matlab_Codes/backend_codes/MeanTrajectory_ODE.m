% This function obtains ODE solutions for all states
function [M] = MeanTrajectory_ODE(Parameters,Nm,Input,Output_Times,x0)
%Begin by defining ODE function of the system parameters(F) as a function of x and t.
F=@(t,x)Model_ODE(t, x, Parameters, Input, Nm);
%Plug into ODE solver and solve
[X Y]=ode15s(F,Output_Times,x0); 
%Obtain solutions and output mean number of mRNA
M=Y(:,4)'; 
end
