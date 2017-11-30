% qBio Chapter Model Code (w/ kinase-dependent k)
% This function analyses the derivative of X at time point T 
function [xdot] = Model_ODE(T, x, Parameters, Input, Nm)

%Bring in initial parameters specified in run code
k12 = Parameters(1); k23 = Parameters(2); 
k21 = Parameters(3); k32 = Parameters(4); 
kr2 = Parameters(5); kr3 = Parameters(6); 
b = Parameters(7); g = Parameters(8);

%Modify whichever k is kinase-dependent based on Model Number, Nm
switch Nm
    case 1
        k12=max(0,-k12+b*Input(T));
    case 2
        k23=max(0,-k23+b*Input(T));
    case 3
        k21=max(0,k21-b*Input(T));
    case 4
        k32=max(0,k32-b*Input(T));
end

%Define ODEs for the gene regulation model
xdot(4,1) = kr2*x(2) + kr3*x(3) - g*x(4); %ODE for #RNA Number
xdot(1,1) = k21*x(2) - k12*x(1); %State 1 (inactive) ODE
xdot(2,1) = k12*x(1) - (k21+k23)*x(2) + k32*x(3); %State 2(active) ODE
xdot(3,1) = k23*x(2) - k32*x(3); %State 3(active) ODE
end

