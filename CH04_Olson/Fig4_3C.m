clear; clc; close all;

f = @(x,t) sin(pi*x).*exp(-pi*pi*t);
dx = 0.05; dt = 0.0013;
beta = dt/dx^2;

x = 0:dx:1; t_final = 0.5;

M = 385; N = length(x);


% Exact solution
plot(x,f(x,t_final),'k','DisplayName','Exact solution');
% Forwad Euler
f0 = sin(pi*x); f1 = f0;
A = zeros(N-2,N-2);
for i = 1:N-2
    A(i,i) = 1-2*beta;
end

for i = 1: N-3
    A(i,i+1) = beta;
    A(i+1,i) = beta;
end
% for t = 0:dt:0.5-dt
%     f1(2:end-1) = beta*f0(1:end-2) + (1-2*beta)*f0(2:end-1) + beta*f0(3:end);
%     f0 = f1;
% end
for i = 1:385
    f1(2:end-1) = A*f0(2:end-1)';
    f0 = f1;
end
hold on;
plot(x, f1,'ko','DisplayName','Forward Euler');



% Backward Euler
f0 = sin(pi*x); f1 = f0;
A = zeros(N-2,N-2);
for i = 1:N-2
    A(i,i) = 1+2*beta;
end

for i = 1: N-3
    A(i,i+1) = -beta;
    A(i+1,i) = -beta;
end

for i = 1:385
    f1(2:end-1) = A\f0(2:end-1)';
    f0 = f1;
end
hold on;
plot(x, f1,'kx','DisplayName','Backward Euler');

legend1 = legend('show');
set(legend1,'Position',[0.35 0.32 0.375 0.22])

xlabel('X');
ylabel('C(x,t)');

set(findall(gcf,'-property','FontSize'),'FontSize',24)
