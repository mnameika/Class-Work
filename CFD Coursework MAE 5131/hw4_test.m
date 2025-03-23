%FEM upwind
clear all; close all;
dx = 0.01;
% x = 0:dx:5;
xend = 5;
x = [0,dx/2:dx:xend-dx/2, xend];


u = zeros(length(x),1);

dt = dx/2;


tend = 2;
c = 1;

uEx = @(t) exp(-100*(x - 3/10 - c*t).^2);

for i = 1:length(x)
    if (0 <= x(i) && x(i) <= 0.6)
       u(i) = exp(-100*(x(i) - 3/10).^2);

    elseif (0.6 < x(i) && x(i) <= 0.8)
        u(i) = 1;
    end
end

nu = dt/(dx);

u(1) = 0;
t0 = 0;

onesVec = ones(length(x)-1,1);
cVec = ones(length(x),1);

%this might be the center difference scheme
A = -nu*(diag(cVec, 0) - diag(onesVec,-1));
% A(1,1) = -1;
% A(end,end) = -1;
% A(end,end-1) = 1;

Nt = tend/dt;


for j=1:Nt
    t0 = t0 + dt;
    unew = u + A*u;
    u = unew;
    plot(x,u)
    drawnow
    % pause(0.1)
    
    
    
end

% plot(x, u)
% plot(x,u,x, uEx(t0))