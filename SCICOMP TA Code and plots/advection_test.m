%1 way wave equation
clear all; close all;
lambda = 1;

dx = 0.01;

dt = lambda*dx;

x = -1:dx:3;

t0 = 0;
tend = 2.4;

Nt = (tend - t0)/dt;

u0 = zeros(length(x),1);

for i=1:length(x)
    if (abs(x(i)) <= 1/2)
        u0(i) = cos(pi*x(i))^2;
    end
    
end

mdiag = ones(length(x),1);
ldiag = ones(length(x)-1,1);

D = lambda*(-diag(mdiag,0) + diag(ldiag,-1));


for j = 1:Nt
    unew = u0 - D*u0;

    plot(x,unew)
    drawnow
    % pause(0.05);

    u0 = unew;
end
