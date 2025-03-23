%pendulum problem
clear all; close all;
N = 100;



a = 0;
b = 2*pi;

h = (b - a)/(N+1);

alpha = 0.7;
beta = 0.7;

x = a:h:b;

x = x.';

theta0 = 0.7 + sin(x/2);

onesVec = ones(N-1,1);
twosVec = -2*ones(N,1);

A0 = (diag(onesVec, -1) + diag(onesVec, 1) + diag(twosVec,0));

J0 = A0 + h^2*diag(cos(theta0));




G0 = A0*theta0 + sin(theta0);
% G0(1) = G0(1) + alpha/h^2;
% G0(end) = G0(end) + beta/h^2;

tol = 1e-10;
maxCount = 100;
err = 10;
count = 0;

% x = a:h:b;

while (err > tol && count < maxCount)
    delta = inv(J0)*(-G0);

    delta(1) = 0;
    delta(end) = 0;
    
    theta0 = theta0 + delta;

    % A0 = 1/h^2*(diag(onesVec, -1) + diag(onesVec, 1) + diag(twosVec));

    J0 = A0 + diag(cos(theta0));
   
    J(1,:) = 0;
    J(1,1) = 1;
    J(end,:) = 0;
    J(end,end) = 1;

    G0 = A0*theta0 + sin(theta0);

    plot(x, theta0)
    drawnow
    % pause(0.5)
    
    
    count = count + 1;
end