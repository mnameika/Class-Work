clear all; close all;
a = 0;
b = 2*pi;

alpha = 0;
beta = 0;

f = @(x) -sin(x);


h = 0.1;

x = h:h:2*pi;

fVec = f(x).';

fVec(1) = fVec(1) - alpha/h^2;
fVec(end) = fVec(end) - beta/h^2;

N = length(x);

onesVec = ones(N-1,1);
twosVec = -2*ones(N,1);

M = 1/h^2*diag(twosVec);
% N = 1/h^2*(diag(onesVec,-1) + diag(onesVec,1));
L = -1/h^2*diag(onesVec,-1);
U = -1/h^2*diag(onesVec,1);
u0 = x.*(2*pi - x);

u0 = u0.';

tol = 1e-5;
maxCount = 10000;
count = 0;

err = 10;
% plot(x,u0, 'b--o', x, sin(x), 'k-')
% axis tight
% drawnow

while (count < maxCount && err > tol)
    
    unew = M \ (fVec +  (U + L) * u0);
    
    count = count + 1;

    err = norm(unew - u0);
    u0 = unew;
    % 
    plot(x, u0,'b--o', x, sin(x),'k-')
    axis tight
    drawnow

    % pause(0.3)
    
end
plot(x, u0,'b--o', x, sin(x),'k-')