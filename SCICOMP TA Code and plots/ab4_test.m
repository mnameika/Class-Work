%4 step adams bashforth
clear all; close all;
% N = 10;

t0 = 0;
tend = 2;



u0 = [-3;-2;2];

f = @(u,t) [u(2); u(3); 4*t.^2 + 8*t - 10 - 4*u(1) - 4*u(2) - u(3)];
% u1 = u0;

% hrk = 0.0001;


% uSave = zeros(3,N);
% uSave(:,1) = u0;

% uInit = zeros(3,3);

uEx = @(t) -sin(2*t) + t.^2 - 3;

uPEx = @(t) -2*cos(2*t) + 2*t;

uP2Ex = @(t) 4*sin(2*t) + 2;


exVec = [uEx(2);uPEx(2);uP2Ex(2)];

numT = 20;
errSave = zeros(numT,1);

hSave = zeros(numT,1);

for j = 1:numT

    N = 2^j;
    t = linspace(0,2,N);

    h = t(2) - t(1);

    hSave(j) = h;

    uSave = zeros(3,N);
    uSave(:,1) = u0;

    for i = 1:N-1
    
        % ti = t(i);
    
        if i <= 3
            k1 = f(uSave(:,i),t(i));
            k2 = f(uSave(:,i) + h/2*k1, t(i) + h/2);
            k3 = f(uSave(:,i) + h/2*k2, t(i) + h/2);
            k4 = f(uSave(:,i) + h*k3, t(i) + h);
    
            unew = uSave(:,i) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    
            u1 = unew;
            uSave(:,i+1) = u1;
    
    
    
        else
            uSave(:,i+1) = uSave(:,i) + h/24*(-9*f(uSave(:,i-3),t(i-3)) + 37*f(uSave(:,i-2), t(i-2)) - 59*f(uSave(:,i-1), t(i-1)) + 55*f(uSave(:,i),t(i)));
        end
        
        
        
        
    end

    errSave(j) = norm(uSave(:,end) - exVec);
end


loglog(hSave, errSave, 'b-', hSave, hSave.^4, 'k-', 'linewidth', 1.1)
legend('error', '$h^4$', 'interpreter', 'latex')
% plot(t, uSave(1,:), t, uEx(t))