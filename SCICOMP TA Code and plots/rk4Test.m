clear all; close all;

u0 = [-3;-2;2];

t0 = 0;
tend = 2;

f = @(u,t) [u(2); u(3); 4*t.^2 + 8*t - 10 - u(3) - 4*u(2) - 4*u(1)];

vEx = @(t) -sin(2*t) + t.^2 - 3;
vPrimeEx = @(t) -2*cos(2*t) + 2*t;
v2PrimeEx = @(t) 4*sin(2*t) + 2;

exVec = [vEx(2); vPrimeEx(2); v2PrimeEx(2)];

% uApprox = zeros(N,1);
% uApprox(1) = u0(1);

errorHold = zeros(15,1);
hVec = zeros(15,1);


for j=1:15
    u1 = u0;
    N = 2^j;



    uApprox = zeros(N,1);
    uApprox(1) = u0(1);
    t = linspace(0,tend,N);

    h = t(2) - t(1);

    hVec(j) = h;

    for i=1:N-1
        t0 = t(i);
    
        % %method listed on Wikipedia
        % k1 = f(u1,t0);
        % k2 = f(u1 + h/2*k1, t0 + h/2);
        % k3 = f(u1 + h/2*k2, t0 + h/2);
        % k4 = f(u1 + h*k3, t0 + h);
        % 
        % un = u1 + h/6*(k1 + 2*k2 + 2*k3 + k4);
        % uApprox(i+1) = un(1);
        % u1 = un;
    
        %method listed in the text
        Y1 = u1;
        Y2 = u1 + 1/2*h*f(Y1,t0);
        Y3 = u1 + 1/2*h*f(Y2,t0 + h/2);
        Y4 = u1 + h*f(Y3, t0 + h/2);

        un = u1 + h/6*(f(Y1,t0) + 2*f(Y2,t0 + h/2) + 2*f(Y3, t0 + h/2) + f(Y4, t0 + h));
        uApprox(i) = un(1);
        u1 = un;
    end

    plot(t,uApprox, 'b--', t, vEx(t),'k', 'linewidth', 1.1)
    grid on

    

    errorHold(j) = norm(u1 - exVec);
    % err = norm(u1 - exVec)

end


loglog(hVec,errorHold, 'b--', hVec, hVec.^4, 'k', 'linewidth', 1.1)
grid on
axis tight
xlabel('$h$', 'fontsize', 25, 'interpreter', 'latex')
ylabel('$\|u(T) - U^T\|_2$', 'fontsize', 25, 'interpreter', 'latex')
legend('$\|u(T) - U^T\|_2(h)$', '$h^4$', 'fontsize', 15, 'interpreter', 'latex', 'location', 'northwest')