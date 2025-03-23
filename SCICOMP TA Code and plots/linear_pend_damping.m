%Numerically solving the damped linear pendulum problem using (1) Midpoint,
%(2) Trapezoidal, and (3) 2-Step Adams-Bashforth
% clear all; close all;

%pendulum constant
a = 100;

%damping coefficient
b = 10;

%Choose LMM
%1 --> Midpoint
%2 --> Trapezoidal
%3 --> Adams-Bashforth 2-Step
method = 2;

%initial and final times
t0 = 0;
tEnd = 10;

%Right-hand-side vector
f = @(u) [u(2), -a*u(1) - b*u(2)];

%initial condition
u0 = [1,0];

%number of time steps to check
N = 15;

%initialize error vector
errVec = zeros(N,1);

%initialize step size vector
hVec = zeros(N,1);

%constants in the solution -- will change if u0(2) =/= 0
c1 = u0(1);
c2 = (2*u0(2) + b*c1)/sqrt(4*a - b^2);

%Initialize matrix to hold eigenvalues of A
eigSave = zeros(N,2);

%operator matrix and its eigenvalues
A = [0 1; -a -b];
A_eigs = eig(A);

for j = 1:N
    %time step
    % dt = 2^(-j);
    m = 2^j;
    dt = 1/(m+1);

    %compute z-values
    eigSave(j,:) = dt*A_eigs;

    %save step size
    hVec(j) = dt;
    
    %number of time steps
    Nt = tEnd/dt + 1;
    
    %time vector
    t = t0:dt:tEnd;
    
    %initialize numerical solution
    uVec = zeros(Nt,2);
    uVec(1,:) = u0;
    
    %time stepping
    for i = 1:Nt-1
        %start with one step of Euler
    
        %one time step of Euler to start the method
        if (i == 1)
            % uVec(i+1,:) = uVec(i,:) + dt*f(uVec(i,:));
            uVec(i+1,:) = uVec(i,:) + dt*f(uVec(i,:) + dt/2*f(uVec(i,:)));
        else
            %step into the methods
            switch method
                case 1
                    mthdStr = "Midpoint";

                    % uVec(i + 1,:) = uVec(i,:) + dt*f(uVec(i,:)+dt/2*f(uVec(i,:)));
                    uVec(i+1,:) = uVec(i-1,:) + 2*dt*f(uVec(i,:));
                case 2
                    mthdStr = "Trapezoidal";
                    %predictor step -- necessary since it's an implicit
                    %method

                    uPred = uVec(i,:) + dt*f(uVec(i,:));
                    uVec(i+1,:) = uVec(i,:) + dt/2*(f(uVec(i,:)) + f(uPred));
                   
                case 3 
                    mthdStr = "AB2";

                    uVec(i+1,:) = uVec(i,:) + dt/2*(-f(uVec(i-1,:)) + 3*f(uVec(i,:)));
            end %end switch statement
        end %end if statement
    end %end for -- time stepping

    %exact solution
    uEx = exp(-b/2*t).*(c1*cos(1/2*sqrt(4*a - b^2)*t) + c2*sin(1/2*sqrt(4*a - b^2)*t));

    %compute error
    errVec(j) = norm(uVec(:,1).' - uEx, inf);

end %end for -- error saving

%plot errors 
loglog(hVec,errVec, 'b--o', hVec, hVec.^2, 'k-', 'linewidth', 1.2)
xlabel('$\Delta t$', 'fontsize', 20, 'interpreter', 'latex')
ylabel('$\|\Theta_{exact}(t) - \Theta_{approx}(t)\|_{\infty}$', 'fontsize', 20, 'interpreter', 'latex')
grid on
legend('$\|\Theta_{exact} - \Theta_{approx}\|_{\infty}$', '$\Delta t^2$' , 'fontsize', 15, 'interpreter', 'latex', 'location', 'northwest')
title(mthdStr + " Method", 'fontsize', 25, 'interpreter', 'latex')

figure(2)
plot(uVec(:,1))