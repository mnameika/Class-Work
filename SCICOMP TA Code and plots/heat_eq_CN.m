clear all; close all;
N = 8;

errVec = zeros(N,1);
dxVec = zeros(N,1);
for j = 1:N
    b = 0.1;
    % Nx = 10*2^j;
    
    xEnd = 10;
    
    x = -xEnd:1/(2^j*10):xEnd;
    x = x';
    Nx = length(x);
    dx = x(2) - x(1);
    dxVec(j) = dx;
    

    dt = dx*b/2;
    tEnd = 1/2;
    
    Nt = tEnd/dt;
    
    
    
    u0 = zeros(Nx,1);
    a = 1;
    
    uEx = @(x,t) 1/2*(erf((a - x)/sqrt(4*b*t)) + erf((a + x)/sqrt(4*b*t)));
    
    
    
    for i = 1:Nx
        if (abs(x(i)) <= a)
            u0(i) = 1;
        end
    end
    
    diagVec = ones(Nx,1);
    upLowDiag = ones(Nx-1,1);
    
    A = diag((1/dt + b/dx^2)*diagVec,0) + diag(-b/(2*dx^2)*upLowDiag, -1) + diag(-b/(2*dx^2)*upLowDiag,1);
    
    B = diag((1/dt - b/dx^2)*diagVec,0) + diag(b/(2*dx^2)*upLowDiag,-1) + diag(b/(2*dx^2)*upLowDiag,1);
    t0 = 0;
    for i=1:Nt
        t0 = t0 + dt;
        
    
        RHS = B*u0;
        uNew = A\RHS;
    
        u0 = uNew;
       
        % plot(x,u0, 'b--')
        % hold on
        % plot(x, uEx(x,t0), '-k')
        % drawnow
        % hold off
        % pause(0.1)
        
    end
    
    errVec(j) = sqrt(dx)*norm(uEx(x,t0) - u0);
end

loglog(dxVec, errVec, 'b--', dxVec, dxVec.^2, 'k-', dxVec, dxVec, 'r-')