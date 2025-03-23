A = [4 2; 2 2];
b = [3;0];
x0 = [2;2];

maxCount = 1000;

for i=1:maxCount
    r = b - A*x0;
    gamma = r'*r/(r'*A*r);

    xnew = x0 + gamma*r;

    x0 = xnew;
    
    
    
end
x0