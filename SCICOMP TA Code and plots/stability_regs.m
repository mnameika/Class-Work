%AB Stability region
clear all; close all;

x0 = -6.5;
x1 = 1.5;

y0 = -4;
y1 = 4;

dx = 0.0025;
dy = 0.0025;

x = x0:dx:x1;

y = y0:dy:y1;

[Y,X] = meshgrid(y,x);

Nx = length(x);
Ny = length(y);

plotMat = zeros(Nx,Ny);



%2-Step AM
% p = @(z) [1-5*z/12, -(1 + 8*z/12), z/12];


%3-Step AM
% p = @(z) [1 - 9/24*z, -(1 + 19/24*z), 5/24*z, -z/24];

%4-Step AM
% p = @(z) [1 - 251/720*z, -(1 + 646/720*z), 264/720*z, -106/720*z, 19/720*z];

%5-Step AM
p = @(z) [1 - 95/288*z, -(1 + 1427/1440*z), 133/240*z, -241/720*z, 173/1440*z, -3/160*z];


%five step AB
% p = @(z) [1, -(1 + 1901/720*z), 2774/720*z, -2616/720*z, 1274/720*z, -251/720*z];

%four step AB
% p = @(z) [1, -(1 + 55/24*z), 59/24*z, -37/24*z, 9/24*z];

%three step AB
% p = @(z) [1 -(1 + 23/12*z) 4/3*z -5/12*z];

%two step AB
% p = @(z) [1 , -(1+3/2*z), z/2];

truVal = false;
truValB = false;

tol = 5*1e-4;
 
for i=1:Nx

    for j = 1:Ny
        z = x(i) + 1i*y(j);
        
        rootVec = abs(roots(p(z)));

        truVecL = rootVec <= 1;

        
        if (sum(truVecL) == length(truVecL))
            truVal = true;
        else 
            truVal = false;
        end

        if (truVal)


            plotMat(i,j) = 1/2;
            
            % plot(x(i), y(j), 'k.', 'markersize', 20);
            % axis([x0 x1 y0 y1])
            % hold on
        end

        if (truValB)
            plotMat(i,j) = 1;
        end
        


        if (x(i) == 0 || y(j) == 0)
            plotMat(i,j) = 1;

        end


        truValB = false;
    end
end

mesh(X,Y,plotMat);
colorMat = [linspace(1,0)', linspace(1,0)', linspace(1,0)'];
colormap(colorMat)
xlabel('Re$(z)$', 'fontsize', 25, 'interpreter', 'latex')
ylabel('Im$(z)$', 'fontsize', 25, 'interpreter', 'latex')
grid on
view(0,90)