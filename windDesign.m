%% Wing calculations

%% Geometry

alpha = 6*pi/180;

critM = 1/1.575; % critical mach number for our airfoil

m0 = 2*pi;

Cd0c = 0.024; % cruise

cs = 3; % avg chord
b = 36; % span in m
S = cs*b; % wing area

AR = b/cs;

xk = [-0.508
-0.499
-0.477
-0.442
-0.397
-0.341
-0.277
-0.205
-0.127
-0.045
0.038
0.121
0.203
0.279
0.348
0.407
0.453
0.484
0.499]';

xk = xk + 0.499; 

yk = [0.004
0.005
0.0065
0.0085
0.011
0.0145
0.018
0.021
0.024
0.0255
0.026
0.025
0.022
0.018
0.013
0.008
0.0035
0.0005
0]';

[a,E,Y] = LSQPoly(xk,yk,2);

plot(xk,yk,xk,Y)

X = @(theta) cs/2*(1-cos(theta)); % transform span distance into theta
% dzdx = @(x) a(2) + 2.*a(3).*X(x) + 3.*a(4).*X(x).^2 + 4.*a(5).*X(x).^3;
% dzdxN = @(x,n) cos(n*x).*( a(2) + 2.*a(3).*X(x) + 3.*a(4).*X(x).^2 + 4.*a(5).*X(x).^3);

dzdx = @(x) a(2) + 2.*a(3).*X(x);
dzdxN = @(x,n) cos(n*x).*(a(2) + 2.*a(3).*X(x) );


dzdxTH = @(x) dzdxN(x,1);

%% Analysis

A0 = @(dzdx) alpha - 1/pi * integral(dzdx,0,pi);
An = @(n,dzdxTH) 2/pi * ( integral(dzdxTH,0,pi)); 
% dzdxTH is the product of dz/dx and cos(N*theta)

gamma = @(V,An,A0,theta) 2*V( A0 * (1 + cos(theta)) / sin(theta) + An );
cl = @(alpha , dzdx ) 2*pi*( alpha + 1/pi * integral(dzdx,0,pi) );

CL = m0*cs*pi*b/4/S * An(1,dzdxTH);
CD = 1.05 * (Cd0c + CL^2/pi/AR );
