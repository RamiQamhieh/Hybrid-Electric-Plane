function [alpha,CL,CLi,CDi] = wingChar(v,rho,m,twist,k)
% Inputs : v [kph] , rho [kg/m^3] , m [kg] , twist [deg] , k [int]
% twist is the geometric twist angle of the wing about the span axis. K is
% the discrete point density. 
%
% Outputs are [alpha,CL,CLi,CDi] = root angle of attack, total CL, CL for
% each discrete segment, and total induced drag coefficient.
%
% For zero args, v = 700, r = 0.2977, m = 38810, twist = -4, k = 12
%
% Also works for 3 or 4 args, using the defaults when needed.
if(nargin < 5)
    k = 12;
end

if nargin < 4
    twist = -4; % geometric angle of twist at wing tip
end

if nargin == 0
    v = 700; % v in m/s
    rho = 0.2977; % air density at 40,000 ft, in kg/m^3
    m = 38810;
end

w = m*9.81;
v = v/3.6;
twist = twist*pi/180;

% analytical solution from textbook

cr = 4; ct = 2; % chord at root and tip
lambda = ct/cr ; % ratio of chords
b = 36;% span
c = 3; % avg chord
S = b*c;
AR = b/c;
j = 1:k;

q = 1/2*rho*v^2;
CLWF = w/q/S ;

th = pi.*j ./ 2./k ;
cosj = cos(th);
sinj = sin(th);
Y = cosj;
C = 1 - (1-lambda) * Y ;

% compute some coeff.

i = 2.*j - 1;
D = zeros(k);
for j = 1:k
    D1 = 1/C(j);
    D2 = pi / ( AR*(1 + lambda) *sinj(j) ) ;
    
    for n = 1:k
        D(j,n) = (D1 + D2*i(n)) * sin(th(j)*i(n));
    end
    
end


al1 = 3*pi/180; al2 = 6*pi/180 ; % two abs alpha vals to solve some eqns

alabs = (al1 + twist.*cosj ) ;

A = D \ alabs' ;

CL1 = zeros(1,k);

CLW1 = A(1)*pi^2 / (1 + lambda) ;
for j = 1:k
    sum = 0;
    for n = 1:k
        sum = sum + A(n) *sin(i(n)*th(j));
    end
    
    CL1(j) = 2*pi/C(j) * sum;
    
end
% repeat for second alpha
alabs = (al2 + twist.*cosj ) ;

A = D \ alabs' ;

CL2 = zeros(1,k);

CLW2 = A(1)*pi^2 / (1 + lambda) ;
for j = 1:k
    sum = 0;
    for n = 1:k
        sum = sum + A(n) *sin(i(n)*th(j));
    end
    
    CL2(j) = 2*pi/C(j) * sum;
    
end

CLA = ( CL2 - CL1 ) ./ ( CLW2 - CLW1 ) ;
CLB = CL1 - CLA.*CLW1;

alf = al1 + (al2 - al1) * (CLWF - CLW1) / (CLW2 - CLW1) ;
% I think this is the alpha needed to fly at given speeed under given load.
alabs = alf + twist.*cosj;

A = D \ alabs' ;

Cl = CLB + CLA.*CLWF;

ALind = zeros(1,k);
CDind = zeros(1,k);
M = zeros(1,k);

for j = 1:k
    sum = 0;
    
    for n = 1:k
        sum = sum + i(n)*A(n)* sin(i(n)*th(j)) / sinj(j);
    end
    
    ALind(j) = sum * -pi / (AR*(1+lambda)) ;
    CDind(j) = -Cl(j) * ALind(j);
    M(j) = Cl(j) / alabs(j);
end

sum = 0;

for j = 1:k
    sum = sum + i(j)*A(j)^2;
end

CDi = pi^3 / (AR*(1 + lambda)^2) *sum;

sum = 0;
for j = 1:k-1
    sum = sum + ( M(j)*C(j) + M(j+1)*C(j+1) ) * ( Y(j) - Y(j+1) ) ;
end

Mbar = sum / (1 + lambda) ;
alabsW = CLWF / Mbar * 180/pi;
alpha = alf * 180/pi;
CLi = Cl;
CL = CLWF;

end
