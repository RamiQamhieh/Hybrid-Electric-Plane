%% aircraft calcs

% using numbers for regional jet currently. Will udpate with private jet as
% soon as I can.

% new numbers from Embraer phenom 100

k = 0.0366; % k,e,Cd0 from research paper linked in propulsion doc. 
e = 0.9213; % span efficiency factor
AR = 9.434; % derrived from k = (pi*e*AR)^-1
Cd0_TO = 0.077; % for takeoff
Cd0c = 0.024; % cruise

b = 31.83/15.913; % area divided by span
% based on cesna Citation II, looks like it has similar AR
a = 12.3; % wing span of Phenom
S = a*b; % wing area in m^2

v_TO = 4633/60; % takeoff speed = 150 knots/hr, converted to m/s;
v_y = 549/60; % steepest climb rate [m/s]
% these numbers are guesstimates from dubious sources.
gamma_TO = asin(v_y/v_TO);

T_TO = @(W) W*gamma_TO + 2*W*sqrt(k*Cd0_TO); % min thrust for maximum climb rate
sf = 1.3; % safety factor to ensure safe amount of takeoff thrust
Tmin = @(w) sf*T_T0(w);

R = 2182*1000; % range in km

rho0 = 1.225; % air density
rhoC = 0.4135; % density @ 35,000 ft, service ceiling is 41,00 ft 

fuel_TOL = 0.15; % portion of fuel used for takeoff and landing
vC = 750*1000/3600; % cruise speed in km/hr -> m/s

e_JA1 = 11.90; % energy density of commercial jet fuel in kWh/kg
e_Li = 0.3; % energy density for batteries = E/m

eff = 0.31; % energy efficiency of jet engine on 737

max_payload = 755; % kg, only 253 @ max fuel
max_fuel = 1272; % kg
fuel_MaxLoad = max_fuel - (755-253); % max fuel given max payload
max_M = 4800; % kg
max_dry = 4030;

W = @(m) 9.81*m; % don't forget to use this eqn for any W values in eqns.

P = @(v,w,S) 0.5*Cd0c*rhoC*S*v^3 + 2*k*w^2/rhoC/S/v; % P_req for cruise

FT = @(e,mf,v,w) e*mf*fuel_TOL / P(v,w); % flight time for electric power
% e = energy density, mf = mass of fuel = [kWhr/kg*kg/kW] = [hr]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% New stuff

PW_shaft = 4.4; % kW/kg , power to weight ratio of turboshaft generator for E-Fan X. Scale to meet power demand for our own designs

PW_Fan = 3.38; % same units, ratio for the motor + inverter + housing and fan. Based on the advertised power of the Siemens SP2000 
                % motor (2.5 MW), and weight of the standard engine is is replacing ( 740 kg) . Scale up for P_req to get engine mass.

e_2025 = 0.5; % projected future energy density for Li batteries

%% BAe 146

w146 = W(38810);

S146 = 77.3;

v146 = 747/3.6; % cruise in m/s

R146 = 500*1.609; % optimal range to capture most flights.

bae146 = [v146;w146;S146;R146];

% fuel tank holds 9430 kg, max payload is 8600 kg, 
% operating empty weight (OEW) = 23820 kg, cruise = 0.7M.

E = @(v, w, S, R) P(v,w,S) / 1000 * R / v / 3.6; % energy required in kWh
% v in m/s, w in N, S in m^2, R in km.

L =@(rho,s,v,cl) 1/2*rho*v^2*s*cl;

 % [~,~,r40] = hw2_1(12192,'SI'); This uses a rho function on my machine.
 r40 = 0.2977; % air density at 40,000 ft, in kg/m^3
v = 700/3.6;
s = 36*3;
LL = 380730;
l = LL/36;
cl = 1.481;

q = 1/2*r40*v^2;

CL = LL/q/s;

Cl = l/q/LL;

k = 1/pi/12/e;

cd0 = 0.021;
CD = @(CL) cd0 + k*CL^2;

D = q*s*CD(CL*1.3);

Preq = @(cd0,r,v,w,S) 0.5*cd0*r*S*v^3 + 2*k*w^2/r/S/v; % P_req for cruise

Pcruise = Preq(cd0,r40,v,LL,s);

e_vol = 1.386*1000; % energy density per metre^3, kWh/m^3, 
                    % converted from density per litre
                    
wingVol = 8.208; % m^3

fuelVol = wingVol*0.85; % usable wing volume to store fuel.

Emax_kWh = fuelVol*e_vol; % in kWh

t_hr = Emax_kWh*1000/Pcruise; % in hours

tChars = ["mass kg","cruise speed kph", "cruise alt m", "Pcruise MW", ...
            "Tcruise kN", "Emax kWh", "t_cruise hrs"];
tVals = [38810,v*3.6,12192,Pcruise/10^6,D/1000,Emax_kWh,t_hr];

m_kg = 38810; Vcruise_kph = v*3.6; alt_m = 12192; Pcruise_MW = Pcruise/10^6;
D_kN = D/1000; 

tt = table(m_kg,Vcruise_kph,alt_m,Pcruise_MW,D_kN,Emax_kWh,t_hr)
