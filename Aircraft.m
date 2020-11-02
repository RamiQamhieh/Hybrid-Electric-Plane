classdef Aircraft < handle
    %Aircraft properties class
    
    properties (Access = public)
        %    STATIC PROPERTIES
        % Bookkeeping
        Name = "bae"
        
        options = optimset('Display','off');
        
        % Wing parameters
        c = 3; %[m] chord length
        b = 36; %[m] wing span
        S ;
        AR ;
        
        % engine parameters
        Ne = 4; % # of engines
        d = 2.7555; % prop diameter in m
        BPR = 21; % bypass ratio TODO!
        bhp = 31; % engine brake horsepower TODO!
        % Dp = 3.75; % [m] propeller diameter TODO: this is expected in [ft]
        gamma_min = 0.030; % for 4-engine
        is_a_jet = false;% jet or prop??? (for takeoff analysis)
        has_turboprops = false; % props or turboprops (if not is_a_jet)
        has_thrust_reversers = false; %(p.677) 
        thrust_reverser_cutoff_speed = 93; %[km/h] typical, varies based on engine
        idle_thrust = 0;
        max_forward_thrust = 31;
        static_foward_thrust = 31;
        Vmax; %[km/h] max vel
        aux = 1.25; % multiplier to run auxiliary systems. IE hydraulics, etc. increase fuel consumption by 25%
        Pmax = 6500; % max capability of power plant in kW    
        
        % Simulation parameters
        tstep = 30/3600; %[h] 15 sec -> hrs; time step size
        
        % Other parameters

        Vcruise = 650; %[km/h] cruise velocity
        cruise_altitude = 35000/3.2808; % [m]
        RHOcruise = 0.3796; % [kg/m^3] cruise air density (35,000 ft)
        Mf0 = 9430; %[kg] initial fuel mass
        fratio = 0.15; % 0 = all battery, 1 = all fuel
        Mfuel; % state var for fuel mass
        MJA1; % state var for jet fuel mass
        Espent = 0;
        Vclimb_vert = 8.894529759;% <- real data %4;%10.16; % [m/s]  vertical climb rate [10.16 = 2000ft/s]       
        Vclimb = 77.78; %[m/s] climb velocity magnitude (this should be 1.2*Vstall)
        
        % Prices, fuel properties
        cost_fuel_mass = 0.56; %[$/kg] fuel cost
        cost_elec = 0.12; %[$/kWh] electricity cost
        dens_batt = 0.5; %[kWh/kg] battery energy density (optimistic for future)
        dens_fuel = 11.90*0.31; %[kWh/kg] fuel energy density (JA-1 jetfuel) 
        cost_fuel;
        flight_cost;
        init_kWh_avail;

        % Airfoil parameters : to be updated by wingChar
        CLmax ;
        CLA ; % CL as a function of wing root abs alpha
        ACL ; % alpha as a function of CL
        alphaA; % wing root abs alpha
        LG = true; % landing gear bool: true means it is delpoyed, 
        % false means it is retracted. call self.CD0_() to update CD0 whenever
        % the landing gear state is changed
        CD0;
        k; % 1 / pi / e / AR
        e; % span efficiency factor
        
        % Environment/Takeoff/Landing parameters
        g = 9.81;
        mu = 0.03; % [] runway rolling friction coefficient (0.03 assumes dry concrete) p.672
        SG; % ground roll distance
        SR; % rotation ground roll distance (p.673)
        STR; % transition distace (p.674)
        SC; % climb distance (p.674)
        h_obstacle = 35*0.3048; % [m] for commercial
        rho_SL = 1.225; %[kg/m3] (rho @ sea level) 
        T_takeoff_static = 31*4; % based on BAE 146 numbers
        
        % Takeoff outputs
        BFL; % balanced field length (p.675)     
        TO_d; % takeoff distance, time, energy
        TO_t;
        TO_E;
        
        % Landing outputs
        E_land; % energy used in landing
        FAR;
        SA;
        SF;
        SGR_b;
        SGR;
        
        %    DYNAMIC PROPERTIES
        fuel_left = 1; %0 = no fuel, 1 = full
        batt_left = 1; %0 = no battery, 1 = full
        kWh_avail; % Fuel kW hours avaliable
        time_elapsed = 0; % [s] total time since plane started moving for takeoff 
        
        eta = 1; % prop efficiency
        x; % horizontal position in meters
        mission; % horizontal flight distance
        T = 0; % thrust kN
        P = 0; % power kW, should never exceed ~4.5ish mega watts
        V = 0; % vel m/s
        CL = 0; 
        CDi = 0; 
        CD;
        M;
        W;
        altitude=0; % m
        rho;%[kg/m3] (sea level)
        pressure;
        temp_env;
        Vstall; % stall velocity m/s
        
        % status
        is_not_stalling;
        stage; % stage of flight
        
        % save data (r_x = record of x)
        r_M = [];
        r_time_elapsed = [];
        r_altitude = [];
        r_T = [];
        r_P = [];
        r_V = [];
        r_eta = [];
        r_kWh_avail = [];
        r_CL = []; 
        r_CDi = []; 
        r_CD = [];
        r_stage = []; % stage of flight
        r_Mfuel = [];
        r_x = [];
      
    end
    
    methods
        
        function init_state(self)
            % call when initializing
            self.Mfuel = self.Mf0; % this var will change over flight
            self.MJA1 = self.Mf0*self.fratio; % mass of jet fuel
            self.S_();
            self.T_();
            self.P_();
            self.AR_();
            self.CD0_();
            self.t_p_rho();
            self.kWh_avail_();
            self.init_kWh_avail = self.kWh_avail;
            self.M = 36584; % initial mass in kg
            self.W_(); % weight in N
            self.cost_fuel = self.cost_fuel_mass/self.dens_fuel; %[$/kWh]
            self.CDi = 0;
            self.CD = self.CD0 + self.CDi;
            self.wingChar();
            self.stall();
            self.x = 0;
            
            self.record_state();
        end
        
        function update_state(self)
            % call at every update
            self.W_();
            self.T_();
            self.eta_();
            self.P_();
            self.alf_(); % sets alpha and CL
            self.CD0_();   
            self.CDi_();
            self.stall();
            self.check_stall();
            self.t_p_rho();
            self.kWh_avail_();
            self.stall();
            self.x_();

            self.record_state();

        end  
        
        function x_(self)
            self.x = self.x + self.V*self.tstep*3600;
        end
        
        function CDi_(self)
            self.CDi = self.k*self.CL^2;
        end
        
        
        function wingChar(self)
            v = self.Vcruise/3.6;
            twist = -2*pi/180;
            K = 12; % density of discrete points along wing

            % numerical solution from textbook

            cr = 4; ct = 2; % chord at root and tip
            lambda = ct/cr;
            j = 1:K;

            q = 1/2*self.RHOcruise*v^2;
            CLWF = self.W/q/self.S ;

            th = pi.*j ./ 2./K ; % range of theta sweeping across the wing
            cosj = cos(th);
            sinj = sin(th);
            Y = cosj;
            C = 1 - (1-lambda) * Y ; % chordlength as a function of theta

            % compute some coeff.

            i = 2.*j - 1;
            D = zeros(K);
            for j = 1:K
                D1 = 1/C(j);
                D2 = pi / ( self.AR*(1 + lambda) *sinj(j) ) ;

                for n = 1:K
                    D(j,n) = (D1 + D2*i(n)) * sin(th(j)*i(n));
                end

            end


            al1 = 3*pi/180; al2 = 6*pi/180 ; % two abs alpha vals to solve some eqns

            alabs = (al1 + twist.*cosj ) ; % absloute angle of attack

            A = D \ alabs' ; % coefficients An

            CL1 = zeros(1,K);

            CLW1 = A(1)*pi^2 / (1 + lambda) ;
            for j = 1:K
                Sum = 0;
                for n = 1:K
                    Sum = Sum + A(n) *sin(i(n)*th(j));
                end

                CL1(j) = 2*pi/C(j) * Sum;

            end
            % repeat for second alpha
            alabs = (al2 + twist.*cosj ) ;

            A = D \ alabs' ;

            CL2 = zeros(1,K);

            CLW2 = A(1)*pi^2 / (1 + lambda) ;
            for j = 1:K
                Sum = 0;
                for n = 1:K
                    Sum = Sum + A(n) *sin(i(n)*th(j));
                end

                CL2(j) = 2*pi/C(j) * Sum;

            end

            CLA = ( CL2 - CL1 ) ./ ( CLW2 - CLW1 ) ;
            CLB = CL1 - CLA.*CLW1;

            alf = al1 + (al2 - al1) * (CLWF - CLW1) / (CLW2 - CLW1) ;
            % I think this is the alpha needed to fly at given speeed under given load.
            alabs = alf + twist.*cosj;

            A = D \ alabs' ;

            Cl = CLB + CLA.*CLWF;

            ALind = zeros(1,K);
            CDind = zeros(1,K);
            M = zeros(1,K);

            for j = 1:K
                Sum = 0;

                for n = 1:K
                    Sum = Sum + i(n)*A(n)* sin(i(n)*th(j)) / sinj(j);
                end

                ALind(j) = Sum * -pi / (self.AR*(1+lambda)) ;
                CDind(j) = -Cl(j) * ALind(j);
                M(j) = Cl(j) / alabs(j);
            end

            Sum = 0;

            for j = 1:K
                Sum = Sum + i(j)*A(j)^2;
            end

            %CDi = pi^3 / (self.AR*(1 + lambda)^2) *Sum;

            sigma = 0;
            for j = 2:K
                sigma = sigma + j*A(j)^2/A(1)^2;
            end

            self.e = 1/(1+sigma);
            self.k = 1/pi/self.e/self.AR;

            Sum = 0;
            for j = 1:K-1
                Sum = Sum + ( M(j)*C(j) + M(j+1)*C(j+1) ) * ( Y(j) - Y(j+1) ) ;
            end

            %Mbar = Sum / (1 + lambda);
            %alabsW = CLWF / Mbar * 180/pi;
            self.alphaA = alf * 180/pi;
            %CLi = Cl;
            self.CL = CLWF;
            
            mr = (CLW2 - CLW1) / 3;
            br = CLW1 - 3*mr;

            % find CL as a function of the absolute angle of attack at the wing root
            % keep in mind the absolute angle of attack is larger than the regular
            % angle of attack
            self.CLA = @(alpha) br + alpha*mr;
            self.ACL = @(CL) (CL - br) / mr; 
            self.CLmax = self.CLA(15);

        end
        
        
        function step_sim(self)
            self.time_elapsed = self.time_elapsed + self.tstep; % in hr
            self.M_tstep();
        end
        
        function cruise(self) 
            % only run this after wingChar so that it has the correct values to pull from

            self.V = self.Vcruise/3.6;
            self.rho = self.RHOcruise;
            self.altitude = self.cruise_altitude;
            
            duration = (self.mission - self.x)/(self.Vcruise/3.6)/3600;
            % cruise duration in hr
            for n = 1:round(duration/self.tstep)
                self.step_sim();
                self.update_state();
            end
            

        end
        
        function mission_(self,d)
            self.mission = d*1609; % assume d is in miles, convert to m
            
            self.init_state();
            fprintf("\n takeoff starting \n");
            self.take_off();
            fprintf("ascent starting\n");
            self.r_stage(end) = 1;
            self.ascend(self.cruise_altitude);
            
            fprintf("cruise starting\n");
            self.r_stage(end) = 2;
            self.cruise();
                        
            fprintf("descent starting\n")
            self.r_stage(end) = 3;
            self.descend(self.h_obstacle);

            fprintf("landing starting\n");
            self.r_stage(end) = 4;
            self.Landing();
            
            kWh_req = self.init_kWh_avail - self.kWh_avail;
            self.flight_cost_(kWh_req);
            
        end
        
        function record_state(self)
            % save data (r_x = record of x)
            self.r_M = [self.r_M self.M];
            self.r_time_elapsed = [self.r_time_elapsed self.time_elapsed];
            self.r_altitude = [self.r_altitude self.altitude];
            self.r_T = [self.r_T self.T];
            self.r_P = [self.r_P self.P];
            self.r_V = [self.r_V self.V];
            self.r_eta = [self.r_eta self.eta];
            self.r_kWh_avail = [self.r_kWh_avail self.kWh_avail];
            self.r_CL = [self.r_CL self.CL]; 
            self.r_CDi = [self.r_CDi self.CDi]; 
            self.r_CD = [self.r_CD self.CD];
            self.r_stage = [self.r_stage 0];
            self.r_Mfuel = [self.r_Mfuel self.Mfuel];
            self.r_x = [self.r_x self.x];

        end
        
        function eta_(self)
            A = (self.d/2)^2*pi;
            ve = @(v) self.rho*A*v*(v-self.V) - self.T*1000; 
            if self.V > 10
                % quadratic equation relating cruise speed, exhaust speed and thrust_req
                Ve = fsolve(ve,self.V,self.options);
                self.eta = 2/(1 + Ve/self.V); 
            else
                self.eta = 0.1; %what should this be?
            end
        end
        
        function alf_(self)
            self.CL = self.W/self.V^2/self.S/self.rho*2;
            self.alphaA = self.ACL(self.CL); % set wing root angle of attack
        end
        
        function M_tstep(self) % calcs energy expend in 1 timestep -> updates M
            
            Ereq = self.aux*self.P*self.tstep; % kW*hr
            self.Espent = self.Espent + Ereq;
            fspent = Ereq*(self.Mf0*self.fratio*self.dens_fuel/self.init_kWh_avail)/self.dens_fuel; % kg of fuel spent
            self.Mfuel = self.Mfuel - fspent;
            self.M = self.M - fspent;
            self.MJA1 = self.MJA1 - fspent;
        end
        
        
        function M_Ereq(self, Ereq) % uses given energy expend -> updates M
            Ereq = Ereq*self.aux;
            fspent = Ereq*(self.Mf0*self.fratio*self.dens_fuel/self.init_kWh_avail)/self.dens_fuel; % kg of fuel spent
            self.Mfuel = self.Mfuel - fspent;
            self.M = self.M - fspent;
            self.MJA1 = self.MJA1 - fspent;
            self.kWh_avail_();
        end
        
        function ascend(self, altitude)
            Pm = self.Pmax; % in kW
            Vy = self.Vclimb_vert;
            while self.altitude < altitude
                % update env.
                self.t_p_rho();
                self.stall();
                self.eta_();
                
                % climb vel calcs
                self.Vclimb = 1.2*self.Vstall;
                self.V = self.Vclimb;
                gamma_climb = sin(Vy/self.Vclimb);
                self.Vclimb = self.V*sin(gamma_climb);

                % power req. Steady Aircraft text p. 150 (where Vclimb is vertical vel)
                Preq = self.W*Vy + 0.5*self.rho*self.V^ 3*self.S*self.CD0 ...
                    + (2*self.k*self.W^2)/(self.rho*self.V*self.S); 
                Preq = Preq/(1000); %W to kW
                
                while Preq/self.eta > Pm
                    Vy = 0.9*Vy;
                    gamma_climb = sin(Vy/self.Vclimb);
                    self.Vclimb = self.V*sin(gamma_climb);
                    Preq = 1/1000*(self.W*Vy + 0.5*self.rho*self.V^ 3*self.S*self.CD0 ...
                        + (2*self.k*self.W^2)/(self.rho*self.V*self.S) ); 
                end

                % update altitude
                self.altitude = self.altitude + self.tstep*3600*Vy;
                if self.altitude > altitude % in case too far, for plot
                    self.altitude = altitude;
                end
                
                % update pos
                V_horiz = sqrt(self.V^2 - Vy^2);
                self.x = self.x + V_horiz*self.tstep*3600;
                
                % updates P usage and time
                self.P = Preq/self.eta;
                self.step_sim();
                self.record_state();
            end
        end

        function descend(self, altitude)
            while self.altitude > altitude
                % update env.
                self.t_p_rho();
                self.stall();
                self.eta_();
                
                % climb vel calcs
                self.Vclimb = 1.2*self.Vstall;
                self.V = self.Vclimb;
                gamma_climb = sin(self.Vclimb_vert/self.Vclimb);
                self.Vclimb = self.V*sin(gamma_climb);

                % power req. Steady Aircraft text p. 150 (where Vclimb is vertical vel)
                Preq = self.W*-1*self.Vclimb_vert + 0.5*self.rho*self.V^ 3*self.S*self.CD0 ...
                    + (2*self.k*self.W^2)/(self.rho*self.V*self.S); 
                Preq = Preq/(1000); %W to kW
                if Preq < 0
                    Preq = 0; % don't allow harvesting
                end

                % update altitude
                self.altitude = self.altitude - self.tstep*3600*self.Vclimb_vert;
                if self.altitude < altitude % in case too far, for plot
                    self.altitude = altitude;
                end
                
                % update pos
                V_horiz = sqrt(self.V^2 - self.Vclimb_vert^2);
                self.x = self.x + V_horiz*self.tstep*3600;
                
                
                % updates P usage and time
                self.P = Preq/self.eta;
                self.step_sim();
                self.record_state();
            end
        end
        
        function S_(self) 
            % Planform area (rectangular - don't do this for final)
            self.S = self.b*self.c;
        end
        
        function T_(self)
            % Thrust required for steady level flight (steady flight text p.95)
            if self.V > 10     
                self.T = 1/1000 * ( 0.5*self.rho*self.V^2*self.S*self.CD0 + ...
                    2*self.k*self.W^2/(self.rho*self.V^2*self.S)  ); % kN
            else 
                self.T = self.T_takeoff_static;
            end
        end
        
        function W_(self)
            self.W = self.M*9.81; % N
        end

        function P_(self)
            % Power - - - - - (- 97)
            self.P = self.T*self.V / self.eta;
            % TODO: add P for climbing, find min climbing rate for takeoff
        end
        
        function stall(self)
            % Stall constraint (sft p.95)
            self.Vstall = sqrt((2*self.W)/(self.rho*self.S*self.CLmax));
        end
        
        function AR_(self)
            self.AR = self.b/self.c;
        end
        
        function CD0_(self)
            if self.LG == true
                self.CD0 = 0.07;
            end
            if self.LG == false
                self.CD0 = 0.021;
            end
        end
        
        
        function check_stall(self)
            if self.Vstall < self.V
                self.is_not_stalling = true;
            else
                self.is_not_stalling = false;
            end
        end
        
        function kWh_avail_(self)
            self.kWh_avail = ( self.MJA1 *self.dens_fuel ) + ...
                ( self.MJA1 /self.fratio * self.dens_batt);
        end
        
        function flight_cost_(self, kWh_req)
            kWh_batt = kWh_req*(1-self.fratio);
            kWh_fuel = kWh_req*self.fratio;

            elec_cost = kWh_batt*self.cost_elec;
            fuel_cost = kWh_fuel*self.cost_fuel;

            self.flight_cost = elec_cost + fuel_cost;
        end
        
        function t_p_rho(self)
            [self.temp_env,self.pressure,self.rho]=T_P_rho(self.altitude,'SI');
        end
        
        
        function take_off(self)
            %Ground Roll (p.672)
            % acceleration
            a_ = @(V, CL) self.g.*(((self.T.*1000./self.W) - self.mu) + ...
                ((self.rho)./(2.*self.W./self.S)).*(1.*self.CD0 - ...
                self.k.*self.CL.^2 + self.mu.*self.CL).*V.^2); % m/s^2
            % assume CL = CLmax
            a = @(V) a_(V, self.CLmax*0.9);
            %a(self.V)
            % ground-roll distance
            fun = @(V) V./a(V);
            
            Vi = self.V;
            Vf = self.Vstall*1.1; % Takeoff velocity (should be this or more...?)
            
            self.SG = integral(fun,Vi,Vf);
            self.V = (1/sqrt(2))*Vf; %avg during ground roll
            self.update_state();
            
            P_G = self.P;
            time_G = self.SG/self.V; % time spent in ground-roll
            E_used_G = P_G*time_G/3600; % kWh
            self.time_elapsed = self.time_elapsed + time_G/3600;
            self.M_Ereq(E_used_G);
            self.update_state();

            % Rotation Ground Roll Distance (p.673)
            % about 1 second of rotation -> SR = VTO
            self.V = Vf;
            self.update_state();
            P_SR = self.P;
            self.SR = self.V;
            time_SR = 1;%self.SR/self.V; % s
            E_used_SR = P_SR*time_SR/3600;% kWh
            self.M_Ereq(E_used_SR);
            self.time_elapsed = self.time_elapsed + time_SR/3600;
            self.update_state();

            % Transition (p.647)
            CL = 0.9*self.CLmax; %avg
            VTR = 1.5*self.Vstall; %avg
            self.V = VTR;
            self.update_state();
            P_TR = self.P;
            R = (VTR^2)/(0.2*self.g); %transition arc radius
            self.STR = sqrt(R^2 - (R - self.h_obstacle)); % transition distance
            time_TR = self.STR/self.V;
            E_used_TR = P_TR*time_TR/3600; % kWh
            self.M_Ereq(E_used_TR);
            self.time_elapsed = self.time_elapsed + time_TR/3600;
        
            % update altitude after transition
            gamma_climb = asin(self.STR/R); % might have misinterpreted [eq 17.109] 
            hTR = R*(1 - cos(gamma_climb));
            self.altitude = hTR;
            self.update_state();
            
            % climb (p.674)
            self.SC = (self.h_obstacle - hTR)/(tan(gamma_climb));
            if self.SC < 0
                self.SC = 0; %obs. was cleared during takeoff
            end
            VC = 1.2*self.Vstall;
            self.update_state();
            P_C = self.P;
            time_C = self.SC/self.V;
            E_used_C =P_C*time_C/3600;
            self.M_Ereq(E_used_C);
            self.time_elapsed = self.time_elapsed + time_C/3600;
            self.altitude = self.h_obstacle;
            self.update_state();
            
            % Balanced Field Length (BFL) (p.675)
            %       total takeoff distance including obstacle clearance
            U = 0.01*self.CLmax + 0.02; % for flaps in takeoff position
            % gamma_climb = asin((self.T-self.D)/self.W); % 1-engine-out, climb speed
            G = gamma_climb - self.gamma_min;
            CL_climb = self.W/(VC)^2/self.S/self.rho*2;
            
            Whp = 0.00134102; % x watts * Whp = y hp
            hp = 1.25e6*Whp*4; 
            % not sure if hp is supposed to be per engine or total. Each engine outputs 1.25MW
            S = self.S*3.28084^2; % convert wing area to ft^2
            rhoSL = 1.225*0.00194032 ; % convert Sea level air density to slugs/ft^3
            g = 32.17; % gravity in ft/s^2
            Dp = self.d*3.28084; % convert prop diam to ft
            RHO = self.rho*0.00194032; 
            % convert rho to slug/ft^3. I believe it needs to be slugs for dimensions to work
            
            % if you see self.M*2.2, that is me converting kg to lbs
            
            if self.is_a_jet
                Tavg = 0.75*self.T_takeoff_static*((5 + self.BPR)/(4 + self.BPR));
            else
                Tavg = 5.75*hp*(((RHO/rhoSL)*self.Ne*(Dp)^2)/(hp))^(1/3) 
                % TODO: this is in hp, FIX! (possible other units issues)
            end 
            
            self.BFL = (   (0.863/(1 + 2.3*G))* ...
                ((self.M*2.2/S)/(RHO*g*CL_climb)  + self.h_obstacle*3.28084)* ...
                ((1/(Tavg/(self.M*2.2)-U)) + 2.7) + ...
                (655/(sqrt(RHO/rhoSL)))   ) / 3.28084 % convert from ft to m
            
            self.TO_d = self.SG + self.SR + self.STR + self.SC;
            
            self.TO_t = time_G + time_SR + time_TR + time_C;
            self.TO_E = E_used_C + E_used_SR + E_used_G + E_used_TR;
            P_C;
            P_SR;
            P_G;
            P_TR;
            
            %self.M_Ereq(self.TO_E); already counted
            self.W_();
            self.LG = 0;
            self.CD0_();
            
            self.x = 0 + self.TO_d; 
            
            
        end
        
        function Landing(self)
            % Landing analysis Pg. 676
            % Trade Study Pg. 724
            self.altitude = self.h_obstacle;
            self.update_state();
            
            % flare (p.677)
            VTD = 1.15*self.Vstall; % touchdown speed
            Vf = 1.23*self.Vstall; % flare velocity (avg)
            Va = 1.3*self.Vstall; % approach speed
            n = 1.2; % typical aircraft
            R = Vf^2/(0.2*self.g); % Pg. 674
            self.V = Va;
            self.update_state();
            gamma_climb = asin((self.T/self.W) - (1/(self.CL/self.CD))); % approach angle
            hTF = R*(1-cos(gamma_climb)); % flare height
            
            % approach distance
            self.SA = (self.h_obstacle - hTF)/(tan(gamma_climb));
            self.V = Va;
            self.update_state();
            P_A = self.P;
            time_A = self.SA/self.V;
            E_used_A =P_A*time_A/3600; %kWh
            Vx = self.V*cos(gamma_climb);
            dt = self.SA/Vx;
            self.time_elapsed = self.time_elapsed + dt/3600;
            self.altitude = hTF;
            self.update_state();
            
            % Pg. 677, horizontal flare distance
            % self.SF = R*sin(gamma_climb);
            self.V = Vf;
            Vx = self.V*cos(gamma_climb);
            % dt = self.SF/Vx;
            dt = Vf;
            self.SF = Vx*dt;
            self.time_elapsed = self.time_elapsed + dt/3600;
            self.altitude = 0;
            self.update_state();
            P_F = self.P;
            time_F = self.SF/self.V;
            E_used_F = P_F*time_F/3600;

            % Note (p.677) - "Although the deceleration from Va to VTD
            % would imply additional energy and thus additional distance, 
            % this is negligible because the pilot usually pulls
            % off all remaining approach power when the flare is begun."
            
            % ground roll (p.677)
            
            % ground roll delay
            GR_delay = 1.5; % [s] ground roll delay (assume 1-3)
            SGR_b = VTD*GR_delay; % ground roll distance before applying brakes
            self.V = VTD;
            dt = GR_delay;
            self.time_elapsed = self.time_elapsed + dt/3600;
            self.update_state();
            P_TD = self.P;
            time_TD = SGR_b/self.V;
            E_used_TD =P_TD*time_TD/3600;
            
            % braking distance
            Vi = VTD;
            Vf = 0;
            mu = 0.5; % rolling resist, for civil (p.678)
            
            T = self.idle_thrust;
            if self.has_thrust_reversers
                if self.is_a_jet
                    T = 0.45*self.max_forward_thrust;
                else
                    if self.has_turboprops
                        T = 0.40*self.static_foward_thrust;
                    else
                        T = 0.60*self.static_foward_thrust;  
                    end
                end
            end
            self.T = T;
            self.V = (1/sqrt(2))*Vi;
            self.P_();
            
            a_ = @(V, CL) self.g.*(((self.T./self.W) - self.mu) + ...
            ((self.rho)./(2.*self.W./self.S)).*(1.*self.CD0 - ...
            self.k.*self.CL.^2 + self.mu.*self.CL).*V.^2);
             
            % assume CL = CLmax
            a = @(V) a_(V, self.CLmax);
            
            % braking ground-roll distance
            fun = @(V) V./a(V);
            self.SGR = integral(fun,Vi,Vf);
            
            P_BGR = self.P;
%             time_BGR = self.SGR/self.V;
            time_BGR = self.BFL/self.V;
            E_used_BGR =P_BGR*time_BGR/3600; %kWh
            self.E_land = E_used_A + E_used_BGR + E_used_F + E_used_TD;
            
            self.V = 0;
            
            self.FAR = 1.666*(self.SA+self.SF+SGR_b+self.SGR); % Pg. 678
            self.altitude = 0; % assume landing at sea level
            self.M_Ereq(self.E_land);
            self.time_elapsed = self.time_elapsed + time_BGR/3600;
            self.update_state();
        end
% Landing end
   
        
        

    end
end
% % Range (p.643)
% % R = @(V,L,C,D,Wi,Wf) = ((V*L)/(C*D))*ln(Wi/Wf);
% % Range (Electric) (p.757)
% R = @(L,D,Esb,Nb2s,Np,mb,g,m) 3.6*(L/D)*(Esb*Nb2s*Np/g)*(mb/m);
% 
% % Endurance (Electric) (p.755)
% E = @(mb, Esb, Nb2s, Pused) (mb*Esb*Nb2s)/(1000*Pused);
% % where:
% % % mb = mass of batteries {kg}
% % % Esb = battery specific energy {wh/kg}
% % % Nb2s = total system efficiency from battery to motor output shaft
% % % Pused = average power used during that period of time {kW}
