% Fundamental Plasma Functions - CGS gaussian units
%
% All inputs are in cgs gaussion units except Temperature which is in [eV]
% Use cgs2SI and SI2cgs to convert between Sia nd cgs units
% Use physicalConstants-SI and physicalConstants-cgs for values of
% constants
%
% EXAMPLE USAGE:
%
% B = SI2cgs(5,'Magnetic Field'); % Magnetic field, Gauss
% n = SI2cgs(1e20, 'Density'); % Electron density [cm^-3]
% F = FundamentalPlasma(); % create object
% F.getElectronGyroFreq(B); % Electron gyro frequency, [s^-1]
% cgs2SI(F.getIonGyroRadius(Ti,B,Z,A),'Length')) % Ion gyro radis in [m]
%
% Reference(s):
%  Shea, J. J. (2003). A plasma formulary for physics, technology, and astrophysics [Book Review].
%  IEEE Electrical Insulation Magazine, 19(1), 52?52. https://doi.org/10.1109/mei.2003.1178121

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2021 Rishabh Datta, MIT
%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef FundamentalPlasma
    properties
    end
    methods
        function obj = FundamentalPlasma()            
        end
        % Frequencies [rad/s]
        function out = getElectronGyroFreq(~,B) % Electron gyro-frequency
            % Returns electron gyro-frequency
            % Inputs:
            %   B = magnetic field [Gauss]
            % w_ce = eB/me.c [rad/s]
            out = 1.76e7 * B; % [rad/s]
        end 
        function out = getIonGyroFreq(~,B,Z,A) % Ion gyro-frequency
            % Inputs:
            %   B [Gauss] = Magnetic field
            %   Z [-] = Ionization
            %   A [-] = Atomic weight, ion mass / proton mass
            out = 9.58e3 .* Z ./ A .* B; % [rad/s]
        end 
        function out = getElecPlasmaFreq(~,ne) % Electron plasma freq.
            % Inputs:
            %       ne [cm^-3] = electron number density
            out = 5.64e4 * sqrt(ne); % [rad/s]
        end
        function out = getIonPlasmaFreq(~,ni,Z,A) % Ion plasma freq.
            % Inputs:
            %   ni [cm^-3] = ion number density
            %   Z [-] = Ionization
            %   A [-] = Atomic weight, ion mass / proton mass
            out = 1.32e3 .* Z .* A .^ (-0.5) .* sqrt(ni); % [rad/s]
        end
        function out = getElElCollFreq(obj,ne,Te) % Electron-Electron Collison Freq
            % ne [cm^-3] = electron number density
            % Te [eV] = electron temperature
            lnA = obj.getCoulombLog(ne,Te); % coulomb logarithm
            out = 2.91e-6 .* ne .* lnA .* Te.^(-3/2); % 1/sec
        end 
        function out = getElIonCollFreq(obj,ne,Te) % Electrin-Ion Collision freq
            out = 2 * obj.getElElCollFreq(ne,Te); % [s^-1]
        end 
        function out = getIonIonCollFreq(obj,ni,Z,A,Ti,Te) % Ion-Ion Collison Freq
            % Inputs:
            %   ni [cm^-3] = Ion numbe density
            %   Z [-] = Ionization
            %   A [-] = Atomic weight, ion mass / proton mass
            %   Ti, Te [eV] = Ion, Electron temperature
            lnA = obj.getCoulombLog(Z*ni,Te); % coulomb logarithm
            out = 4.8e-8 .* Z.^4 .* A^(-1/2) .* ni .* lnA .* Ti.^(-3/2); 
        end
        function out = getIonElCollFreq(obj,ne,Te,A) % Ion-Electron Collison Freq
            out = A * obj.getElIonCollFreq(ne,Te); 
        end
        function out = getIonElecEqmFreq(obj,ne,Te,A,Z)
                % ne = electron denisty [cm^-3]
                % Te = electron temp. [eV]
                % A = ion mass [-]
                % Z = ionization [-]
                lnA = obj.getCoulombLog(ne,Te); % coulomb logarithm
                out = 3.2e-9 .* ne .* Z.^2 .* lnA ./ (A .* Te.^(3/2)); % [s^-1]
        end
        function out = getElecIonEqmFreq(obj,ni,Te,A,Z)
                % ni = electron denisty [cm^-3]
                % Te = electron temp. [eV]
                % A = ion mass [-]
                % Z = ionization [-]
                ne = Z * ni;
                lnA = obj.getCoulombLog(ne,Te); % coulomb logarithm
                out = 3.2e-9 .* ni .* Z.^2 .* lnA ./ (A .* Te.^(3/2)); % [s^-1]
        end
        % Characteristic Times
        function out = getElElCollTime(obj,ne,Te)
            out = 1 ./ obj.getElElCollFreq(ne,Te); % [s]
        end
        function out = getIonIonCollTime(obj,ni,Z,A,Ti,Te)
            out = 1 ./ obj.getIonIonCollFreq(ni,Z,A,Ti,Te); % [s]
        end
        function out = getEnergyEqmTime(obj,ne,Te,A,Zbar) % Energy Equilibriation time
            % INPUTS:
            %       ne [cm^-3] = electron number density
            %       Te [eV] = electron temperature
            %       A [-] = atomic mass, ion mass / proton mass
            %       Zbar [-] = ionization
            load physicalConstants-cgs me mp
            ni = ne / Zbar;
            out = 1 / obj.getElecIonEqmFreq(ni,Te,A,Zbar);
            % out = Zbar^-1 * 1/2 * A * mp / me .* obj.getElElCollTime(ne,Te); 
        end
        % Lengths [cm]
        function out = getElecDeBrogLen(~,Te) % Electron de Broglie length
            % Inputs:
            %   Te [eV] = electron temp.
            out = 2.76e-8 .* Te.^(-0.5); % [cm]
        end
        function out = getClasMinApproach(~,Te) % Classic dist. of minimum approach
            % Inputs:
            %   Te [eV] = electron temp.
            out = 1.44e-7 .* (1./Te); % [cm], electron temperature
        end
        function out = getElecGyroRadius(~,Te,B) % Electron Gyro-radius
            % Inputs:
            %   Te [eV] = electron temperature
            %   B [Gauss] = magnetic field
             out = sqrt(2) * 2.38 .* Te .^ 0.5 * (1./B); % [cm]
        end
        function out = getIonGyroRadius(~,Ti,B,Z,A) % Ion Gyro-Radius
            % Inputs:
            %   Ti [eV] = ion temperature
            %   B [Gauss] = magnetic field
            %   Z [-] = Ionization
            %   A [-] = Atomic weight, ion mass / proton mass
            out = sqrt(2) * 1.02e2 * A.^(0.5) * Z.^-1 * Ti.^0.5 * B.^-1; % [cm]
        end 
        function out = getDebyeLen(~,n,T) % Debye length
            % Inputs:
            % n [cm^-3] = elec. number density
            % T [eV] = temperature
            out = sqrt(1) * 7.43e2 * T.^0.5 .* n.^(-0.5); % [cm]
        end
        function out = getElecSkinDepth(obj,ne) % Electron Skin Depth
            % Inputs:
            %       ne [cm^-3] = electron number density
            load physicalConstants-cgs c
            out = c ./  obj.getElecPlasmaFreq(ne); % cm
        end
        function out = getIonSkinDepth(obj,ni,Z,A) % Ion Skin depth
            % Inputs:
            %   ni [cm^-3] = ion number density
            %   Z [-] = Ionization
            %   A [-] = Atomic weight, ion mass / proton mass
            load physicalConstants-cgs c
            out = c /  obj.getIonPlasmaFreq(ni,Z,A); % [cm]
        end
        function out = getElElMfp(obj,ne,Te) % Electron-Electron mean free Path
            % Inputs:
            %   ne [cm^-3] = electron number density
            %   Te [eV] = electron temperature
            out = obj.getElecThermalVel(Te) .* obj.getElElCollTime(ne,Te); % [cm]
        end
        function out = getIonIonMfp(obj,ni,Z,A,Ti,Te) % Ion-Ion mfp
            % Inputs:
            %   ni [cm^-3] = Ion numbe density
            %   Z [-] = Ionization
            %   A [-] = Atomic weight, ion mass / proton mass
            %   Ti, Te [eV] = Ion, Electron temperature
            out = obj.getIonThermalVel(Ti,A) .* obj.getIonIonCollTime(ni,Z,A,Ti,Te);  
        end
        function out = getElecInerLen(~,ne) % Electron inertial length
            % Inputs:
            %   ne [cm^-3] = electron number density
            out = 5.31e5 * ne.^(-0.5); % [cm]
        end
        function out = getIonInerLen(~,ni,Z,A) % Ion inertial length
            % Inputs:
            %   ni [cm^-3] = ion number density
            %   Z [-] = Ionization
            %   A [-] = Atomic weight, ion mass / proton mass
            out = 2.28e7 * Z.^(-1) .* (A ./ ni)^0.5; % [cm] 
        end
        % Velocities [cm/s]
        function out = getElecThermalVel(~,Te) % Electron thermal velocity
            % Inputs:
            %   Te [eV] = Electron temp.
            out = 1 * 4.19e7 * Te.^0.5; % [cm/s]
        end
        function out = getIonThermalVel(~,Ti,A) % Ion thermal velocity
            % Inputs:
            %   Ti [eV] = ion temp.
            %   A [-] = Atomic weight, ion mass / proton mass
            out = 1 * 9.79e5 .* A^-0.5 .* Ti.^0.5; % [cm/s]
        end
        function out = getIonSoundSpeed(~,Te,Z,A,gamma) % Ion sound speed
            % Inputs:
            %   Te [eV] = electron temp.
            %   Z [-] = Ionization
            %   A [-] = Atomic weight, ion mass / proton mass
            %   gamma [-] = Adiabatic const.
            out = 9.79e5 .* (gamma .* Z .* Te / A).^0.5; % [cm/s]
        end
        function out = getAlfvenSpeed(~,B,A,ni) % Alfven speed
            % Inputs:
            %   B [Gauss] = Magnetic field
            %   A [-] = Atomic weight, ion mass / proton mass
            %   ni [cm^-3] = Ion number denisty
            out = 2.18e11 * A.^(-0.5) .* ni.^(-0.5) .* B; % [cm/s]
        end
        % Misc. / Transport
        function out = getBohmDiffConst(~,T,B) % Bohm diffusion const.
            % Inputs:
            %   T [eV] = Temperature
            %   B [gauss] = Magnetic field
            out = 6.25e6 .* T .* (1./B); % [cm^2/s]
        end
        function out = getTransSpitzerRes(obj,Z,T,n) % eta_perp, Transverse spitzer resitivity
            % Inputs:
            %   Z [-] = Ionization
            %   T [eV] = Temperature
            %   n [cm^-3] = electron number density
            lnA = obj.getCoulombLog(n,T); 
            fprintf('lnA = %2.1f\n',lnA);
            out = 1.03e-2 .* Z .* lnA .* T.^(-3/2); % [Ohm-cm]
        end
        function out = getParSpitzerRes(obj,Z,T,n) % eta_par, Parallel spitzer resitivity
            % Inputs:
            %   Z [-] = Ionization
            %   T [eV] = Temperature
            %   n [cm^-3] = electron number density
            out = 0.5 * obj.getTransSpitzerRes(Z,T,n);
        end
        function out = getParIonViscosity(obj,ni,Ti,Te,Z,A) % Parallel viscosity [SI output]
             % Inputs:
            %   Z [-] = Ionization
            %   Ti [eV] = Ion Temperature
            %   ni [cm^-3] = Ion number density
            
            tau_i = obj.getIonIonCollTime(ni,Z,A,Ti,Te); % s
            % convert to SI
            ni = cgs2SI(ni,'Density'); % m^-3
            load physicalConstants-SI.mat E0
            Ti = Ti * E0; % [J]
            
            out = 0.96 * ni * Ti * tau_i; % [kg/m^3 * m^2/s]
            
        end
        function out = getThermalConductivity(obj,ne,T)
            tau_ee =  obj.getElElCollTime(ne,T); % [s], electron collison time
            load physicalConstants-SI.mat kb me T_eV
            kappa = 3.2 * cgs2SI(ne,'Density') * kb^2 * (T * T_eV) * tau_ee / me;
            out = kappa; % W / m. K
        end
        function [kappa_par,kappa_perp] = getIonConductivities(obj,ni,Ti,Te,A,Z,B)
            % Inputs:
            %   B [Gauss] = Magnetic field
            %   Z [-] = Ionization
            %   A [-] = Atomic weight, ion mass / proton mass
            %   ni [cm^-3] = Ion density
            %  Ti [eV] = Ion temp.
            
            load physicalConstants-SI.mat kb me T_eV mp
            mi = A * mp; % kg
            tau_i =  obj.getIonIonCollTime(ni,Z,A,Ti,Te); % [s], ion collison time
            
            % Parallel
            kappa_par = kb * 3.9 * cgs2SI(ni,'Density') * kb * (Ti * T_eV) * tau_i / mi; % [W m^-1 K^-1]
            
           % Perpendicular
            w_ci = obj.getIonGyroFreq(B,Z,A); % per second
            kappa_perp = kb * 2 * cgs2SI(ni,'Density') * kb * (Ti * T_eV) * tau_i / (mi * w_ci^2 * tau_i^2); % [W m^-1 K^-1]
        end
        function [kappa_par,kappa_perp] = getElectronConductivities(obj,ne,Te,B)
            % Inputs:
            %   B [Gauss] = Magnetic field
            %   ne [cm^-3] = electron density
            %  Te [eV] = elec. temp.
            
            load physicalConstants-SI.mat kb me T_eV 
            tau_e =  obj.getElElCollTime(ne,Te); % [s], electron collison time
            
            % Parallel
            kappa_par = kb * 3.2 * cgs2SI(ne,'Density') * kb * (Te * T_eV) * tau_e / me; % [W m^-1 K^-1]
            
           % Perpendicular
            w_ce = obj.getElectronGyroFreq(B); % per second
            kappa_perp = kb * 4.7 * cgs2SI(ne,'Density') * kb * (Te * T_eV) * tau_e / (me * w_ce^2 * tau_e^2); % [W m^-1 K^-1]
        end
        function out = getBremsPowerDensity(obj,ne,ni,Z,Te) % [W/cm^3]
            % ne = electron denisty [cm^-3]
            % ni = ion denisty [cm^-3]
            % Z = ionization [-]
            % Te = electron temp. [eV]
           out =  1.69e-32 .* ne .* Te .^ (1/2) .* Z.^2 .* ni;
        end
        function out = getCoulombLog(obj,n,T) % Coulomb logarithm
            % n [cm^-3] = electron number density
            % T [eV] = temperature
            lamD = obj.getDebyeLen(n,T); % Debye length, [cm]
            out = log(12 * pi .* n .* lamD.^3); 
        end
        function out = getGamma(~,ni,Te)
            % ni = ion density [cm^-3]
            % Te = electron temp. [eV]
            p = 1.6e-12 * ni .* Te .* (1 + 0.63 .* sqrt(Te) - 2.76e-8 .* ni.^(1/3));
            re = 1.6e-12 .* ni .* ( 1.43 .* sqrt(Te) + 4.2 .* Te + Te.^1.5 .*...
                (1.3 - 0.315.* log(ni * 1e-23 ./ Te)));
            out = 1 + p ./ re;
        end % adiabatic index for ionizing plasma
        function out = getPlasmaBeta(obj,B,ni,Te,A,Z)
            gam = obj.getGamma(ni,Te);
            out = 2 / gam * (obj.getIonSoundSpeed(Te,Z,A,gam) ./ ...
                obj.getAlfvenSpeed(B,A,ni)).^2;
        end
        function [Zbar,ne] = getZbarTF(~,z,a,rho,tev)
            % Approximation to the LTE TF model
            % z = nuclear charge [-]
            % a = atomic mass [-]
            % rho = mass denisty [g/cc]
            % tev = electron temperaure [eV]

            load physicalConstants-cgs.mat mp
            M = a * mp; 
            ni = rho / M; % Ion number denisty (/cc)

            tzero = tev./(z.^(4/3));
            rho = min(rho,10/1e3); 
            r = rho./(z * a); 
            tf = tzero./(1+tzero); 

            % Constants
            a1=0.003323467;                                              
            a2=0.97183224;
            a3=9.26148e-05;
            a4=3.1016524;
            b0=-1.762999;
            b1=1.4317567;
            b2=0.31546338;
            c1=-0.36666667;
            c2=0.98333333;
            alpha=14.3139316;
            beta=0.66240046;

            % caluclate 

            x1=log(tzero);
            x2=a4*x1;
            x1=a2*x1;
            x1=exp(x1);
            x2=exp(x2);

            aa=a1.*x1+a3.*x2;
            arg=b0+b1.*tf+b2.*tf.^7;                                          
            b=-exp(arg);                                            
            c=c2+(c1.*tf);

            x1=log(r);
            x2=c.*x1;
            x1=b.*x1;
            x1=exp(x1);
            q1=aa.*x1;
            x1=c.*log(q1);

            q=exp(x2)+exp(x1);
            cm=1.0./c;                                                                    
            x1=cm.*log(q);

            q=exp(x1);

            x1=log(q);
            x2=beta.*x1;

            x=alpha.*exp(x2);

            %

            f= x./(1.0+x+(sqrt(1.0+(2.0*x))));
            Zbar = f .* z; ne = Zbar .* ni; 
            out = Zbar;
        end
        function [Zbar,ne] = getZbarSpK_Al(~,rho,tev)
            % Approximation to the LTE TF model
            % rho = mass denisty [g/cc]
            % tev = electron temperaure [eV]

            z = 13;
            a = 17;

            load physicalConstants-cgs.mat mp
            M = a * mp; 
            ni = rho / M; % Ion number denisty (/cc)

            fname = 'zb_Al.txt'; 
            data = dlmread(fname,' ',2,0); 
            T = data(:,1); % Temperature [eV]
            Ni = [1e17,5e17,1e18,5e18,1e19]; % ion density [/cc]

            T_grid = repmat(T,1,numel(Ni));
            Ni_grid = repmat(Ni,numel(T),1);

            [T_grid2, Ni_grid2] = meshgrid(T,Ni);

            Zbar = interp2(T_grid2,Ni_grid2,data(:,2:end)',tev,ni); 
            ne = ni * Zbar;
        end
        function out = getOwavencrit(~,lambda) % get O-wave cut-off denisty
            % lamda - Laser wavelength [cm]
            % out = 10^21 / lambda.^2;
            load physicalConstants-cgs me e c
            w = 2 * pi * c / lambda;
            out = w ^ 2 * me / (4 * pi * e^2);
        end
        % Kinetic theory
        function out = getMaxwellDist1D(~,v,vth) % 1D Normalized Maxwellian Distribution
            % Inputs:
            %   v [cm/s] = velocity
            %   vth [cm/s] = sqrt(2 * kb * T / m) = thermal velocity
            out = (pi * vth^2)^(-1/2) .* exp(-v.^2./vth^2); 
        end
        function out = getMaxwellDist(~,v,vth) % 3D Normalized Maxwellian Distribution
            % Inputs:
            %   v [cm/s] = velocity
            %   vth [cm/s] = sqrt(2 * kb * T / m) = thermal velocity
            out = (pi * vth^2)^(-3/2) .* exp(-v.^2./vth^2); 
        end
        % MHD Theory
        function out = collisionalityCriterion(obj,ni,Z,A,Ti,Te,a) % Check if plasma is highly collisional
            % Inputs:
            %   ni [cm^-3] = Ion numbe density
            %   Z [-] = Ionization
            %   A [-] = Atomic weight, ion mass / proton mass
            %   Ti, Te [eV] = Ion, Electron temperature
            %   a [cm] = Chaacteristic plasma length
            load physicalConstants-cgs mp me
            out = sqrt(A * mp / me) .* (1./a) .* obj.getIonIonMfp(ni,Z,A,Ti,Te); % << 1 for highly collisonal
        end
        function out = getPlankianDistribution(obj,T,lam)
            % T = Temperature [K]
            % lam = wavelength [m]
            load physicalConstants-SI.mat c h kb
            out = 2 * h * c^2 ./ lam.^5 .* 1 ./ (exp(h*c./(lam*kb*T)) - 1);
        end
        % Radiation
        function out = getIonizationEnergy(obj,Z,n,l)
            del = 0.75 * l.^(-5);
            Eh = 13.6; % eV
            out =  Z.^2 * Eh ./ (n - del).^2; 
        end
        function out = getRecombBremsPowerDensity(obj,ne,ni,Z,Te,delE) % W/cm^3
            out = 1.69e-32 .* ne .* Te.^(1/2) .* Z.^2 .* ni .* (1 + delE./Te); 
        end
        function out = getInverseBremsAbsorptivity(obj,ne,T,Z,w) % cm^-1
            % ne = electron denisty cm^-3
            % T = electron temp. eV
            % Z = ionization
            % w = [rad/s] angular freq.
%             load physicalConstants-cgs.mat e me E0 hbar
            wp = obj.getElecPlasmaFreq(ne); % rad/s
%             V = max(wp,w) .* max(Z.*e^2/(E0*T),hbar/(me*T*E0)^(1/2));
%             A = obj.getElecThermalVel(T) / V; 
            out = 3.1e-7 .* Z .* ne.^2 * T.^(-3/2) .* ...
                w.^-2 .* (1 - wp.^2/w.^2).^-0.5 .* obj.getCoulombLog(ne,T);
        end
         function out = getSpkEmi(obj,rho_q,Teq) % W/m^3
            % rho_q = mass DENISTY [g per CC]
            %  Teq = electron temp. [eV]
            inDir = checkDir('/Users/rishabhdatta/Dropbox (MIT)/PUFFIN/Codes/Other/atomic_rad_data/');
            rho = dlmread([inDir 'rho.csv'],","); % kg/m^3
            Te = dlmread([inDir 'Te.csv'],","); % eV
            emi = dlmread([inDir 'total_emissivity.csv'],","); % W/m^3

            rho_q = cgs2SI(rho_q,'Mass Density'); % kg/m^3
            out = interp2(Te,rho,emi,Teq,rho_q);% W/m^3
        end
        % Fast MHD Shock Relations
        function out = FastShock_get_r(obj,Ms,gamma,beta) % dimensionless
            % Density, Magnetic field compression ratio r [-]
            % Ms = Sonic Mach no. 
            % gamma = adibatic index
            % beta = Plasma beta
            a = 2 * (2 - gamma);
            b = 2 .* gamma .* (beta + 1) + beta .* gamma .* (gamma -1) .* Ms.^2;
            c = -beta .* gamma .* (gamma+1).* Ms.^2; 
            out1 = 1/(2.*a) .* (-b + sqrt(b.^2 - 4 * a.*c));
            out2 = 1/(2.*a) .* (-b - sqrt(b.^2 - 4 * a.*c));
            out = max([out1,out2]);
        end
        function out = FastShock_get_R(obj,Ms,gamma,beta) % dimensionless
            % Pressure compression ratio r [-]
            % Ms = Sonic Mach no. 
            % gamma = adibatic index
            % beta = Plasma beta
            r = obj.FastShock_get_r(Ms,gamma,beta); 
            R = 1 + (1-r^2)/beta + gamma * Ms^2 * (r -1) / r;
            out = R;
        end
        function out = FastShock_get_M2(obj,Ms,gamma,beta) % dimensionless
            % Pressure compression ratio r [-]
            % Ms = Sonic Mach no. 
            % gamma = adibatic index
            % beta = Plasma beta
            r = obj.FastShock_get_r(Ms,gamma,beta); 
            R = obj.FastShock_get_R(Ms,gamma,beta); 
            out = sqrt(Ms^2 / (r * R));
        end
       
    end
end