function out = cgs2SI(val,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert from cgs to SI units
% Example usage:
% v = cgs2SI(10,'Velocity'); % Velocity, [cm/s] to [m/s]
% B = cgs2SI(1,'Magnetic Field'); % Magnetic Field, [T] to [gauss]
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2021 Rishabh Datta, MIT
%%%%%%%%%%%%%%%%%%%%%%%%%%
alph = 10^2; % [cm/m]
beta = 1e7; % [erg/J]
gam = 1e3; % [g/kg]
load physicalConstants-SI.mat ep0 mu0 
switch type
    case 'Capacitance'
        k = alph * (4 * pi * ep0)^-1; 
    case 'Charge'
        k = (alph * beta / (4 * pi * ep0))^0.5; 
    case 'Charge Density'
        k = (beta / (4 * pi * alph^5 * ep0))^0.5; 
    case 'Conductance'
        k = alph * (4 * pi * ep0)^-1; 
    case 'Conductivity'
        k = (4 * pi * ep0)^-1;
    case 'Current'
        k = (alph * beta / (4 * pi * ep0))^0.5;
    case 'Current Density'
        k = 1/ alph^2 * (alph * beta / (4 * pi * ep0))^0.5; 
    case 'Density'
        k = alph^-3;     
    case 'Electric Field'
        k = (4 * pi * beta * ep0 / alph^3)^0.5;
    case 'EMF'
        k = alph * (4 * pi * beta * ep0 / alph^3)^0.5;
    case 'Electric Potential'
        k = (4 * pi * beta * ep0 / alph)^0.5;
    case 'Electric Conductivity'
        k = (4 * pi * ep0)^-1; 
    case 'Energy'
        k = beta; 
    case 'Energy Density'
        k = beta / alph^3; 
    case 'Force'
        k = beta / alph; 
    case 'Impedance'
        k = 4 * pi * ep0 / alph;    
    case 'Inductance'
        k = 4 * pi * ep0 / alph; 
    case 'Length'
        k = alph;
    case 'Magnetic Induction'
        k = (4 * pi * beta / (alph^3 * mu0))^0.5; 
    case 'Magnetic Field'
        k = (4 * pi * beta / (alph^3 * mu0))^0.5; 
    case 'Magnetic Flux'
        k = (4 * pi * beta / (alph^3 * mu0))^0.5 * alph^2; 
    case 'Magnetic Intensity'
        k = (4 * pi * mu0 *  beta / alph^3)^0.5; 
    case 'Magnetic Moment'
        k = 1e3;
    case 'Magnetization'
        k = (4 * pi * mu0 *  beta / alph^3)^0.5; 
    case 'Magnetomotance'
        k = alph * (4 * pi * mu0 *  beta / alph^3)^0.5;     
    case 'Mass'
        k = beta / alph^2;
    case 'Mass Density' 
        k = gam * alph^-3;    
    case 'Momentum'
        k = beta / alph;
    case 'Momentum Density'
        k = beta / alph / alph^3;
    case 'Polarization'
        k =  (beta / (4 * pi * alph^5 * ep0))^0.5;    
    case 'Power'
        k = beta;
    case 'Power Density'
        k = beta / alph^3;
    case 'Pressure'
        k = beta / alph^3;
    case 'Reluctance'
        k = 4 * pi * 1e-9;    
    case 'Resistance'
        k = 4 * pi * ep0 / alph;
    case 'Resistivity'
        k = 4 * pi * ep0 / alph * alph;    
    case 'Velocity'
        k = alph;
    otherwise
        error('Type not supported');
end
       out = val ./ k; 
end