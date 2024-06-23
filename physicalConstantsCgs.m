% Physical Costants - cgs
%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Returns value of physical constants in cgs
%
% Example Usage:
%
% load physicalConstants-cgs.mat % Load all variables
% load physicalConstants-cgs.mat me c % Load electron mass and light speed
%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2021 Rishabh Datta, MIT
%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear;

kb = 1.3807e-16; % [erg/K], Boltzmann Constant
e = 4.8032e-10; % [statcoul], Elementary charge
me = 9.1094e-28; % [g], electron mass
mp = 1.6726e-24; % [g], proton mass
G = 6.6726e-8; % [dyne-cm^2/g^2], Gravitational const.
h = 6.6261e-27; % [erg-sec], Plank const.
hbar = 1.0546e-27; % [erg-sec], Dirac const. 
c = 2.9979e10; % [cm/s], Light speed in vacuum
pe_mass_rat = 1.8362e3; % [-], Proton/electron mass ratio
e_char_mass_rat = 5.2728e17; % [statcoul/g], Electron charge/mass ratio
R_inf = 1.0974e5; % [cm^-1], Rydberg const.
a0 = 5.2918e-9; % [cm], Bohr radius
at_xsec = 8.7974e-17; % [cm^2], Atomic x-section
re = 2.8179e-13; % [cm], Classic electron radius
T_xsec = 6.6525e-25; % [cm^2], Thompson x-section
lam_cp = 3.38616e-11; % [cm], compton wavelength, [m]
alph_fine = 7.2972e-3; % [-], fine structure constant
c1 = 3.7418e-5; % [erg-cm^2/s], First radiation const. 
c2 = 1.4388; % [cm-K], 2nd radiation const.
sig_sb = 5.6705e-5; % [erg/cm^2-sec-K^4], Stefan-Boltzmann const.
lam0 = 1.2398e-4; % [cm], Wavelength associated w/ 1eV
nu0 = 2.4180e14; % [Hz], Freq. associated w/ 1eV
k0 = 8.0655e3; % [cm^-1], Wave no. associated w/ 1eV
E0 = 1.6022e-12; % [erg], Energy associated w/ 1eV
E_ryd = 13.606; % [eV], Energy associated w/ 1 Ryd
E_kel = 8.6174e-5; % [eV], Energy associated w/ 1 Kelvin
T_eV = 1.1604e4; % [K], Temp. associated w/ 1 eV
N_A = 6.0221e23; % [mol^-1], Avogadro const.
F = 2.8925e14; % [statcoul/mol], Faraday const. 
R = 8.3145e7; % [erg/K-mol], Gas cont.
n0 = 2.6868e19; %[cm^-3], Loschmidt no. 
amu = 1.6605e-24; % [g], atomic mass unit
T0 = 273.15; % [K], Standard temp.
p0 = 1.0133e6; % [dyne/cm^2], atmospheric pressure
p_Hg = 1.3332e3; % [dyne/cm^2], 1 mm Hg
V0 = 2.2414e4; % [cm^3], molar volume at STP
Mair = 28.971; %[g/mol], molar wirght of air
cal = 4.1868e7; % [erg], calorie
g = 980.67; % [cm/s^2], gravitational acceleration
R0 = 376.73; % [Ohm], Resistance of free space

save('physicalConstants-cgs.mat'); 