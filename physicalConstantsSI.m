% Physical Costants - SI
%%%%%%%%%%%%%%%%%%%%%%%%
% Returns value of physical constants in SI
% Example Usage:
%
% load physicalConstants-SI.mat % Load all variables
% load physicalConstants-SI.mat me c % Load electron mass and light speed
%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright 2021 Rishabh Datta, MIT
%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear;

kb = 1.3807e-23; % [J/K], Boltzmann Constant
e = 1.6022e-19; % [C], Elementary charge
me = 9.1094e-31; % [kg], electron mass
mp = 1.6726e-27; % [kg], proton mass
G = 6.6726e-11; % [m^3s^-2kg^-1], Gravitational const.
h = 6.6261e-34; % [Js], Plank const.
hbar = h / (2 * pi); % [Js], Dirac const. 
c = 2.9979e8; % [m/s], Light speed in vacuum
ep0 = 8.8542e-12; % [F/m], Permittivity of free space
mu0 = 4 * pi * 1e-7; % [H/m], Premeability of free space
pe_mass_rat = mp / me; % [-], Proton/electron mass ratio
e_char_mass_rat = e / me; % [C/kg], Electron charge/mass ratio
R_inf = 1.0974e7; % [m^-1], Rydberg const.
a0 = 5.2918e-11; % [m], Bohr radius
at_xsec = pi * a0 ^ 2; % [m^2], Atomic x-section
re = 2.8179e-15; % [m], Classic electron radius
T_xsec = (8 * pi / 3) * re^2; % [m^2], Thompson x-section
lam_cp = hbar / (me * c); % [m], compton wavelength, [m]
alph_fine = e^2 / (2 * ep0 * h * c); % [-], fine structure constant
c1 = 2 * pi * h * c^2; % [Wm^2], First radiation const. 
c2 = h * c / kb; % [mK], 2nd radiation const.
sig_sb = 5.6705e-8; % [Wm^-2K^-4], Stefan-Boltzmann const.
lam0 = 1.2398e-6; % [m], Wavelength associated w/ 1eV
nu0 = 2.4180e14; % [Hz], Freq. associated w/ 1eV
k0 = 8.0655e5; % [m^-1], Wave no. associated w/ 1eV
E0 = 1.6022e-19; % [J], Energy associated w/ 1eV
E_ryd = 13.606; % [eV], Energy associated w/ 1 Ryd
E_kel = 8.6174e-5; % [eV], Energy associated w/ 1 Kelvin
T_eV = 1.1604e4; % [K], Temp. associated w/ 1 eV
N_A = 6.0221e23; % [mol^-1], Avogadro const.
F = 9.6485e4; % [C/mol], Faraday const. 
R = 8.3145; % [J/k/mol], Gas cont.
n0 = 2.6868e25; %[m^-3], Loschmidt no. 
amu = 1.6605e-27; % [kg], atomic mass unit
T0 = 273.15; % [K], Standard temp.
p0 = 1.0133e5; % [Pa], atmospheric pressure
p_Hg = 1.3332e2; % [Pa], 1 mm Hg
V0 = 2.2414e-2; % [m^3], molar volume at STP
Mair = 2.8971e-2; %[kg/mol], molar wirght of air
cal = 4.1868; % [J], calorie
g = 9.8067; % [m/s^2], gravitational acceleration
R0 = 376.73; % [Ohm], Resistance of free space


save('physicalConstants-SI.mat'); 