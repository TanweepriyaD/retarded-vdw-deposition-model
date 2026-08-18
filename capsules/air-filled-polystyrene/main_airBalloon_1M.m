clc
clear all
close all

global u Q kB h T mu rho height width Length roughness salt nstatmax dt r r_coating N c
global L vdwr zetap zetaw
global N_A e epsilon0 epsilon

vdata =  xlsread('NvdW_AirBalloon.xlsx','1M');
L =vdata(:,1);
vdwr =vdata(:,2);

%% Change
salt=1;                                 % in M (mol/L) M to mM occurs inside the function

c= salt*1000;

zetap = -0.008;                                                 %in V
zetaw = -0.001;                                                 %in V

%% Particles parameters
N = 4e5;                             % no of particles
r_core = [7.5e-7 4e-7 2e-7 1.5e-7];    % core radius, m

r_coating(1:4) = 2e-9;                  % coating thickness, m
r = r_core + r_coating;                 % particle radius with coating, m

phi1=3e-3;                               %particle concentration (e.g. 1.0% is 1e-2)
phi = phi1*(r.^3)./r(1)^3;

%% Parameters for EDL

epsilon = 70;                       %dielectric water at high salt concentration
epsilon0 = 8.854e-12;               %F/m permittivity of free space
e = 1.602e-19;                      %Coulomb electron charge
N_A = 6.02e23;                      %Avogadro's number

u = 1e-02;                          %average velocity in m/s
kB = 1.3806488e-23;                 %Boltzmann constant
h  = 1.054e-34 ;                    %Reduced planck constant(H/2pi),J.s

T = 293.15;                         %Temperature in K
mu = 1.002e-03;                            
rho = 1e03;                         %density in kg/m3

Rep = 2.*r*rho*u/mu;
D = (kB*T/6/pi/mu)./r;            %Stokes-Einstein diffusion coefficient
vp = 4/3*pi*r.^3;                 %volume of spherical particle

%% channel dimension
height = 10e-06;                          %channel height in m
width = 100e-06;                          %channel width in m
Length = 4e-04;                          %channel length in m

Q = u*height*width;

roughness = 0;                            %wall roughness

%% CALCULATION PARAMETERS
nstatmax = 10;                          % statistical mean on nstatmax processes
dt = 1e-05;                             % temporal increment - typical value dt=0.5e-05

%% Run function

for i = 1:length(r)
    i
    [timestore,xstore,ystore,zstore,U_r,adsorption_particle] = fun_airBalloon_2nmPS_1M(r(i),phi(i),D(i),vp(i));
    
    TimeStore{i} = timestore
    XStore{i} = xstore;
    YStore{i} = ystore;
    ZStore{i} = zstore;
    U_radius{i} = U_r;
    deposited{i} = adsorption_particle;
    
end

save('G:\My Drive\PhD\deposition\simulation\results\simulation\airBalloon\1000mM\nststat10_N4e5_2nmPScoatingAirBalloon_1M_r_interp')
