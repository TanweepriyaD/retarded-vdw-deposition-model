clc
clear all
close all

global u Q kB h T mu rho height width Length roughness salt nstatmax dt N c
global L vdwr zetap zetaw
global N_A e epsilon0 epsilon

vdata =  xlsread('NvdW_PDAcoatedSiO2_5mM.xlsx','5mM');
L = vdata(:,1);
vdwr  =vdata(:,2);

%% Change
salt=0.005;                                 % in M (mol/L) M to mM occurs inside the function

c = salt*1000;

zetap = -0.0452;                                                 %in V
zetaw = -0.048;                                                  %in V

%% Particles parameters
N = 4e5;                             % no of particles
phi1=7e-4;                               %particle concentration (e.g. 1.0% is 1e-2)

r_core = 1e-6;                        % particle radius, m
r_coating = 5e-9;  
r = r_core + r_coating;                 % particle radius with coating, m
phi = phi1;

%% Parameters for EDL

epsilon = 70;                       %dielectric water at high salt concentration
epsilon0 = 8.854e-12;               %F/m permittivity of free space
e = 1.602e-19;                      %Coulomb electron charge
N_A = 6.02e23;                      %Avogadro's number

u = 5.5e-3;                          %average velocity in m/s
kB = 1.3806488e-23;                 %Boltzmann constant
h  = 1.054e-34 ;                    %Reduced planck constant(H/2pi),J.s

T = 293.15;                         %Temperature in K
mu = 1.002e-03;                            
rho = 1e03;                         %density in kg/m3

Rep = 2*r*rho*u/mu;
D = (kB*T/6/pi/mu)/r;                    %Stokes-Einstein diffusion coefficient
vp = 4/3*pi*r^3;                         %volume of spherical particle

%% channel dimension
height = 10e-06;                          %channel height in m
width = 100e-06;                          %channel width in m
Length = 10e-03;                          %channel length in m

Q = u*height*width;

roughness = 0;                            %wall roughness

%% CALCULATION PARAMETERS
nstatmax = 10;                          % statistical mean on nstatmax processes
dt = 1e-05;                             % temporal increment - typical value dt=0.5e-05

%% Run function

for i = 1:length(r)
    i
    
    [timestore,xstore,ystore,zstore,U_r,adsorption_particle] = fun_SiO2particle_2nmPDA_5mM(r(i),phi(i),D(i),vp(i));
    
    TimeStore{i} = timestore;
    XStore{i} = xstore
    YStore{i} = ystore;
    ZStore{i} = zstore;
    U_radius{i} = U_r;
    deposited{i} = adsorption_particle;
    
end

save('G:\My Drive\PhD\deposition\simulation\results\simulation\PDA coated SiO2\5mM\nststat10_N4e5_2nmPDAcoatedSiO2_5mM_5')
