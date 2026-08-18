clc
clear all
close all

global r Q kB h T mu rho height width Length roughness salt nstatmax dt N c phi D vp
global L vdwr zetap zetaw 
global N_A e epsilon0 epsilon

vdata =  xlsread('NvdW_SiO2coatedTiO2_10mM.xlsx','10mM');
L = vdata(:,1);
vdwr =vdata(:,2);

%% Change
salt=0.01;                                 % in M (mol/L) M to mM occurs inside the function

c = salt*1000;

zetap = -0.028;                                                  %in V
zetaw = -0.028;                                                  %in V

%% Particles parameters
N = 4e5;                             % no of particles
phi1=3e-3;                               %particle concentration (e.g. 1.0% is 1e-2)

r_core = 1e-6;                        % particle radius, m
r_coating = 2e-9;  
r = r_core + r_coating;                 % particle radius with coating, m
phi = phi1*(r/(1e-6+r_coating))^3;

%% Parameters for EDL

epsilon = 70;                       %dielectric water at high salt concentration
epsilon0 = 8.854e-12;               %F/m permittivity of free space
e = 1.602e-19;                      %Coulomb electron charge
N_A = 6.02e23;                      %Avogadro's number

u = [1e-02 0.278 0.55 0.14];                          %average velocity in m/s
kB = 1.3806488e-23;                 %Boltzmann constant
h  = 1.054e-34 ;                    %Reduced planck constant(H/2pi),J.s

T = 293.15;                         %Temperature in K
mu = 1.002e-03;                            
rho = 1e03;                         %density in kg/m3

Rep = 2*r*rho*u/mu;
D = (kB*T/6/pi/mu)/r;            %Stokes-Einstein diffusion coefficient
vp = 4/3*pi*r^3;                 %volume of spherical particle

%% channel dimension
height = 10e-06;                          %channel height in m
width = 100e-06;                          %channel width in m
Length = 4e-04;                          %channel length in m

Q = u.*height*width;

roughness = 0;                            %wall roughness

%% CALCULATION PARAMETERS
nstatmax = 10;                          % statistical mean on nstatmax processes
dt = 1e-05;                             % temporal increment - typical value dt=0.5e-05

%% Run function

for i = 1:length(Q)
    i
    [timestore,xstore,ystore,zstore,U_r,adsorption_particle] = fun_TiO2particle_2nmSiO2_10mM_velvar(Q(i));
    
    TimeStore{i} = timestore
    XStore{i} = xstore
    YStore{i} = ystore;
    ZStore{i} = zstore;
    U_radius{i} = U_r;
    deposited{i} = adsorption_particle;
    
end

save('G:\My Drive\PhD\deposition\simulation\results\simulation\SiO2 coated TiO2\10mM\nststat10_N4e5_2nmSiO2coated2000nmTiO2_10mM_velvar')
