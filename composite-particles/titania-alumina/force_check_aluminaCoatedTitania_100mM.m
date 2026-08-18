clc
clear all
close all

kB = 1.3806488e-23;        %Boltzmann constant
T = 293.15;                %Temperature in K

r = 1e-7 ;
salt =0.1;
c=salt*1000;                %conversion M to mM

Er = 70;                    %dielectric water at high salt concentration
Eo = 8.854e-12;              %F/m permittivity of free space
e = 1.602e-19;                     %Coulomb electron charge
N_A = 6.02e23;                     %Avogadro's number

lambdad = sqrt((Er*Eo*kB*T)/(2*N_A*c*(e^2))); %Debye length
inv_Dl = 1/lambdad;

zp = 0.04;                                                 %in V
zs = -0.0108;  
C = 4*pi*Eo*Er*zp*zs*r/lambdad;

vdata =  xlsread('NvdW_Al2O3coatedTiO2_100mM.xlsx','100mM');
L = vdata(:,1);
vdwr = vdata(:,2);

min_dist = 5e-10;
max_dist = 50e-9;
dist_y = (min_dist:0.005e-9:max_dist);

for i= 1:length(dist_y)
    if (dist_y(i)<=0.5e-9)
        F_vdw_y(i) = -2.47e-22*exp(-2.084*log(dist_y(i)))*r;
    elseif (0.5e-9<dist_y(i))&&(dist_y(i)<=0.7e-9)
        F_vdw_y(i) = -2.052e-22*exp(-2.093*log(dist_y(i)))*r;
    elseif (0.7e-9<dist_y(i))&&(dist_y(i)<=1e-9)
        F_vdw_y(i) = -1.843e-22*exp(-2.098*log(dist_y(i)))*r;
    elseif (1e-9<dist_y(i))&&(dist_y(i)<=20e-9)
        F_vdw_y(i) = (interp1(L,vdwr,dist_y(i)))*r;
    elseif(20e-9<dist_y(i))&&(dist_y(i)<=25e-9)
        F_vdw_y(i) = -4.219e-25*exp(-2.44*log(dist_y(i)))*r;
    elseif (25e-9<dist_y(i))&&(dist_y(i)<=40e-9)
        F_vdw_y(i) = -6.901e-26*exp(-2.544*log(dist_y(i)))*r;
    elseif (40e-9<dist_y(i))&&(dist_y(i)<=50e-9)
        F_vdw_y(i) = -3.296e-07*exp(-0.184*log(dist_y(i)))*r;
    elseif (50e-9<dist_y(i))&&(dist_y(i)<=70e-9)
        F_vdw_y(i) = -5.104e-27*exp(-2.697*log(dist_y(i)))*r;
    elseif (70e-9<dist_y(i))&&(dist_y(i)<=102.1e-9)
        F_vdw_y(i) = -2.108e-27*exp(-2.751*log(dist_y(i)))*r;
    elseif (102.1e-9<dist_y(i))&&(dist_y(i)<=154.5e-9)
        F_vdw_y(i) = -1.028e-27*exp(-2.795*log(dist_y(i)))*r;
    elseif (154.5e-9<dist_y(i))&&(dist_y(i)<=250.6e-9)
        F_vdw_y(i) = -4.127e-27*exp(-2.854*log(dist_y(i)))*r;
    elseif (250.6e-9<dist_y(i))&&(dist_y(i)<=397.2e-9)
        F_vdw_y(i) = -8.009e-29*exp(-2.961*log(dist_y(i)))*r;
    elseif (397.2e-9<dist_y(i))
        F_vdw_y(i) = -8.818e-30*exp(-3.111*log(dist_y(i)))*r;
    end
end

figure
semilogx(dist_y*1e09,F_vdw_y)

xcheck = log(dist_y);
