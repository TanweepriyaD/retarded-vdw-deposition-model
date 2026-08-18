clc
clear all
close all

kB = 1.3806488e-23;        %Boltzmann constant
T = 293.15;                   %Temperature in K

r = 1e-7;
salt = 0.1;
c = salt*1000;              %conversion M to mM

Er = 70;                    %dielectric water at high salt concentration
Eo = 8.854e-12;              %F/m permittivity of free space
e = 1.602e-19;                     %Coulomb electron charge
N_A = 6.02e23;                     %Avogadro's number

lambdad = sqrt((Er*Eo*kB*T)/(2*N_A*c*(e^2))); %Debye length
inv_Dl = 1./lambdad;

zp = -0.04;                                                 %in V
zs = -0.0108;  
C=4*pi*Eo*Er*zp*zs*r/lambdad;

vdata =  xlsread('NvdWE_AirBalloon_100mM.xlsx','100mM');
L = vdata(:,1);
vdwr = vdata(:,2);

min_dist= 5e-10;
max_dist = 500e-9;
dist_y = (min_dist:0.005e-9:max_dist);

for i= 1:length(dist_y)
    E_vdw_y(i) = (interp1(L,vdwr,dist_y(i)))*r/(kB*T);
    E_edl_y (i)= 16*Eo*Er*r*kB*T/e^2*(tanh(e*zp/(4*kB*T)))*(tanh(e*zs/(4*kB*T)))*exp(-dist_y(i)/lambdad);
    E_total(i)= E_vdw_y(i)+ E_edl_y (i);
end
scatter(dist_y*1e09,E_edl_y,'o','filled')
hold on
scatter(dist_y*1e09,E_vdw_y,'o','filled')
scatter(dist_y*1e09,E_total,'o','filled')
hold on
xlabel('Separation Distance (nm)')
ylabel('Energy (KT)')
set(gca,'fontname','Times New Roman','FontWeight','bold','fontsize',26)

legend('EDL','vdW','Total DLVO')
