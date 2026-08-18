clc
clear all
close all

kB=1.3806488e-23;        %Boltzmann constant
T=293.15;                   %Temperature in K

r = 5e-7 ;
salt = 1;
c=salt*1000;                %conversion M to mM

Er = 70;                    %dielectric water at high salt concentration
Eo = 8.854e-12;              %F/m permittivity of free space
e=1.602e-19;                     %Coulomb electron charge
N_A=6.02e23;                     %Avogadro's number

lambdad=sqrt((Er*Eo*kB*T)/(2*N_A*c*(e^2))); %Debye length
inv_Dl = 1./lambdad;

zp = -0.001;                                                 %in V
zs = -0.001;  
C=4*pi*Eo*Er*zp*zs*r/lambdad;

vdata =  xlsread('NvdW_TiO2_1000mM.xlsx','1000mM');
L =vdata(:,1);
vdwr =vdata(:,2);

min_dist= 5e-10;
max_dist = 50e-9;
dist_y = (min_dist:0.005e-9:max_dist);

for i= 1:length(dist_y)
        F_vdw_y(i) = (interp1(L,vdwr,dist_y(i)))*r;
    F_edl_y_DA (i) = C*exp(-dist_y(i)/lambdad);          %solution to screened Poisson equation
    
    F_total(i)= F_vdw_y(i)+ F_edl_y_DA (i);
end

figure
scatter(dist_y*1e9,F_edl_y_DA*1e9,'filled')
hold on
scatter(dist_y*1e9,F_vdw_y*1e9,'filled')
hold on

scatter(dist_y*1e09,F_total*1e9,'.')
xlabel('Separation Distance, nm')
ylabel('Force, nN')
set(gca,'fontname','Cambria','FontWeight','bold','fontsize',12)

legend('EDL','vdW','Total DLVO')
