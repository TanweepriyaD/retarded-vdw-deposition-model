clc
clear all
close all

kB=1.3806488e-23;        %Boltzmann constant
T=293.15;                   %Temperature in K

r = 5e-07+2e-8 ;
salt = 1;
c=salt*1000;                %conversion M to mM

Er = 70;                    %dielectric water at high salt concentration
Eo = 8.854e-12;              %F/m permittivity of free space
e = 1.602e-19;                     %Coulomb electron charge
N_A = 6.02e23;                     %Avogadro's number

lambdad = sqrt((Er*Eo*kB*T)/(2*N_A*c*(e^2))); %Debye length
inv_Dl = 1./lambdad;

zp = -0.045;                                                 %in V
zs = -0.012;  
C=4*pi*Eo*Er*zp*zs*r/lambdad;

vdata =  xlsread('NvdW_PScoatedSilica_75mM.xlsx','75mM');
L = vdata(:,1);
vdwr = vdata(:,2);

min_dist= 5e-10;
max_dist = 100e-9;
dist_y = (min_dist:0.005e-9:max_dist);

for i= 1:length(dist_y)
    if (dist_y(i)<=0.5e-9)
        F_vdw_y(i) = -1.355e-23*exp(-2.2*log(dist_y(i)))*r;
       elseif (0.5e-9<dist_y(i))&&(dist_y(i)<=12.6e-9)   
         F_vdw_y(i) = (interp1(L,vdwr,dist_y(i)))*r;
    elseif (12.6e-9<dist_y(i))&&(dist_y(i)<=20e-9)
        F_vdw_y(i) = -3.565e-26*exp(-2.509*log(dist_y(i)))*r;
    elseif(20e-9<dist_y(i))&&(dist_y(i)<=25e-9)
        F_vdw_y(i) = -2.68e-27*exp(-2.655*log(dist_y(i)))*r;
    elseif (25e-9<dist_y(i))&&(dist_y(i)<=40e-9)
        F_vdw_y(i) = -3.318e-28*exp(-2.774*log(dist_y(i)))*r;
    elseif (40e-9<dist_y(i))&&(dist_y(i)<=50e-9)
        F_vdw_y(i) = -6.471e-29*exp(-2.87*log(dist_y(i)))*r;
    elseif (50e-9<dist_y(i))&&(dist_y(i)<=70e-9)
        F_vdw_y(i) = -3.231e-29*exp(-2.911*log(dist_y(i)))*r;
    elseif (70e-9<dist_y(i))&&(dist_y(i)<=100e-9)
        F_vdw_y(i) = -2.538e-29*exp(-2.926*log(dist_y(i)))*r;
    elseif (100e-9<dist_y(i))&&(dist_y(i)<=155e-9)
        F_vdw_y(i) = -3.22e-29*exp(-2.911*log(dist_y(i)))*r;
    elseif (155e-9<dist_y(i))&&(dist_y(i)<=251e-9)
        F_vdw_y(i) = -3.704e-29*exp(-2.902*log(dist_y(i)))*r;
    elseif (251e-9<dist_y(i))&&(dist_y(i)<=397.2e-9)
        F_vdw_y(i) = -1.363e-29*exp(-2.968*log(dist_y(i)))*r;
    elseif (400e-9<dist_y(i))
        F_vdw_y(i) = -1.699e-30*exp(-3.109*log(dist_y(i)))*r;
    end
    
    F_edl_y_DA (i) = C*exp(-dist_y(i)/lambdad);          %solution to screened Poisson equation
    F_total(i) = F_vdw_y(i)+F_edl_y_DA (i);
end

figure
scatter(dist_y*1e9,F_edl_y_DA*1e9,'filled')
hold on
scatter(dist_y*1e9,F_vdw_y*1e9,'filled')
hold on
scatter(dist_y*1e09,F_total*1e9,'.')
xlabel('Separation Distance (nm)')
ylabel('Force (nN)')
set(gca,'fontname','Cambria','FontWeight','bold','fontsize',12)

legend('EDL','vdW','Total DLVO')

xcheck = log(dist_y);
