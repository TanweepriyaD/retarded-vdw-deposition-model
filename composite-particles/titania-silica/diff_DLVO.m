clc
clear all
close all

kB=1.3806488e-23;        %Boltzmann constant
T=293.15;                   %Temperature in K

r = [1e-6 5e-7 2.5e-7 1e-7]+(2e-9) ;
salt = 0.05;
c = salt*1000;                %conversion M to mM

Er = 70;                    %dielectric water at high salt concentration
Eo = 8.854e-12;              %F/m permittivity of free space
e = 1.602e-19;                     %Coulomb electron charge
N_A = 6.02e23;                     %Avogadro's number
mu = 1.002e-03;
D = (kB*T/6/pi/mu)./r;
mobility = 1./(6*pi*r*mu);

lambdad = sqrt((Er*Eo*kB*T)/(2*N_A*c*(e^2))); %Debye length
inv_Dl = 1/lambdad;

zp = -0.015;                                                 %in V
zs = -0.015;
C = 4*pi*Eo*Er*zp*zs.*r/lambdad;
u = 1e-02;

vdata =  xlsread('NvdW_SiO2coatedTiO2_50mM.xlsx','50mM');
L = vdata(:,1);
vdwr = vdata(:,2);

min_dist= 5e-10;
max_dist = 50e-9;
dist_y = (min_dist:0.005e-9:max_dist);
for j = 1:length(r)
    for i= 1:length(dist_y)
        
        F_vdw_y(i) = (interp1(L,vdwr,dist_y(i)))*r(j);
        
        F_edl_y_DA (i) = C(j)*exp(-dist_y(i)/lambdad);          %solution to screened Poisson equation
        
        F_total(j,i) = F_vdw_y(i)+ F_edl_y_DA (i);
        
        D_normal(i) = D(j)*(((6*dist_y(i))+(2*r(j)*dist_y(i))/((6*dist_y(i)^2)+(9*r(j)*dist_y(i))+(2*r(j)^2))));
                D_F(j,i) = kB*T/D_normal(i)*u;
        ratio(j,i) = abs(F_total(j,i))/D_F(j,i);
    end
    scatter(dist_y*1e09,ratio(j,:),'o','filled')
    hold on
    xlabel('Separation Distance (nm)')
    set(gca,'fontname','Times New Roman','FontWeight','bold','fontsize',26)
    
end
legend('2000 nm','1000 nm','500 nm','200 nm')

xcheck = log(dist_y);
