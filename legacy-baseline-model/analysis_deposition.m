clc
clear all
close all

folder = ('G:\My Drive\PhD\deposition\simulation\results\simulation\SiO2 coated TiO2\1mM');

fullMatFileName = fullfile(folder,'nststat10_N4e5_2nmSiO2coated_1mM.mat');

x=importdata(fullMatFileName);
Channel_Length =4e-04;
height=10e-06;                          %channel height in m

Q = x.Q;
phi = x.phi;
vp=x.vp;
nstatmax= x.nstatmax;
U_r = x.U_radius;

kB=1.3806488e-23;                       %Boltzmann constant
T=300;             

mu=1.09e-03;

r = x.r;
D=x.D;

for i = 1:length(x.TimeStore)

    timestore = x.TimeStore{i};
    xstore = x.XStore{i};
    ystore = x.YStore{i};
    zstore = x.ZStore{i};
    
    xstore(xstore<0)=0;
    
    tstarprime=mean(timestore);             %average adsorption time of one particle
    stdtstarprime=std(timestore);
    
    [nhistr,xhist]=hist(xstore,100);
    nhistr=nhistr/nstatmax;
    
    hold all
    
    deltaxhist=xhist(3)-xhist(2);
    pdf=sum(nhistr)*deltaxhist*0.5/Channel_Length^0.5./xhist.^0.5;
    
    [nhist,thist]= hist(timestore,100);
    
    deltat=thist(3)-thist(2);
    N_adsorbed_t=cumtrapz(thist,nhist/deltat/nstatmax);
    
    collector_eff_S(i) = mean(vp(i)/(phi(i)*Q)*(N_adsorbed_t/thist)/(r/(height/2-r)));
    collector_eff_S_t  = (vp(i)/(phi(i)*Q)*(N_adsorbed_t./thist)/(r/(height/2-r)));
    
    scatter(thist, collector_eff_S_t,'filled')
    hold on
    xlabel('time, s')
    ylabel('collector efficiency')
    title('Collector efficiency with time')
    set(gca,'fontname','Cambria','FontWeight','bold','fontsize',12);
    
    [r(i)   collector_eff_S(i)]
end
legend('2000 nm','1000 nm','500 nm','200 nm','100 nm')

figure
scatter(2*r*1e9,collector_eff_S,'filled')
xlabel('Particle Size, nm')
ylabel('collector efficiency')
title('Collector efficiency with change in Particle Size')
set(gca,'fontname','Cambria','FontWeight','bold','fontsize',12);
