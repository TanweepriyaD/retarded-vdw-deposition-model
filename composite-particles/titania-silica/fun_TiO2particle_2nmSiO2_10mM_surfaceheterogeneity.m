function[timestore,xstore,ystore,zstore,U_r,adsorption_particle] = fun_TiO2particle_2nmSiO2_10mM_surfaceheterogeneity(r,phi,D,vp)

global Q kB T mu  height width Length roughness nstatmax dt N c
global N_A e epsilon0 epsilon L vdwr zetap zetaw

debye_length = sqrt((epsilon*epsilon0*kB*T)/(2*N_A*c*(e^2)));     %Debye length

facteur=0;
for n=1:10
    betan=(2*n-1)*pi/width;
    facteurn=1/betan^4*(1-2/betan/height*tanh(betan*height/2));
    facteur=facteur+facteurn;
end

dp_by_dl =Q*(mu*width/8/height)/facteur;

U_r=0;
for n=1:10
    betan=(2*n-1)*pi/width;
    U_r_n=2/width/betan^4*(1-cosh(betan*(r-height/2))/cosh(betan*height/2));
    U_r=U_r+U_r_n;
end
U_r=4*dp_by_dl/mu/width*U_r; %fluid velocity near wall

diff_length=sqrt(2*D*dt); %Brownian motion particles

epsz=0.5e-09;   % depletion near the wall at the entrance = Born layer
epsy=epsz;

xstore=[];
timestore=[];
controle=[];

x=zeros(1,1);
z=x;
y=x;

adsorption_particle = 0;
exit_particle =0;

zstore=[];
ystore=[];

time=0;

for nstat=1:nstatmax
    n_particle = 0;
    nstat;
    while (n_particle < N)
        z_check = [];
        n_particle = n_particle+1;
        x=zeros(1,1);
        z=x;
        y=x;
        
        z(1)=(height/2-r-epsz)*rand+r+epsz;
        y(1)=(width/2-r-epsy)*rand+r+epsy;
        x(1)=0;
        
        controle=[controle,z(1)];
        time=(n_particle-1)*vp/Q/phi;
        
        i=1;
        exit_flag=0;
        adsorption_flag=0;
        
        while (adsorption_flag==0)&(exit_flag==0)
            i=i+1;
            
            eta=z(i-1)/r;
            
            beta_y = (6*y(i-1)^2-10*y(i-1)*r+4*r^2)/(6*y(i-1)^2-3*y(i-1)*r-r^2);
            dbeta_y = ((6*y(i-1)^2-3*y(i-1)*r-r^2)*(12*y(i-1)-10*r)-(6*y(i-1)^2-10*y(i-1)*r+4*r^2)*(12*y(i-1)-3*r)) /(6*y(i-1)^2-3*y(i-1)*r-r^2)^2;
            
            beta_z = (6*z(i-1)^2-10*z(i-1)*r+4*r^2)/(6*z(i-1)^2-3*z(i-1)*r-r^2);
            dbeta_z = (((6*z(i-1)^2-3*z(i-1)*r-r^2)*(12*z(i-1)-10*r))-((6*z(i-1)^2-10*z(i-1)*r+4*r^2)*(12*z(i-1)-3*r))) /(6*z(i-1)^2-3*z(i-1)*r-r^2)^2;
            
            beta_x = 1 - 9/16/z(i-1)*r + 1/8/(z(i-1))^3*r^3 - 45/256/(z(i-1))^4*r^4 - 1/16/(z(i-1))^5*r^5; % diffusion parallel to the wall
            beta_x = beta_x *(1 - 9/16/(y(i-1)/r) + 1/8/((y(i-1)/r)^3) - 45/256/((y(i-1)/r)^4) - 1/16/((y(i-1)/r)^5)); % diffusion parallel to the wall
            
            if(eta>=2)
                gamma=1-5/16/eta^3;
            elseif(eta>1+1e-04)
                gamma=1/eta*exp(0.68902+0.54756*log(eta-1)+0.072332*(log(eta-1))^2+0.0037644*(log(eta-1))^3);
            else
                gamma=0.7431/(0.6376-0.2*log(eta-1));
            end
            
            sautx=rand-0.5;
            sauty=rand-0.5;
            sautz=rand-0.5;
            
            dist_y = y(i-1)- r + roughness;
            dist_z = z(i-1)- r + roughness;
            
            if (dist_z<=0.5e-9)
                F_vdw_z = -1.197e-21*exp(-1.993*log(dist_z))*r;
            elseif (0.5e-9<dist_z)&&(dist_z<=0.7e-9)
                F_vdw_z = -1.238e-21*exp(-1.991*log(dist_z))*r;
            elseif (0.7e-9<dist_z)&&(dist_z<=1e-9)
                F_vdw_z = -1.273e-21*exp(-1.99*log(dist_z))*r;
            elseif (1e-9<dist_z)&&(dist_z<=2.6e-9)
                F_vdw_z = -1.028e-21*exp(-2*log(dist_z))*r;
            elseif (2.6e-9<dist_z)&&(dist_z<=5e-9)
                F_vdw_z = -9.493e-22*exp(-2.004*log(dist_z))*r;
            elseif (5e-9<dist_z)&&(dist_z<=12.6e-9)
                F_vdw_z = (interp1(L,vdwr,dist_z))*r;
            elseif (12.6e-9<dist_z)&&(dist_z<=20e-9)
                F_vdw_z = -3.043e-23*exp(-2.193*log(dist_z))*r;
            elseif(20e-9<dist_z)&&(dist_z<=25e-9)
                F_vdw_z = -2.123e-24*exp(-2.343*log(dist_z))*r;
            elseif (25e-9<dist_z)&&(dist_z<=40e-9)
                F_vdw_z = -2.195e-25*exp(-2.473*log(dist_z))*r;
            elseif (40e-9<dist_z)&&(dist_z<=50e-9)
                F_vdw_z = -2.936e-26*exp(-2.591*log(dist_z))*r;
            elseif (50e-9<dist_z)&&(dist_z<=70e-9)
                F_vdw_z = -8.935e-27*exp(-2.661*log(dist_z))*r;
            elseif (70e-9<dist_z)&&(dist_z<=100e-9)
                F_vdw_z = -3.052e-27*exp(-2.727*log(dist_z))*r;
            elseif (100e-9<dist_z)&&(dist_z<=155e-9)
                F_vdw_z = -1.3e-27*exp(-2.78*log(dist_z))*r;
            elseif (155e-9<dist_z)&&(dist_z<=251e-9)
                F_vdw_z = -4.759e-28*exp(-2.844*log(dist_z))*r;
            elseif (251e-9<dist_z)&&(dist_z<=397.2e-9)
                F_vdw_z = -8.749e-29*exp(-2.955*log(dist_z))*r;
            elseif (400e-9<dist_z)
                F_vdw_z = -9.38e-30*exp(-3.106*log(dist_z))*r;
            end
            
            if (dist_y<=0.5e-9)
                F_vdw_y = -1.197e-21*exp(-1.993*log(dist_y))*r;
            elseif (0.5e-9<dist_y)&&(dist_y<=0.7e-9)
                F_vdw_y = -1.238e-21*exp(-1.991*log(dist_y))*r;
            elseif (0.7e-9<dist_y)&&(dist_y<=1e-9)
                F_vdw_y = -1.273e-21*exp(-1.99*log(dist_y))*r;
            elseif (1e-9<dist_y)&&(dist_y<=2.6e-9)
                F_vdw_y = -1.028e-21*exp(-2*log(dist_y))*r;
            elseif (2.6e-9<dist_y)&&(dist_y<=5e-9)
                F_vdw_y = -9.493e-22*exp(-2.004*log(dist_y))*r;
            elseif (5e-9<dist_y)&&(dist_y<=12.6e-9)
                F_vdw_y = (interp1(L,vdwr,dist_y))*r;
            elseif (12.6e-9<dist_y)&&(dist_y<=20e-9)
                F_vdw_y = -3.043e-23*exp(-2.193*log(dist_y))*r;
            elseif(20e-9<dist_y)&&(dist_y<=25e-9)
                F_vdw_y = -2.123e-24*exp(-2.343*log(dist_y))*r;
            elseif (25e-9<dist_y)&&(dist_y<=40e-9)
                F_vdw_y = -2.195e-25*exp(-2.473*log(dist_y))*r;
            elseif (40e-9<dist_y)&&(dist_y<=50e-9)
                F_vdw_y = -2.936e-26*exp(-2.591*log(dist_y))*r;
            elseif (50e-9<dist_y)&&(dist_y<=70e-9)
                F_vdw_y = -8.935e-27*exp(-2.661*log(dist_y))*r;
            elseif (70e-9<dist_y)&&(dist_y<=100e-9)
                F_vdw_y = -3.052e-27*exp(-2.727*log(dist_y))*r;
            elseif (100e-9<dist_y)&&(dist_y<=155e-9)
                F_vdw_y = -1.3e-27*exp(-2.78*log(dist_y))*r;
            elseif (155e-9<dist_y)&&(dist_y<=251e-9)
                F_vdw_y = -4.759e-28*exp(-2.844*log(dist_y))*r;
            elseif (251e-9<dist_y)&&(dist_y<=397.2e-9)
                F_vdw_y = -8.749e-29*exp(-2.955*log(dist_y))*r;
            elseif (400e-9<dist_y)
                F_vdw_y = -9.38e-30*exp(-3.106*log(dist_y))*r;
            end
            
            C = 4*pi*epsilon*epsilon0*zetaw(n_particle)*zetap(n_particle)*r/debye_length;
            F_edl_y = C*exp(-dist_y/debye_length);
            F_edl_z = C*exp(-dist_z/debye_length);
            
            Vfluid=0;
            for n=1:5
                betan = (2*n-1)*pi/width;
                Vfluidn = (2*mod(n,2)-1)/betan^3*(1-cosh(betan*(z(i-1)-height/2))/cosh(betan*height/2))*cos(betan*(y(i-1)-width/2));
                Vfluid = Vfluid+Vfluidn;
            end
            Vfluid = 4*dp_by_dl/mu/width*Vfluid;
            mobility = 1/(6*pi*r*mu);         %particle mobility
            
            x(i)=x(i-1)+gamma*Vfluid*dt+ diff_length*sqrt(beta_x)*sautx/abs(sautx);
            y(i)=y(i-1)+diff_length*sqrt(beta_y)*sauty/abs(sauty)+dbeta_y*D*dt + D/kB/T*beta_y*(F_vdw_y+F_edl_y)*dt;
            z(i)=z(i-1)+diff_length*sqrt(beta_z)*sautz/abs(sautz)+dbeta_z*D*dt + D/kB/T*beta_z*(F_vdw_z+F_edl_z)*dt;
            
            if (z(i)>height/2)
                z(i)=height-z(i);
            end
            if (y(i)>width/2)
                y(i)=width-y(i);
            end
            
            z_check = [z_check, z(i)];
            if (z(i)<r+epsz)
                adsorption_flag = 1;
                adsorption_particle = adsorption_particle+1;
                xstore=[xstore,x(i)];
                ystore = [ystore,y(i)];
                x(i);
                timestore=[timestore,time];
                zstore=[zstore,z(1)];
                controlezstore=zstore-r;
                adsorption_flag=1;
                
            end
            
            if (y(i)<r+epsy)
                adsorption_flag=1;
            end
            
            if (x(i)>Length)
                exit_flag=1;
                exit_particle=exit_particle+1;
            end
        end
    end
end
end
