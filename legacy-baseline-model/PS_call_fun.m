function[Q,timestore,xstore,ystore,zstore,U_r,adsorption_particle] = PS_call_fun(r,D)

global Q kB T mu  height width Length roughness nstatmax dt N phi vp c
global N_A e epsilon0 epsilon zetap zetaw v_a v_b

debye_length=sqrt((epsilon*epsilon0*kB*T)/(2*N_A*c*(e^2)));     %Debye length

C = 4*pi*epsilon*epsilon0*zetaw*zetap*r/debye_length;

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
        
        n_particle = n_particle+1;
        x=zeros(1,1);
        z=x;
        y=x;
        
        z(1)=(height/2-r-epsz)*rand+r+epsz;
        y(1)=(width/2-r-epsy)*rand+r+epsy;
        x(1)=0;
        
        controle=[controle,z(1)];
        time=(n_particle-1)*vp/Q./phi;
        
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
            
            if (dist_z<1e-9)
                F_vdw_z = -1.54e-22*exp(-2.064*log(dist_z))*r;
            elseif (1e-9<dist_z)&&(dist_z<2.51e-9)
                F_vdw_z = -1.156e-22*exp(-2.078*log(dist_z))*r;
            elseif (2.51e-9<dist_z)&&(dist_z<5e-9)
                F_vdw_z=(3.762e-23*exp(-2.315*log(dist_z)))*r;
            elseif (5e-9<dist_z)&&(dist_z<11e-9)
                F_vdw_z = (-4.498e-24*exp(-2.242*log(dist_z)))*r;
            elseif(11e-9<dist_z)&&(dist_z<12.3e-9)
                F_vdw_z = -4.14e-24*exp(-2.251*log(dist_z))*r;
            elseif (12.3e-9<dist_z)&&(dist_z<30e-9)
                F_vdw_z = -2.989e-26*exp(-2.521*log(dist_z))*r;
            elseif(30e-9<dist_z)&&(dist_z<50e-9)
                F_vdw_z = -4.597e-28*exp(-2.761*log(dist_z))*r;
            elseif(50e-9<dist_z)&&(dist_z<100e-9)
                F_vdw_z = -4.882e-29*exp(-2.895*log(dist_z))*r;
            else
                F_vdw_z = 0;
            end
            
            if (dist_y<1e-9)
                F_vdw_y = -1.54e-22*exp(-2.064*log(dist_y))*r;
            elseif (1e-9<dist_y)&&(dist_y<2.51e-9)
                F_vdw_y = -1.156e-22*exp(-2.078*log(dist_y))*r;
            elseif (2.51e-9<dist_y)&&(dist_y<5e-9)
                F_vdw_y=(3.762e-23*exp(-2.315*log(dist_y)))*r;
            elseif (5e-9<dist_y)&&(dist_y<11e-9)
                F_vdw_y = (-4.498e-24*exp(-2.242*log(dist_y)))*r;
            elseif(11e-9<dist_y)&&(dist_y<12.3e-9)
                F_vdw_y = -4.14e-24*exp(-2.251*log(dist_y))*r;
            elseif (12.3e-9<dist_y)&&(dist_y<30e-9)
                F_vdw_y = -2.989e-26*exp(-2.521*log(dist_y))*r;
            elseif(30e-9<dist_y)&&(dist_y<50e-9)
                F_vdw_y = -4.597e-28*exp(-2.761*log(dist_y))*r;
            elseif(50e-9<dist_y)&&(dist_y<100e-9)
                F_vdw_y = -4.882e-29*exp(-2.895*log(dist_y))*r;
            else
                F_vdw_y = 0;
            end
            
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
            
            if (z(i)<r+epsz)
                adsorption_flag = 1;
                adsorption_particle = adsorption_particle+1;
                xstore=[xstore,x(i)];
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
