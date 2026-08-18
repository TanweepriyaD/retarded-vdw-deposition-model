clc
clear all
close all

x =[1.359 1.93 3.04 9.61];
ln_x = log(x);
x_z = [x;x;x;x];
ln_x_z=log(x_z);
y = [200 500 1000 2000];
ln_y = log(y);
z = [1.221 0.9039 0.4055 0.1113; 1.393 0.9309 0.4228 0.1461; 1.391 0.9719 0.4891 0.135;1.626 1.155 0.5623 0.2009];
ln_z = log(z);

x1 = [50 25 10 1];
ln_x1 = log(x1);
z1 = [1.221 0.9039 0.4055 0.1113];
z2 = [1.393 0.9309 0.4228 0.1461];
z3 = [1.391 0.9719 0.4891 0.135];
z4 = [1.626 1.155 0.5623 0.2009];
ln_z1 = log(z1);
ln_z2 = log(z2);
ln_z3 = log(z3);
ln_z4 = log(z4);

y1(1:4) = 200;
y2(1:4) = 500;
y3(1:4) = 1000;
y4(1:4) = 2000;
ln_y1 = log(y1);
ln_y2 = log(y2);
ln_y3 = log(y3);
ln_y4 = log(y4);

ln_y_z = [ln_y1;ln_y2;ln_y3;ln_y4];

plot (x,z1,'Linewidth',4)
hold on
plot (x,z2,'Linewidth',4)
plot (x,z3,'Linewidth',4)
plot (x,z4,'Linewidth',4)

xlabel('Debye length (nm)')
zlabel('Deposition rate')
legend('200 nm','500 nm','1000 nm','2000 nm')
set(gca,'fontname','Cambria','FontWeight','bold','fontsize',26);

figure
plot (x(1:3),ln_z1(1:3),'-.','Linewidth',4)
hold on
plot (x(1:3),ln_z2(1:3),'-.','Linewidth',4)
plot (x(1:3),ln_z3(1:3),'-.','Linewidth',4)
plot (x(1:3),ln_z4(1:3),'-.','Linewidth',4)

xlabel('Debye length (nm)')
zlabel('Deposition rate')
legend('200 nm','500 nm','1000 nm','2000 nm')
set(gca,'fontname','Cambria','FontWeight','bold','fontsize',26);
