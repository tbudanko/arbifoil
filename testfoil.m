close all; clear all;

t1 = linspace(0,pi/2,10);
t2 = linspace(pi/2,pi,10); t2 = t2(2:end); 


zeta1 = (2-0.5*sin(t1).*(pi-t1)).*exp(j*t1);
zeta2 = (2-0.5*sin(t2).*(pi-t2)).*exp(j*t2);

% Inverse Joukowsky
z1 = zeta1/2 + sqrt((zeta1/2).^2-1);
z2 = zeta2/2 - sqrt((zeta2/2).^2-1);
z3 = zeta2/2 + sqrt((zeta2/2).^2-1);
z4 = zeta1/2 - sqrt((zeta1/2).^2-1);


zeta = [zeta1 zeta2];
z = [z1 z2 fliplr(z3) fliplr(z4)];


figure()
%plot(real(zeta), imag(zeta))
hold on;
plot(real(z), imag(z))
axis equal






