close all; clear all;

m = 0.1;
c = 4;
p = 0.4;

x1 = 0:0.05:0.4*c;
x2 = 0.4*c:0.05:c; x2 = x2(2:end);
x = [x1 x2];
y1 = m/p^2*(2*p*x1/c-(x1/c).^2);
y2 = m/(1-p)^2*((1-2*p)+2*p*x2/c-(x2/c).^2);
x1 = x1-2; x2 = x2-2;

zeta = [complex(fliplr(x2), fliplr(y2)) complex(fliplr(x1)), fliplr(y1)))];

[~, ind] = find(angle(zeta)>0.5*pi, 1);

zeta1 = zeta(1:ind-1); zeta2 = zeta(ind:end);

% Inverse Joukowsky
z1 = zeta1/2 + sqrt((zeta1/2).^2-1);
z2 = zeta2/2 - sqrt((zeta2/2).^2-1);
z3 = zeta2/2 + sqrt((zeta2/2).^2-1);
z4 = zeta1/2 - sqrt((zeta1/2).^2-1);

z = [z1 z2 fliplr(z3(1:end-1)) fliplr(z4(2:end))];

r = abs(z);
t = angle(z);

psi = ln(r);



figure()
plot(real(zeta), imag(zeta))
hold on;
plot(real(z), imag(z))
axis equal



