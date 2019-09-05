function [x,z] = naca4 (M,P,XX,n)
% NACA4 Geometrija NACA profila serije 4.
%	 Funkcija NACA4 definira koordinate profila NACA serije 4 za zadanu 
%	 maksimalnu zakrivljenost M (u stotinama tetive), apcisa maksimalne
%	 zakrivljenosti P (u desetinama tetive) te najvece debljine XX u stotinama 
%    tetive). Te su velicine dane u relativnom odnosu prema tetivi (c=1). 
%    Tako da prva tri ulazna podatka predstavljaju oznaku zeljenog NACA
%	 profila: NACA MPXX.
%    Pri tome funkcija generira ukupno 2*N+1 tocka (N tocaka na gornjaci i
%	 N tocaka na donjaci uz to da je tocka na izlaznom rubu prva, a ujedno i
%	 zadnja). Tocke su slozene, pocev od izlaznog ruba preko donjake pa nazad
%	 preko gornjake do izlaznog ruba. U pozivu oblika:
%
%		[X,Z] = NACA4(M,P,XX,N)
%
%	 u vektoru X spremljene su x koordinate točaka, a z koordinate točaka u
%	 vektoru Z.

% M.Vrdoljak, FSB, 2005, 2014.

c=1;

% najveca zakrivljenost srednje linije (relativno prema duljini tetive)
f=M/100;
% polozaj najvece zakrivljenosti srednje linije (rel.prema dulj.tetive)
xf=P/10;
% najveca debljina (rel.prema dulj.tetive)
t=XX/100;

% raspodjela xc koordinate duz tetive
beta = pi/2:pi/2/n:pi;
xc   = c*(1+cos(beta)); % cos-raspodjela: gusce oko LE

% promjena debljine = f(xc)
zt = 5*t*(0.2969.*sqrt(xc)-0.126.*xc-0.3537.*(xc).^2 +0.2843.*(xc).^3-0.1015.*(xc).^4);

% srednja linija = f(xc) i nagib srednje linije = f(xc)
for i=1:n+1
	if ((f==0) || (xf==0))
		% simetrican profil!
		% koordinate tocaka na gornjaci
		x_g(i) = xc(i);
		z_g(i) = zt(i);
		% koordinate tocaka na donjaci
		x_d(i) = xc(i);
		z_d(i) =-zt(i);
	else
		if (xc(i)/c<xf)
			zc(i)    = c*f/xf^2*(2*xf*xc(i)/c-(xc(i)/c)^2);
			theta(i) = atan(2*f/xf^2*(xf-xc(i)/c));
		else 
			zc(i)    = c*f/(1-xf)^2*(1-2*xf+2*xf*(xc(i)/c)-(xc(i)/c)^2);
			theta(i) = atan(2*f/(1-xf)^2*(xf-xc(i)/c));
		end
		% koordinate tocaka na gornjaci
		x_g(i) = xc(i) - zt(i)*c*sin(theta(i));
		z_g(i) = zc(i) + zt(i)*c*cos(theta(i));
		% koordinate tocaka na donjaci
		x_d(i) = xc(i) + zt(i)*c*sin(theta(i));
		z_d(i) = zc(i) - zt(i)*c*cos(theta(i));
	end

end

% preslagivanje vektora x i z koordinata
x = [x_d'; x_g(n:-1:1)'];
z = [z_d'; z_g(n:-1:1)'];
