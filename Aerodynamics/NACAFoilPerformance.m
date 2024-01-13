clear; clc; close all;

% Expansion of AirfoilGeometry.m, done as a personal project as I felt the
% previous did not capture a lot of the essence of what was being taught.
% This uses vortex panel method, which enforces the Kutta Condtion, and
% approximates lift coefficient, but not drag

%% Settings (N2: Panel #, alpha: AoA, XS: non-sim resolution of airfoil

N2 = 30;
alpha = 0;
XS = 0.005;

%% Do Not Edit

vinf = 1;
CC = rem(N2,2);
N2 = N2 - CC;

%% Input for 4SERIES

NACA4 = input('Enter NACA foil #: NACA ','s');
display(['NACA ', num2str(NACA4)]);

m = str2double(NACA4(1))/100; p = str2double(NACA4(2))/10;
t = str2double(NACA4(3))*10+str2double(NACA4(4));

%% SYMMETRIC 4SERIES
if m == 0 || p == 0
x = 0:XS:1;
yt = 5*t*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2 + 0.2843*x.^3 - 0.1015*x.^4)/100;
yU = yt; yL = -yt;
NACAX = [x, flip(x(1:end-1))]; NACAY = [yL,flip(yU(1:end-1))];
tag = "4SS";
end

%% Assymetric 4Series 
if m ~= 0
% Calculation of Camber Line 
xc1 = 0:XS:p; xc2 = (p+XS):XS:1; x = [xc1,xc2];
yc1 = (m/p^2) * (2*p*xc1 - xc1.^2);
yc2 = (m/(1-p)^2) * ((1-2*p)+2*p*xc2 - xc2.^2);
ycam =  [yc1,yc2];

% Calculation of Camber Line Slope
dycdx1 = (2*m/p^2) * (p - xc1);
dycdx2 = (2*m/(1-p)^2) * (p - xc2);
dycdx = [dycdx1, dycdx2];
CLS = atan(dycdx);

% Calculation of Symmetric Airfoil
yt = 5*t*(0.2969*sqrt(x)-0.1260*x-0.3516*x.^2 + 0.2843*x.^3 - 0.1015*x.^4)/100;

nCLS = length(CLS);
xU = zeros(1,nCLS); xL = zeros(1,nCLS); yU = zeros(1,nCLS); yL = zeros(1,nCLS);
for i = 1:length(CLS)
    xU(i) = x(i) - yt(i)*sin(CLS(i));
    xL(i) = x(i) + yt(i)*sin(CLS(i));
    yU(i) = ycam(i) + yt(i)*cos(CLS(i));
    yL(i) = ycam(i) - yt(i)*cos(CLS(i));
end

NACAX = [xL, flip(xU(1:end-1))]; NACAY = [yL,flip(yU(1:end-1))];
tag = "4SA";
end

%% PANEL METHOD | Primary Geometries

xstep = 2/N2;
xb = 1; xb(N2/2) = 0;
for i = 1:N2/2
    xb(i+1) = 1-xstep*i;
end
ybt = interp1(x,yU,xb);
ybl = interp1(x,yL,xb);
xb = [xb,flip(xb(1:end-1))]';
yb = [ybl,flip(ybt(1:end-1))]';


%% Secondary Geometries
xc=zeros(N2,1); yc=zeros(N2,1); S =zeros(N2,1); phi=zeros(N2,1);
for i = 1:N2
    xc(i,1) = (xb(i)+xb(i+1))/2;
    yc(i,1) = (yb(i)+yb(i+1))/2;
    DelX = xb(i+1)-xb(i);
    DelY = yb(i+1)-yb(i);
    S(i,1) = sqrt(DelX^2+DelY^2);
    phi(i) = atan2d(DelY,DelX);
    if phi(i) < 0
        phi(i) = 360 + phi(i);
    end
end
deln = phi + 90;
beta = deln - alpha;
beta(beta>360) = beta(beta>360) - 360;
xb = [xb;xb(1)]; yb = [yb;yb(1)];

%% Geometric Integrals I,J,K,L
I = zeros(N2,N2); J = zeros(N2,N2);
K = zeros(N2,N2); L = zeros(N2,N2);
for i= 1:N2
    for j = 1:N2
    if (j ~= i)
    A = -(xc(i)-xb(j))*cosd(phi(j)) - (yc(i)-yb(j))*sind(phi(j));
    B = (xc(i)-xb(j))^2 + (yc(i)-yb(j))^2;
    E = sqrt(B-A^2);

    C = sind(phi(i)-phi(j));
    CJ = -cosd(phi(i)-phi(j));
    CL = sind(phi(j)-phi(i));

    D = -(xc(i)-xb(j))*sind(phi(i)) + (yc(i)-yb(j))*cosd(phi(i));
    DJ = (xc(i)-xb(j))*cosd(phi(i)) + (yc(i)-yb(j))*sind(phi(i));
    DL = (xc(i)-xb(j))*sind(phi(i)) - (yc(i)-yb(j))*cosd(phi(i));

    if (~isreal(E))
        E = 0;
    end
    I(i,j) = (C/2) * log((S(j)^2 + 2*A*S(j) + B)/B) + ((D-A*C)/E)*(atan2((S(j)+A),E)-atan2(A,E));
    J(i,j) = (CJ/2) * log((S(j)^2 + 2*A*S(j) + B)/B) + ((DJ-A*CJ)/E)*(atan2((S(j)+A),E)-atan2(A,E));
    K(i,j) = (CJ/2) * log((S(j)^2 + 2*A*S(j) + B)/B) + ((DJ-A*CJ)/E)*(atan2((S(j)+A),E)-atan2(A,E));
    L(i,j) = (CL/2) * log((S(j)^2 + 2*A*S(j) + B)/B) + ((DL-A*CL)/E)*(atan2((S(j)+A),E)-atan2(A,E));
    end
    end
end
clearvars A B C CJ D DJ E
I(isnan(I))=0; I(isinf(I))=0; I(~isreal(I))=0;
J(isnan(J))=0; J(isinf(J))=0; J(~isreal(J))=0; 
K(isnan(K))=0; K(isinf(K))=0; K(~isreal(K))=0;
L(isnan(L))=0; L(isinf(L))=0; L(~isreal(L))=0; 

%% Source Strength Calculation
Ax = I + pi*eye(N2,N2);
Bx = -vinf*2*pi*cosd(beta);
lam = Ax\Bx;

%% Vortex Strength Calculation
Omx = -K;  Bxg = Bx;
Omx(N2,:) = 0; Omx(N2,1) = 1; Omx(N2,end) = 1; Bxg(N2) = 0; 
gamma = Omx\Bxg;

%% Velocity and Coefficent Pressure Calculations (Source and Vortex Panel Method)
vt = zeros(N2,1);
cp = zeros(N2,1);
for i = 1:N2
    Sum = 0;
    for j = 1:N2
        Sum = Sum + (lam(j)/(2*pi))*(J(i,j));
    end
    vt(i) = vinf*sind(beta(i)) + Sum;
    cp(i) = 1 - (vt(i)/vinf)^2;
end

vt_v = zeros(N2,1);
cp_v = zeros(N2,1);
for i = 1:N2
    Sum = 0;
    for j = 1:N2
        Sum = Sum - (gamma(j)/(2*pi))*L(i,j);
    end
    vt_v(i) = vinf*sind(beta(i)) + Sum + gamma(i)/2;
    cp_v(i) = 1 - (vt_v(i)/vinf)^2;
end

%% Coefficient Lift, Drag Calculations (Source and Vortex Panel Method)
cn = -cp.*S.*sind(beta);
cc = -cp.*S.*cosd(beta);

cl = sum(cn.*cosd(alpha)) - sum(cc.*sind(alpha));
cd = sum(cn.*sind(alpha)) + sum(cc.*cosd(alpha));

cn_v = -cp_v.*S.*sind(beta);
cc_v = -cp_v.*S.*cosd(beta);
cl_v = sum(cn_v.*cosd(alpha)) - sum(cc_v.*sind(alpha));
cd_v = sum(cn_v.*sind(alpha)) + sum(cc_v.*cosd(alpha));

%% Mxy Calculations (Source Panel Method)
ngX = 150; ngY = 150;
xlimit = [-0.5, 1.5]; ylimit = [-0.5, 0.5];
dss = 0.01; Pmax = ngX*ngY*10; SLP = 70;
Yst = linspace(ylimit(1), ylimit(2),floor((SLP/100)*ngY))';
Xg = linspace(xlimit(1), xlimit(2),ngX)'; Yg = linspace(ylimit(1), ylimit(2),ngY)';
[X,Y] = meshgrid(Xg,Yg);

Vx = zeros(ngX,ngY); Vy = zeros(ngX,ngY);
for m = 1:ngX
    for n = 1:ngY
        Px = X(m,n);
        Py = Y(m,n);
        Mx = zeros(N2,1); My = zeros(N2,1);
        for j = 1:N2
            A = -(Px-xb(j))*cosd(phi(j)) - (Py-yb(j))*sind(phi(j));
            B = (Px-xb(j))^2 + (Py-yb(j))^2;
            Cx = -cosd(phi(j)); Cy = -sind(phi(j));
            Dx = (Px-xb(j)); Dy = (Py-yb(j));
            E = sqrt(B-A^2);
            if (~isreal(E))
                E = 0;
            end

            Mx(j) = (Cx/2) * log((S(j)^2 + 2*A*S(j) + B)/B) + ((Dx-A*Cx)/E)*(atan2((S(j)+A),E)-atan2(A,E));
            My(j) = (Cy/2) * log((S(j)^2 + 2*A*S(j) + B)/B) + ((Dy-A*Cy)/E)*(atan2((S(j)+A),E)-atan2(A,E));
        end
        Mx(isnan(Mx))=0; Mx(isinf(Mx))=0; Mx(~isreal(Mx))=0; 
        My(isnan(My))=0; My(isinf(My))=0; My(~isreal(My))=0;  
        [in,on] = inpolygon(Px,Py,xb,yb);
        if (in == 1 || on == 1)
            Vx(m,n) = 0;
            Vy(m,n) = 0;
        else
            Vx(m,n) = vinf*cosd(alpha) + sum(lam.*Mx./(2*pi));
            Vy(m,n) = vinf*sind(alpha) + sum(lam.*My./(2*pi));
        end
    end
    Vxy = sqrt(Vx.^2 + Vy.^2);
    Cpxy = 1 - (Vxy./vinf).^2;
end

%% Nxy Calculations (Vortex Panel Method)
Vvx = zeros(ngX,ngY); Vvy = zeros(ngX,ngY);
for m = 1:ngX
    for n = 1:ngY
        Px = X(m,n);
        Py = Y(m,n);
        Nx = zeros(N2,1); Ny = zeros(N2,1);
        for j = 1:N2
            A = -(Px-xb(j))*cosd(phi(j)) - (Py-yb(j))*sind(phi(j));
            B = (Px-xb(j))^2 + (Py-yb(j))^2;
            E = sqrt(B-A^2);
            if (~isreal(E))
                E = 0;
            end
            Cx = sind(phi(j)); Cy = -cosd(phi(j));
            Dx = -(Py-yb(j)); Dy = (Px-xb(j));

            Nx(j) = (Cx/2) * log((S(j)^2 + 2*A*S(j) + B)/B) + ((Dx-A*Cx)/E)*(atan2((S(j)+A),E)-atan2(A,E));
            Ny(j) = (Cy/2) * log((S(j)^2 + 2*A*S(j) + B)/B) + ((Dy-A*Cy)/E)*(atan2((S(j)+A),E)-atan2(A,E));
        end
        Nx(isnan(Nx))=0; Nx(isinf(Nx))=0; Nx(~isreal(Nx))=0; 
        Ny(isnan(Ny))=0; Ny(isinf(Ny))=0; Ny(~isreal(Ny))=0;  
        [in,on] = inpolygon(Px,Py,xb,yb);
        if (in == 1 || on == 1)
            Vvx(m,n) = 0;
            Vvy(m,n) = 0;
        else
            Vvx(m,n) = vinf*cosd(alpha) + sum(-gamma.*Nx./(2*pi));
            Vvy(m,n) = vinf*sind(alpha) + sum(-gamma.*Ny./(2*pi));
        end
    end
    Vxy_v = sqrt(Vvx.^2 + Vvy.^2);
    Cpvxy = 1 - (Vxy_v./vinf).^2;
end

%% PLOTTING

figure(1)
hold on; grid on; axis equal;
title("NACA " + NACA4 + " Airfoil")
plot(NACAX,NACAY,'Color',[0 0 0])
plot(xb,yb,'Color',[0 0 1])
plot(xc,yc,'o','Color',[1 0 0])
legend("Airfoil", "Airfoil Approximation", "Airfoil Approx. Midpoints");
if exist("ycam",'var')
    plot(x,ycam,'Color',[0 1 0])
    legend("Airfoil", "Airfoil Approximation", "Airfoil Approx. Midpoints", "Camber Line");
end

figure(2)
hold on; grid on; 
title("NACA " + NACA4 + " Airfoil Flow under Source Panel Method")
subtitle("Angle of Attack: " + alpha)
% quiver(X,Y,Vx,Vy,'r');
for i = 1:length(Yst)
    sline = streamline(X,Y,Vx,Vy,xlimit(1),Yst(i),[dss,Pmax]);
    set(sline,'LineWidth',1)
end
plot(xb,yb,'k');
xlim([xlimit(1),xlimit(2)]); ylim([ylimit(1),ylimit(2)]); zoom reset;

figure(3)
hold on; grid on; 
title("NACA " + NACA4 + " Airfoil Flow under Source Panel Method")
subtitle("Angle of Attack: " + alpha)
%quiver(X,Y,Vvx,Vvy,'r');
for i = 1:length(Yst)
    sline = streamline(X,Y,Vvx,Vvy,xlimit(1),Yst(i),[dss,Pmax]);
    set(sline,'LineWidth',1)
end
plot(xb,yb,'k');
xlim([xlimit(1),xlimit(2)]); ylim([ylimit(1),ylimit(2)]); zoom reset;

disp(cl_v)