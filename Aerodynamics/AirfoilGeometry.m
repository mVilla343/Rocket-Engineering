clear
clc
close all

% Edit 1.11.24: This code shows the source panel method for Van der Walls 
% airfoils, completed as an assignment for an aerodynamics course.
% To the best of knowledge, this course no longer uses this assignment,
% so I post this knowing that nobody can use it as such, and as a precursor
% to NACAFoilPerformance.m, which is an alteration and major extension of
% this script.

%% EDIT BELOW

En = .0222; %.0222, .04725, .0721 
aln = 1; %3,9,15
k = 1.889;
CL = 1;
S = 360;
theta = 0;
an = 2*CL*(1.00+En).^(k-1.00)/(2.00.^k); %.574; %.551, .563, .574
n = 50; %10, 20, 30, 50, 100, 250
fpcount = 2*n;

% NO EDITING
%% PART I
for R = 1:length(aln)
for L = 1:length(En)
%Variable Call of respective Epsilon and a
E = En(L);    
a = an(L);

%Computation of the airfoil geometry
for i = 1:S

%Module I1
%from Handout p1
r1 = sqrt((a*cos(theta)-a)^2 + a^2*sin(theta)^2);
r2 = sqrt((a*cos(theta)-E*a)^2 + a^2*sin(theta)^2);
O1 = atan(a*sin(theta)/(a*cos(theta)-a))+pi;
O2 = atan((sin(theta))/(cos(theta)-E));
dr = cos(theta)-E;

n1 = 0;
if (dr < 0 && O2 < pi()/2)
    n1 = 1;
elseif (dr > 0 && theta > pi()/2)
    n1 = 2;
end

O2 = O2+n1*pi();

%from Handout p2
x(i,L) = (r1^k)/(r2^(k-1))*(cos(k*O1)*cos((k-1)*O2)+sin(k*O1)*sin((k-1)*O2)) + 1;
z(i,L) = (r1^k)/(r2^(k-1))*(sin(k*O1)*cos((k-1)*O2)-cos(k*O1)*sin((k-1)*O2));

%MODULE I2
al = aln(R)*pi()/180;
A = cos(O1*(k-1))*cos(k*O2)+sin(O1*(k-1))*sin(k*O2);
B = sin(O1*(k-1))*cos(k*O2)-cos(O1*(k-1))*sin(k*O2);
D0 = a*(1-k+k*E);
D1 = A*(a*cos(theta)-D0) - B*a*sin(theta);
D2 = B*(a*cos(theta)-D0) + A*a*sin(theta);

%From Handout
u = 2 * (r2^k)/(r1^(k-1)) * (sin(al)-sin(al-theta))/(D1^2+D2^2) * (D1*sin(theta)+D2*cos(theta));
w = -2 * (r2^k)/(r1^(k-1)) * (sin(al)-sin(al-theta))/(D1^2+D2^2) * (D1*cos(theta) - D2*sin(theta));
cp(i,L) = 1 - (u^2 + w^2);
cp(isnan(cp)) = 1;

theta = 2*pi()*i/S;
end
end
end

%% START OF PART II

%Panel Method for X and Z
x(isnan(x)) = 1;
linX = CL;
dlinX = 2*CL/n;
dx(1) = linX;
for i = 1:n
dx(i+1) = linX - dlinX*i; 
end
dz = interp1(x(1:length(x)/2),z(1:length(x)/2),dx);
dz(isnan(dz))=0; %Variable Correction

l = length(dx);
dx((l+1):(2*l-1)) = flip(dx(1:l-1));
dz((l+1):(2*l-1)) = -flip(dz(1:l-1));  %Bottom Side
clearvars l

%Zeta and Eta calculations, along with DeltaS and Orientation Theta
for i = 1:2*n
    zeta(i) = (dx(i+1)+dx(i))/2; 
    eta(i) = (dz(i+1)+dz(i))/2;
    X = dx(i+1)-dx(i);
    Z = dz(i+1)-dz(i);
    Dsj(i) = sqrt(X^2+Z^2);
    thetai(i) = atan2d(Z,X);
    if thetai(i) < 0
        thetai(i) = 360 + thetai(i);
    end
    clearvars X Z
end 

%Determining of Phi, R, B, C, and CBAR, along with correspond A
for j = 1:2*n
    for i = 1:2*n
        X = zeta(i)-zeta(j);
        Z = eta(i)-eta(j);
        phii(i,j) = atan2d(Z,X);
        if phii(i,j) < 0
            phii(i,j) = 360 + phii(i,j);
        end
        dD = phii + 90;
        R(i,j) = sqrt(X^2+Z^2);
        B(i) = -sind(thetai(i));
        C(i,j) = Dsj(i)*sind(thetai(i)-phii(i,j))/(2*pi()*R(i,j));
        C_(i,j) =  Dsj(j)*cosd(thetai(i)-phii(i,j))/(2*pi()*R(i,j));
        A(i,j) = C(i,j);
        if i == j
            A(i,j) = 1/2;
        end
        C_(isinf(C_))=0;
        clearvars X Z;
    end
end

q = A\B';
%Finding of Numerical CoefficientP
for i = 1:2*n
    for j = 1:2*n
    qC_(i,j)=(q(j)*C_(i,j));
    end
    vti(i) = cosd(thetai(i)-aln)+sum(qC_(i,:));
    Cp(i) = 1 - (vti(i))^2;
    zeta(i);
end

%% PLOTTING OF FOUND GEOMETRIES

%Vanfoil
figure(1)
hold on
grid on
xlabel('X', FontSize=14)
ylabel('Z', FontSize=14)
plot(x(:,L),z(:,L))
title('Airfoil Geometry')
subtitle(['AoA = ' num2str(al), ' E = ' num2str(E)])
set(gcf,'units','normalized','position',[0 0 1 1])
axis equal
hold off

%Theoretical
figure(2)
hold on
grid on
xlabel('X', FontSize=14)
ylabel('Cp', FontSize=14)
plot(x,cp)
title('Coefficient of Pressure for Airfoil')
subtitle(['AoA = ' num2str(aln), ', E = ' num2str(E)])
set(gcf,'units','normalized','position',[0 0 1 1])
hold off

%Panel Method
figure(3)
grid on, hold on
title(['Actual Airfoil vs Source Panel Airfoil, ' num2str(fpcount), ' Panels'])
plot(x,z)
plot(dx,dz, ':.', "Color", [0 0 0])
plot(zeta,eta,'o', "Color", [.8 .2 .3])
legend('Van de Vooren', 'Source Panel', 'Panel Midpoints')
set(gcf,'units','normalized','position',[0 0 1 1])
axis equal
hold off
    
%Numerical CoefficientP
figure(4)
hold on
grid on
xlabel('X', FontSize=14)
ylabel('Cp', FontSize=14)
plot(x,cp,zeta,Cp,':square')
title(['Numerical Center of Pressure with ', num2str(fpcount), ' Panels'])
subtitle(['AoA = ' num2str(aln), ', E = ' num2str(E)])
set(gcf,'units','normalized','position',[0 0 1 1])
legend('Theoretical', 'Numerical')
hold off


