clear; clc; close all;

Pref = 4.875;
pik = 0.20/100;
a = 0.561;
n = 0.35;
Ro = 1760/100^3;
c = 1511;

t = 0;
dt = 0.1;
W = 0;
i = 1;

At = 0.9676*100^2;
L = 40*100;
ri = 0.6*100;
ws = 0.2*100;
rg = 3.7/2*100;
as = (pi/2*rg - ws)/(2+pi/2); 


while W < as

    t = t+dt;
%     pc(i) = ((a*Ro*c/At) * (2*pi*L/100*(ri+W)))^(1/(1-n));

    D = 2*pi*(ri+W) - 4*(ws+2*W);
    E = 8*(rg-as-ri-W);
    F = 4*pi*W;
    G = ws;
    pc(i) = ((a*Ro*c/At*L/100) * (D+E+F+G))^(1/(1-n));

    rb(i) = a*pc(i)^n;
    W = W + rb(i)*dt; 
    i = i+1;

end

t = linspace(0,dt*length(pc),length(pc));

figure(1)
plot(t,pc)

figure(2)
plot(t,rb)




