close all; clc; clear;
format shortG

Isp = 250;
g = 32.2;
T = 300;
tp = 5.5;

mdot = T/Isp;
mp = T/Isp * tp;
mf = 30;
mi = mf + mp;

v = 0;
h = 0;
m = mi;
t = 0;
tstep = 0.01;
Cd = 0.2;
rho = .002377;
A = pi*4;

while (t < tp) || (v > 0)

a = - g - (1/2*rho*Cd*A*v^2)/m;
if m > (mf + mdot*tstep)
    a = (T/m)*g - g - (1/2*rho*Cd*A*v^2)/m;
    m = m - mdot*tstep;
end

v = v + a*tstep;
h = h + v*tstep;

t = t + tstep;

end

t = 0;
gloss = g*tp;
vn = 0;

while (t<tp)
    a = (T/m)*g - g + 1/2*rho*Cd*A*(v^2)/m;
    vn = vn + a*tstep;
    m = m - mdot*tstep;
    t = t+tstep;
end


disp(h)