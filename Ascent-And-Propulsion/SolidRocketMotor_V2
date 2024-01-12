clear; close all; clc;

WAb = readtable("WebAreaRelations.xls");
WAb = table2array(WAb);

Wref = WAb(:,1);
Abref = WAb(:,2);

a = 0.561;
n = 0.35;
Ro = 1760/100^3;
c = 1511;
At = 0.9676*100^2;
g = 1.4;
Erat = 10;

ti = 0;
dt = 0.01;
Wi = 0;
i = 1;
L = 40*100;
rg = 123;

while Wi < (0.8*rg)

    ti = ti+dt;
    Abi = interp1(Wref,Abref,Wi, "spline");
    pc(i) = ((a*Ro*c/At*L/100) * (Abi))^(1/(1-n));
    rb(i) = a*pc(i)^n;
    Wi = Wi + rb(i)*dt; 
    W(i) = Wi;
    t(i) = ti;
    Ab(i) = Abi;
    i = i+1;

end

for i = 1:length(pc)                                                                                                                  
    if pc(i) < 0.1013
        pc(i) = 0.1014;
    end 
    p2(i) = ExitPress(pc(i),Erat,g);
end

pt = pc / ((1 + 0.5*(g-1))^ (g/(g-1)));

mdot = pc*(10^6)*(At/100^2)/(c);
CF = sqrt(2*g^2 /(g-1) * (2/(g+1))^((g+1)/(g-1)) * (1 - p2./pc).^((g-1)/g)) + (p2 - 0.1013)./pc * Erat;
F = CF.*pc*10^6*(At/100^2);
Isp = mean(F)/(9.81*mean(mdot));
cstar(1:length(t)) = c;

figure(1)
hold on; grid on;
title('Chamber Pressure vs. Time')
plot(t,pc)
xlabel("Time (s)")
ylabel("Pressure (MPa)")

figure(2)
hold on; grid on;
title('Burning Rate vs Time')
plot(t,rb)
xlabel("Time (s)")
ylabel("r_b (cm/s)")

figure(3)
hold on; grid on;
title('Characteristic Velocity vs. Time')
plot(t, cstar)
xlabel("Time (s)")
ylabel("C^* (m/s)")

figure(4)
hold on; grid on;
title('Thrust vs. Time')
plot(t,F/10^6)
xlabel("Time (s)")
ylabel("F (MN)")

figure(5)
hold on; grid on;
title('Burn Area vs. Time')
plot(t,L*Ab/100^3)
xlabel("Time (s)")
ylabel("Ab (m^3)")

Pcmax = max(pc);
Fnmax = max(F);
% Z = [t; pc; rb; cstar; Ab]'
