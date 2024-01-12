clear; clc; close all;
format shortG
tic;

%% EDITABLE

%diam = Diameter of Parachute
%coefDrag = Coefficient of Drag. Usually from 0.5 to 2.2, provided by parachutes manufacturers usually, if unknown, assume a lower number rather than higher.
%dryMass = Weight of the Rocket (lbs) or (kg)
%angInc = Inclination in degrees (usually 0)
%lateralVelocity = Velocity of rocket at apogee in combined XZ plane
%delay on parachute deployment and inflation
%I really hope everything else is self explanatory

coefDragDrogue = 2.2;
coefDragMain = 2.2;
diamDrogue = 4; 
diamMain = 14;
timeFinal = 25;
timeStepsPerSecond = 100;
dryMass = 100;
angInc = 0;
lateralVelocity = 0;
delay = 1;
inflationDelay = 4;

%% Mode of Selection

tag = "both";
SI = false;
clc;

%% DO NOT EDIT

T = 0;                          %Start at T = 0, for simplicity
dT = 1/timeStepsPerSecond;      %Step of DT, as now defined by the caller of the function
g = 32.2;                       %ft/s2
gc = 1;                         % used to balance out SI units if chosen
densAir = .002377*g;            %slug/ft3 * ft/s2 == lb/ft3

if SI == true
    g = 9.81;
    gc = 9.81;
    densAir = 1.225;
end

coefDgArMain = coefDragMain*pi*(diamMain/2)^2;                              %CdS for Main
coefDgArDrogue = coefDragDrogue*pi*(diamDrogue/2)^2;                        %CdS for Drogue
drogueTerminalVelocity = sqrt((2*dryMass*g)/(densAir*coefDgArDrogue));      %Terminal Velocity using known equation, only used for internal debugging for if drogue event has been properly modeled
mainTerminalVelocity = sqrt((2*dryMass*g)/(densAir*coefDgArMain));

ratioCoefMD = coefDgArDrogue/coefDgArMain;                              

CdS = coefDgArMain * ratioCoefMD;                                       
X = ratioCoefMD;
Fs = .5*densAir*(drogueTerminalVelocity^2)*coefDgArMain;

XM = (2*dryMass)/(densAir*g*coefDgArMain*drogueTerminalVelocity*timeFinal); %unitless factor
J = 6;
n = 0;
drogueData = 0;
mainData = 0;

%% Drogue Parachute Deployment

if (contains(tag,"drogue",'IgnoreCase',true) == 1) || (contains(tag,"both",'IgnoreCase',true) == 1)

Fs = .5*densAir*(drogueTerminalVelocity^2)*coefDgArDrogue;
sPerDrogue = 0;
Vel = 0;
areaDragDrogueRatio = 0;
F = 0;

for i = 1:length(T:dT:timeFinal)
    T(i) = i*dT;

    if (T(i) > delay && T(i) < (delay + inflationDelay)) 
        Tx = T(i) - delay;
        sPerDrogue = (1-areaDragDrogueRatio)*((Tx/inflationDelay)^J) + areaDragDrogueRatio; %How open the parachute is at certain moments
    end

    CdS = coefDgArDrogue * sPerDrogue;
    dV(i,1) = - densAir * Vel^2 * CdS * sind(angInc) / (2*dryMass);
    dV(i,2) = (g - densAir * Vel^2 * CdS * cosd(angInc) / (2*dryMass));

    if i == 1
        V(i,1) = lateralVelocity;
        V(i,2) = 0;
        R(i,1:2) = 0;
        DS(i,1) = 0;
        continue
    end

    V(i,:) = V(i-1,:) + dV(i,:)*dT;
    R(i,:) = R(i-1,:) + V(i,:)*dT;
    DS(i) = DS(i-1) + (norm(R(i,:) - R(i-1,:)));
    Vel = norm(V(i,:));

    X = ((Vel/drogueTerminalVelocity)^2)*(CdS/coefDgArDrogue);  %Shock Factor
    F(i) = X*Fs*gc/g;
end
    
    drogueData = [T', V, R, F', F'/(gc*dryMass)];

    clearvars Tx V R Vel F X Ds dV CdS Tx T
    n = 6;
end

%% Main Function

if (contains(tag,"main",'IgnoreCase',true) == 1) || (contains(tag,"both",'IgnoreCase',true) == 1)

T = 0;
Vel = drogueTerminalVelocity;
sPerMain = 0;
Fs = .5*densAir*(drogueTerminalVelocity^2)*coefDgArMain;
F = 0;

for i = 1:length(T:dT:timeFinal)  
    T(i) = i*dT;
    if T(i) < (inflationDelay)
        Tx = T(i);
        sPerMain = (1-ratioCoefMD)*((Tx/inflationDelay)^J) + ratioCoefMD;            %How open the parachute is at certain moments
    end

    CdS = coefDgArMain * sPerMain;
    dV(i,1) = - densAir * Vel^2 * (CdS + coefDgArDrogue) * sind(angInc) / (2*dryMass);
    dV(i,2) = (g - densAir * Vel^2 * (CdS + coefDgArDrogue) * cosd(angInc) / (2*dryMass));

    if i == 1
        V(i,1) = 0;
        V(i,2) = drogueTerminalVelocity;
        R(i,1:2) = 0;
        DS(i,1) = 0;
        continue
    end

    V(i,:) = V(i-1,:) + dV(i,:)*dT;
    R(i,:) = R(i-1,:) + V(i,:)*dT;
    DS(i) = DS(i-1) + (norm(R(i,:) - R(i-1,:)));
    Vel = norm(V(i,:));

    X = ((Vel/drogueTerminalVelocity)^2)*(CdS/(coefDgArMain + coefDgArDrogue));  %Shock Factor
    F(i) = X*Fs*gc/g;
end

    mainData = [T', V, R, F', F'/(gc*dryMass)];

end

ParachuteShockPlot(drogueData,mainData,tag,n,delay,inflationDelay, SI);
toc;