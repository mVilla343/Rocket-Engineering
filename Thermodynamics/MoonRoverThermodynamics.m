clc; close all; clear;

% Thermodynamic Analysis for a Spacecraft Design class, using notes from SMAD, NEW SMAD, and AIAA BROWN. Presented only with sanitizations of comments.

%% Stuff Engineers would reasonably be in control of, Edit:

SurfArea = 10.7416;
DissipationMax = 184;
DissipationMin = 92;

SolAbsorp = 0.163;
IRemis = 0.8;

TempLimit = [-20 50];
bodyTag = "Moon";

case1tag = 0;
case2tag = 1;

SetAlt = 0.1;

%% View Settings

n = 1;
FaceAlpha = 1;

%% DNE, immutable vars ahead

switch bodyTag
    case "Moon"
        bodyRadius = 1740;
        rmax = 5e3;
        Flux = [216 258];
        H = [0 rmax-bodyRadius];
        Albedo = 0.3;
    case "Earth"
        bodyRadius = 6380;
        rmax = 1e5;
        Flux = [216 258];
        H = [0 rmax-bodyRadius];
        Albedo = 0.07;
    otherwise
        disp('No Body');
end

SolFluxDir = 1378;
Alt = min(H):10^n:max(H);
SCdiam =  2*(SurfArea / (4*pi))^0.5;

ViewFactor = zeros(length(Alt));
WCH = zeros(length(Alt), 1);
WCC = zeros(length(Alt), 1);

% bunch of nerd stuff from brown book, Section 7.5, starts page 395

if case1tag == true

    for i = 1:length(Alt)
        ViewFactor(i) = 0.5 * (1 - ((Alt(i)^2 + 2 * Alt(i) * bodyRadius)^0.5 / (Alt(i) + bodyRadius)));
        AlbedoCorrection = 0.657 + 0.54 * (bodyRadius/(bodyRadius + Alt(i))) - 0.196 * (bodyRadius/(bodyRadius + Alt(i)))^2;
        WCH(i) = ((SolFluxDir * SolAbsorp/4 + max(Flux) * IRemis * ViewFactor(i) + SolFluxDir * Albedo * SolAbsorp * AlbedoCorrection(i) * ViewFactor(i) + DissipationMax/(pi * SCdiam^2)) / (0.0000000567 * IRemis))^0.25 - 273.15;
        WCC(i) = ((min(Flux) * IRemis * ViewFactor(i) + DissipationMin / (pi * SCdiam^2)) / (0.0000000567 * IRemis))^0.25 - 273.15;
    end
    
    figure(1)
    grid on; hold on;
    plot(Alt,WCH)
    title('Altitude of S/C and hot conditions')
    xlabel('Altitude (km)'); ylabel('Temperature (C)');
    
    figure(2)
    grid on; hold on;
    plot(Alt,WCC)
    title('Altitude of S/C and cold conditions')
    xlabel('Altitude (km)'); ylabel('Temperature (C)');

end

if case2tag == true

    SolAbs2 = linspace(0.01,1);
    IRe2 = SolAbs2;

    ViewFactor2 = 0.5 * (1 - ((SetAlt^2 + 2 * SetAlt * bodyRadius)^0.5 / (SetAlt + bodyRadius)));
    AlbedoCorrection2 = 0.657 + 0.54 * (bodyRadius/(bodyRadius + SetAlt)) - 0.196 * (bodyRadius/(bodyRadius + SetAlt))^2;

    WCH2 = zeros(length(SolAbs2),length(IRe2));
    WCC2 = zeros(length(SolAbs2),length(IRe2));

    for i = 1:length(SolAbs2)
        for j = 1:length(IRe2)
        WCH2(i,j) = ((SolFluxDir * SolAbs2(i)/4 + max(Flux) * IRe2(j) * ViewFactor2 + SolFluxDir * Albedo * SolAbs2(i) * AlbedoCorrection2 * ViewFactor2 + DissipationMax/(pi * SCdiam^2)) / (0.0000000567 * IRe2(j)))^0.25 - 273.15;
        WCC2(i,j) = ((min(Flux) * IRe2(j) * ViewFactor2 + DissipationMin / (pi * SCdiam^2)) / (0.0000000567 * IRe2(j)))^0.25 - 273.15;
        end
    end
    
    figure(3)
    hold on; grid on;
    CH = surf(IRe2, SolAbs2, WCH2, 'FaceAlpha', FaceAlpha);
    CH.EdgeColor = 'none';
    colorbar;
    view(2); axis tight;
    title('Worst Case Hot vs IR Emissions and Solar Absorption')
    xlabel('IR Emissitivity'); ylabel('Solar Absorption'); zlabel('Temperature (C)');

    figure(4)
    hold on; grid on;
    CC = surf(IRe2, SolAbs2, WCC2, 'FaceAlpha', FaceAlpha);
    CC.EdgeColor = 'none';
    colorbar;
    view(2); axis tight;
    title('Worst Case Cold vs IR Emissions and Solar Absorption')
    xlabel('IR Emissitivity'); ylabel('Solar Absorption'); zlabel('Temperature (C)');
end

