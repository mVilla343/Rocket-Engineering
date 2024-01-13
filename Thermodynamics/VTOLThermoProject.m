clear; clc;

% This script was a collaborative effort with Stanley Ossyra (BSAE CPP
% 2024, future MAE, 6'1", german/argentinian, real nice guy) for his E-VTOL 
% project which uses controllers for this. While he set the framework, he
% fell ill and passed it onto me to be completed on time. This script is
% the brunt of our analysis, visualized with VTOLThermoProjectPlots.m

T_inlet = 25;           % Inlet air temperature (°C)
T_ESC_desired = 95;     % Desired ESC operating temperature (°C)
T_ESC_start = 20;       %Starting Temp of ESC, which can ambient or something in between
A_ESC = 0.01;           % ESC surface area (m^2)
emissivity = 0.85;      % ESC surface emissivity
d_hydraulic = 0.05;     % Hydraulic diameter of duct (m)
L_duct = 1;             % Length of duct (m)
nu = 15.6e-6;           % Kinematic viscosity of air (m^2/s)
k = 6.763e-2;           % Thermal conductivity of air (W/mK)
v_air = 20;             % Initial guess for velocity of air (m/s)
C = 290;                % Specific heat of ESC (J/kgC or Ws/kgC)
m = .11;                % Mass of ESC (kg)
v_air_max = 40;         % Max Speed of Air through duct possible
efficiency = 97.5;      %Efficiency of the ESC

% Calculate Reynolds number and friction factor for initial guess
Re = (v_air * d_hydraulic) / nu;
f = 0.0791 / Re^(1/4);

% Iterate until T_ESC matches the desired temperature or maximum number of iterations is reached
iter = 1;
tstep = .1;
max_iter = 1000/tstep;
t = 0;
T_ESC = T_ESC_start;
Q = C*m*T_ESC;
dQ = 1;
while abs(dQ) > 0.01 || (T_ESC_desired-T_ESC) < 0
    % Calculate heat transfer and pressure drop
    h = 18.2+8.19*log(v_air);
    qConvection = A_ESC * (T_ESC - T_inlet) * (h); %change in heat from convection (Watts or J/S)
    qRadiation = A_ESC * emissivity * 5.67e-8 * (T_ESC^4 - T_inlet^4); %Radiation heat away (Watts or J/s)
    dq_total = (1-efficiency/100)*1040 - (qConvection + qRadiation); %1040 Watts power the ESC, however most of that is not heating. Assuming that all function are electrical, most of tha
    pressureDrop = 4 * f * L_duct * (v_air^2 / (2 * d_hydraulic));
    q_flow = v_air * A_ESC * 1000; % Calculate required air flow rate
    dQ = tstep*dq_total;    %Change in Heat over Tstep
    dT = dQ/(m*C);          %Change in Temp
    Q = Q + dQ;             %New Heat
    T_ESC = Q/(m*C);        %Temperature

    % Check if T_ESC matches desired temperature
    if iter>max_iter % Tolerance of 0.1°C
        break
    end

    % Adjust velocity and Reynolds number for next iteration
    if v_air > v_air_max
        v_air = v_air_max;
    else
        if T_ESC > T_ESC_desired
            v_air = v_air + 0.1;
        else
            v_air = v_air + 0;
        end
    end
    Re = (v_air * d_hydraulic) / nu;
    f = 0.0791 / Re^(1/4);

    t = t+tstep;
    EphQ(iter+1,:) = [iter;t;T_ESC;dT;dq_total;dQ;h;v_air];
    % Increment iteration counter
    iter = iter + 1;
end
EphQ(1,:) = [0,0,T_ESC_start,0,0,0,0,20];

% Output results
EphQT = array2table(EphQ);
% Default heading for the columns will be A1, A2 and so on. 
% You can assign the specific headings to your table in the following manner
EphQT.Properties.VariableNames(1:8) = {'Iteration','Time','Temperature','dTemp','Heat flow total', 'Heat Added','h', 'Velocity'};

if abs(T_ESC_desired-T_ESC) > 1 && (T_ESC-T_ESC_desired) > 0  
    fprintf('Error: Desired temperature cannot be reached with given parameters.\n')
    fprintf('Final Temperature after %0.0f seconds: %0.2f°C\n',t,T_ESC)
else
    fprintf('Required air flow rate: %0.2f L/min\n', q_flow)
    fprintf('Pressure drop in duct: %0.2f Pa\n', pressureDrop)
    fprintf('Reynolds number: %0.2f\n', Re)
    fprintf('ESC operating temperature: %0.2f°C\n', T_ESC)
end
