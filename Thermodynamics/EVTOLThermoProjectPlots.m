EVTOLThermoProject;
close all;


figure(1)
hold on; grid on;
title('ESC Temperatue over time')
txt = {['Initial Temperature: ' num2str(T_ESC_start) , '°C'], ['Ambient Temperature: ' num2str(T_inlet), '°C']}
subtitle(txt)
xlabel('Time (s)');
yyaxis left; ylabel('Temperature (°C)');
plot(EphQ(:,2),EphQ(:,3))
yyaxis right; ylabel('Change in Temp (°C/s)')
plot(EphQ(:,2),EphQ(:,4)/tstep)
xlim([min(EphQ(:,2)) max(EphQ(:,2))])

% figure(2)
% hold on; grid on;
% plot(EphQ(:,2),EphQ(:,3))
% title('ESC Temperatue over time')
% txt = {['Initial Temperature: ' num2str(T_ESC_start) , '°C'], ['Ambient Temperature: ' num2str(T_inlet), '°C']}
% subtitle(txt)
% xlabel('Time (s)'); ylabel('Temperature (°C)');
% xlim([min(EphQ(:,2)) max(EphQ(:,2))])

figure(3)
hold on; grid on;
v = 0:.1:30;
k = 18.2 + 8.19*log(v);
vd = [2.5, 5, 7.25, 15, 20];
kd = [25, 32, 35, 40, 42.5];
plot(vd,kd,'o',v,k)
title('Relationship between Air Velocity and Thermal Heat Coefficient')
subtitle('Trend: k = 18.2 + 8.19ln(v)')
xlabel('Velocity (m/s)');ylabel('Thermal Heat Coefficient (W/°Cm^2)')
legend('Data Points','Trendline','location','northwest')
ylim([0,50])
