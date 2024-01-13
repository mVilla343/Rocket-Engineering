function [A] = PrandtlMeyerFunc(var1,var2,tag)
%Initialization of the PM Function, relating a turn of Mach to an expansion corner
g = 1.4;
Mach = 1:0.01:10;
Nu = sqrt((g+1)/(g-1)).*atand(sqrt( ((g-1)/(g+1)) * (Mach.^2-1) )) - atand(sqrt(Mach.^2-1));
% Theta = nu(Mach2) - nu(Mach1)
% What are we solving for?
if contains(tag,'O') == 1 || contains(tag,'theta') == 1
    M1 = var1; M2 = var2;
    Nu1 = interp1(Mach,Nu,M1);
    Nu2 = interp1(Mu,Nu,M2);
    O = Nu2 - Nu1; A = O;
%     fprintf('Value of Theta: %0.2f', O);
elseif contains(tag,'Mach1') == 1 || contains(tag,'M1') == 1
    M2 = var1; O = var2; 
    Nu2 = interp1(Mach,Nu,M2);
    Nu1 = Nu2 - O;
    M1 = intepr(Nu,Mach,Nu1); A = M1;
%     fprintf('Value of Mach 1: %0.2f', M1);
elseif contains(tag,'Mach2') == 1 || contains(tag,'M2') == 1
    M1 = var1; O = var2;
    Nu1 = interp1(Mach,Nu,M1);
    Nu2 = O + Nu1;
    M2 = interp1(Nu,Mach,Nu2); A = M2;
%     fprintf('Value of Mach 2: %0.2f', M2);
end
end
