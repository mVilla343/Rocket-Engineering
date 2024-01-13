function [Bc,Oc,P2P1,R2R1, T2T1] = ObliqueShockFunc(M1,B,O,g,tag)
Beta = 0:0.5:65;
Theta = atand((((M1^2)*sind(Beta).^2 - 1)./(M1^2*(g+cosd(2*Beta))+2))*2.*cotd(Beta));
% plot(Theta,Beta)
if contains(tag,'B') == 1
Bc = interp1(Theta,Beta,O);
Oc = O;
elseif contains(tag,'O') == 1
Oc = interp1(Beta,Theta,B);
Bc = B;
end
Mn1 = M1*sind(Bc);
P2P1 = 1 + 2*g/(g+1)* Mn1^2-1;
R2R1 = (g+1)*Mn1^2/((g-1)*Mn1^2 + 2);
T2T1 =  P2P1/R2R1;
CP = 2/(g*M1^2)*(P2P1-1);
end
% M1 = 3
% Beta = 0:0.5:90;
% Theta = atand(2*cotd(Beta)*((M1^2.*sind(Beta).^2-1)/(M1^2.*(g+cosd(2.*Beta))+2)));