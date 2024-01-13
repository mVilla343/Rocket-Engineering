function [ExitPressure] = ExitPress(p1, Erat, g)

pyp1 = 0.001:0.001:0.9;
E = 1./(((g+1)/2 ^ (1/(g-1))) .* (pyp1).^(1/g) .* sqrt(((g+1)/(g-1)) .* ( 1- (pyp1).^((g-1)/g))));

PRatio = interp1(E, pyp1, Erat);
ExitPressure = p1 * PRatio;
end
