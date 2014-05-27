function F=fKan(y)
% y=[y y' y'']
T=InvStateDiff(y);

CB1 = -y(2)-k2(T)/k1(T) * y(2)^2 / y(1)^2;

F = -CB1 + 2*k2(T)/k1(T) * (-y(2))/y(1)^2 * ( CB1 + y(2)^2/y(1));
