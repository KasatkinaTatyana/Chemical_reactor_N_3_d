function F=gKan(y)
% y=[y y' y'']
global R E1 E2

T=InvStateDiff(y);

F= y(2)^2 * k2(T) * (E2- E1) / R/T^2 /( k1(T) * y(1)^2 );
