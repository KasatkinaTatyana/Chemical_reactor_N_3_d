function F=s(x)
% ������ ������� s(x)=dt/d*tau
% x = [cA cB T]

F=x(2)/(k1(x(3))*x(1)^2);