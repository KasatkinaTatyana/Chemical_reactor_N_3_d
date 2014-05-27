function F=InvStateDiff(y)
%T=InvStateDiff(y), ��� y=[y y' y''];
%������ ������� ������ ����������� �� ������������ ���������� � ��������
global R E1 E2 A10 A20

T=(E1-E2)/R/(log(A10*(y(3)-y(2))*y(1)^2/A20/y(2)^2));

% if ((y(3)-y(2))<=0)
%     disp(y(3)-y(2));
% end

F=T;