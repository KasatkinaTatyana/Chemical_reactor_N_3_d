% �������� ������� � �������� �����������
close all
clear all
clc
global x_0 x_End
global y_0 y_End
global t_0
global R E1 E2 A10 A20
global k b alpha betta
global d
global a0 a1 b0 b1 g1 g2
%--------------------------------------------------------------------------
%-------------------------��������� ��������-------------------------------
R=8.3;
E1=2.09e4;
E2=4.18e4;
A10=1.1;
A20=172.2; 
%----------------------������������ ����������-----------------------------
a0=4.3145;
a1=-0.1099;
b0=1.4962;
b1=0.0515;
g1=41.8;
g2=83.6;
%--------------------------------------------------------------------------
global N    %����� ���������
N=50000;
t_0=0;
%--------------------------------------------------------------------------


x_0=[0.993 0.007 298];
% x_End=[0.35 0.55 333];

x_End=[0.693 0.3 328];
condition = 1;
while (condition)
%--------------------------------------------------------------------------
%��������� �������� ������������ ����������
y_0=YTrans(x_0);
y_End=YTrans(x_End);
%--------------------------------------------------------------------------
k=(y_End(2)-y_0(2))/(y_End(1)-y_0(1));
b=(y_0(2)*y_End(1)-y_End(2)*y_0(1))/(y_End(1)-y_0(1));
alpha = (y_0(3)-k*y_0(2))/y_0(2)/(y_0(1)-y_End(1));
betta = (y_0(2)*y_End(3)+y_0(3)*y_End(2)-2*k*y_0(2)*y_End(2))/y_0(2)/y_End(2)/(y_End(1)-y_0(1))^2;
%%
conf=DConf();
if (conf==0)
    break
else
    x_End
    x_End(1)=x_End(1)-0.045;
    x_End(2)=x_End(2)+0.01;
end

end