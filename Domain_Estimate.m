% ������ ������� ������������ ������� 
% ���� ������ ������ 12.10.2014
% �������� ������� � �������� �����������
close all
clear all
clc
global x_0 x_End
global y_0 y_End
global t_0
global R E1 E2 A10 A20
global k b alpha betta
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
t_0=0;
%--------------------------------------------------------------------------


x_0=[0.993 0.007 298];
% x_End=[0.35 0.55 333];

x_End=[0.4 0.28 333];
%--------------------------------------------------------------------------
%��������� �������� ������������ ����������
y_0=YTrans(x_0);
y_End=YTrans(x_End);
%--------------------------------------------------------------------------
   
N=10000;           %����� ���������
% dT = (0.01 - y_0(1)) / N;
dT = -0.0001;
T=zeros(1,N+1);
Y=zeros(1,N+1);
T(1) = y_0(1);
Y(1) =  y_0(2);
t = y_0(1);
y_sing = 1e-1; % �������� y ��� ������� ������ ����� ��������� y' = \Psi_1(y) ���������� �����������
i=1;
while (t > 1e-1)
    Y(i+1) = Y(i) + (1 + A20*Y(i)/A10/(T(i))^2)*dT;
    T(i+1) = T(i) + dT;
    i=i+1;
    t = t + dT;
    % if (Y(i+1) >= 0) 
    %     break;
    % end
end

x_p = (T(1:i))';
y_p = (Y(1:i))';

Mas = [x_p y_p];

subplot(2,1,1);
hold on 
axis([0 1.1 -0.008 0])
grid on
plot(y_0(1),y_0(2),'o','MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',10);
legend('P_0');
plot(x_p,y_p,'Linewidth',1);
xlabel('y');
ylabel('\Psi_1(y)');

%--------------------------------------------------------------------------
subplot(2,1,2);
hold on
grid on
axis([0 1.1 -1 0.1])
plot(y_0(1),y_0(2),'o','MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',10);
plot(y_End(1),y_End(2),'s','MarkerEdgeColor','k',...
                       'MarkerFaceColor','g',...
                       'MarkerSize',10);
legend('P_0','P_*','Location','SouthEast');
plot(x_p,y_p,'Linewidth',2);
T=[];
Y=[];

T=y_0(1):dT:y_sing;
Y=T + y_0(2) - y_0(1);
len = length(T);

plot(T,Y,'Linewidth',2);
xlabel('y');
h=ylabel('\frac{dy}{d\tau}');
set(h,'Interpreter','tex')
%------------------------------���������-----------------------------------
di = (len-1)/10;
for i=1:di:len
    t=T(i):dT:y_sing;
    y_t = -t + Y(i) + T(i);
    for j=1:length(t)
        if (y_t(j)<0)
            plot(t(j),y_t(j));
        end
    end
end
%--------------------���������, �������� �����-----------------------------
