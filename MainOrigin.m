% Исходная система с исходным управлением
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
%-------------------------Константы реактора-------------------------------
R=8.3;
E1=2.09e4;
E2=4.18e4;
A10=1.1;
A20=172.2; 
%----------------------Коэффициенты управления-----------------------------
a0=4.3145;
a1=-0.1099;
b0=1.4962;
b1=0.0515;
g1=41.8;
g2=83.6;
%--------------------------------------------------------------------------
global N    %число разбиений
N=50000;
t_0=0;
%--------------------------------------------------------------------------


x_0=[0.993 0.007 298];
% x_End=[0.35 0.55 333];

x_End=[0.584 0.4 333];
%--------------------------------------------------------------------------
%Граничные значения канонических переменных
y_0=YTrans(x_0);
y_End=YTrans(x_End);
%--------------------------------------------------------------------------
k=(y_End(2)-y_0(2))/(y_End(1)-y_0(1));
b=(y_0(2)*y_End(1)-y_End(2)*y_0(1))/(y_End(1)-y_0(1));
alpha = (y_0(3)-k*y_0(2))/y_0(2)/(y_0(1)-y_End(1));
betta = (y_0(2)*y_End(3)+y_0(3)*y_End(2)-2*k*y_0(2)*y_End(2))/y_0(2)/y_End(2)/(y_End(1)-y_0(1))^2;
%%
% [dPos dNeg]=DFind();
%%
[dUp dD]=DAnalyses();


d=0.09;
%%
dy=(y_End(1)-y_0(1))/N;

tau_End=quad(@y4Polynom,y_0(1),y_End(1));
dtau=tau_End/N;
%% Определение t_End
t_End=t_0;
for y=y_0(1):dy:y_End(1)
    y0=y;
    y1=k*y0+b+alpha*(y0-y_0(1))*(y0-y_End(1))+betta*(y0-y_0(1))^2*(y0-y_End(1))+...
        d*(y0-y_0(1))^2*(y0-y_End(1))^2;

    y2=(k+alpha*(y0-y_0(1))+alpha*(y0-y_End(1))+2*betta*(y0-y_0(1))*(y0-y_End(1))+betta*(y0-y_0(1))^2+...
        2*d*( (y0-y_0(1))*(y0-y_End(1))^2 + (y0-y_End(1))*(y0-y_0(1))^2 ) )*y1;

    y3=2*(alpha+3*betta*y0-2*betta*y_0(1)-betta*y_End(1)+...
        2*d*( (y0-y_0(1))^2 + 4*(y0-y_0(1))*(y0-y_End(1)) + (y0-y_End(1))^2 ))*y1^2+...
        (k+alpha*(y0-y_0(1))+alpha*(y0-y_End(1))+2*betta*(y0-y_0(1))*(y0-y_End(1))+betta*(y0-y_0(1))^2+...
        2*d*( (y0-y_0(1))*(y0-y_End(1))^2 + (y0-y_End(1))*(y0-y_0(1))^2 ) )^2*y1;

    V_y=[y0 y1 y2];
    
    T=InvStateDiff(V_y);
    t_End=t_End-dy/(k1(T)*y^2);
end
t_End

%% Моделирование исходной системы в реальном времени t
h_t = t_End/N;

x_r=x_0; tau=t_0;
i=1;

c=0.1;

% x_r=zeros(N,3);
% mas_u=zeros(N);

for t=t_0:h_t:t_End-h_t
    y_r=YTrans(x_r(i,:));
    
    U=Control(tau);

    y_tau=U(4:6,1);
    
    f_r=fKan(y_r);
    g_r=gKan(y_r);

    v=(-f_r+U(7,1) - 3*c*(y_r(3) - y_tau(3)) - 3*c^2*(y_r(2) - y_tau(2)) - c^3*(y_r(1) - y_tau(1)))/g_r;    % стабилизирующее управление
    % u - реальное управление, которое должно лежать в диапазоне [0; 1]
    u=( v / s(x_r(i,:)) - g1*k1(x_r(i,3))*x_r(i,1)^2 - g2*k2(x_r(i,3))*x_r(i,2)...
        - a0 - a1*(x_r(i,3)-273) ) / (b0+b1*(x_r(i,3)-273));
    
    mas_u(i) = u;
    
%     if (u < 0)
%         u = 0;
%     end
%     if (u > 1)
%         u = 1;
%     end
    
    h_tau = h_t / s(x_r(i,:));
    
    x_r(i+1,1)=x_r(i,1)+(-k1(x_r(i,3))*x_r(i,1)^2)*h_t;
    x_r(i+1,2)=x_r(i,2)+( k1(x_r(i,3))*x_r(i,1)^2 - k2(x_r(i,3))*x_r(i,2) )*h_t;
    x_r(i+1,3)=x_r(i,3)+( g1*k1(x_r(i,3))*x_r(i,1)^2+g2*k2(x_r(i,3))*x_r(i,2)...
        +a0+a1*(x_r(i,3)-273) + (b0+b1*(x_r(i,3)-273))* u )* h_t;
    i=i+1;
    tau=tau+h_tau;
end

disp(x_r(end,:));
disp(U(4:6,1));

tx = t_0:h_t:t_End;

% figure(1);
% title('Зависимость концентрации C_A от времени t');
% plot(tx',x_r(:,1));
% xlabel('t, c');
% ylabel('C_A');
% 
% figure(2);
% title('Зависимость концентрации C_B от времени t');
% plot(tx',x_r(:,2));
% xlabel('t, c');
% ylabel('C_B');
% 
% figure(3);
% title('Зависимость температуры T от времени t');
% plot(tx',x_r(:,3));
% xlabel('t, c');
% ylabel('T, K');
% 
% figure(4);
% title('Зависимость управления u от времени t');
% txu = t_0:h_t:t_End-h_t;
% plot(txu',mas_u);
% xlabel('t, c');
% ylabel('u');               % в обозначениях статьи u это v - исходное
% управление

subplot(2,2,1); 
title('Зависимость концентрации C_A от времени t');
plot(tx',x_r(:,1));
xlabel('t, c');
ylabel('C_A');

subplot(2,2,2);
title('Зависимость концентрации C_B от времени t');
plot(tx',x_r(:,2));
xlabel('t, c');
ylabel('C_B');

subplot(2,2,3);
title('Зависимость температуры T от времени t');
plot(tx',x_r(:,3));
xlabel('t, c');
ylabel('T, K');

subplot(2,2,4);
title('Зависимость управления v от времени t');
txu = t_0:h_t:t_End-h_t;
plot(txu',mas_u);
xlabel('t, c');
ylabel('v');