% ������ ��������� ����������� ��� ����� ��������� d, ������� � ���������
% �� dD �� dUp �������� ���������� ������� ����� � ��������� [0; 1].

close all
clear all
clc
global x_0 x_End
global y_0 y_End
global t_0
global R E1 E2 A10 A20
global k b alpha betta
global d
global g1 g2 b1 b0 a1 a0
%% -------------------------��������� ��������------------------------------
R=8.3;
E1=2.09e4;
E2=4.18e4;
A10=1.1;
A20=172.2; 
%% ----------------------������������ ����������-----------------------------
a0=4.3145;
a1=-0.1099;
b0=1.4962;
b1=0.0515;
g1=41.8;
g2=83.6;
%% --------------------------------------------------------------------------
t_0=0;
h_t=0.05;       %��� �� ������� t
c_st=[1 3 3];
%% --------------------------------------------------------------------------
global N    %����� ���������
N=1000;
N_tau=100;

x_0=[0.993 0.007 298];
x_End=[0.35 0.55 333];

x_End=[0.693 0.3 333];
%% --------------------------------------------------------------------------
%��������� �������� ������������ ����������
y_0=YTrans(x_0);
y_End=YTrans(x_End);
%% --------------------------------------------------------------------------
k=(y_End(2)-y_0(2))/(y_End(1)-y_0(1));
b=(y_0(2)*y_End(1)-y_End(2)*y_0(1))/(y_End(1)-y_0(1));
alpha = (y_0(3)-k*y_0(2))/y_0(2)/(y_0(1)-y_End(1));
betta = (y_0(2)*y_End(3)+y_0(3)*y_End(2)-2*k*y_0(2)*y_End(2))/y_0(2)/y_End(2)/(y_End(1)-y_0(1))^2;
%% --------------------------------------------------------------------------
d = 0;
[dUp dD]=DAnalyses();
%% ----------------------------Test DAnalyses--------------------------------
N_y=1000;
dy=(y_End(1)-y_0(1))/N_y;
YUp=zeros(1,N_y+1);
dYUp=zeros(1,N_y+1);
YD=zeros(1,N_y+1);
dYD=zeros(1,N_y+1);
j=1;
for y=y_0(1):dy:y_End(1)
    YUp(j)=k*y+b+alpha*(y-y_0(1))*(y-y_End(1))+betta*(y-y_0(1))^2*(y-y_End(1))+...
                dUp*(y-y_0(1))^2*(y-y_End(1))^2;
    dYUp(j)=(k+alpha*(y-y_0(1))+alpha*(y-y_End(1))+2*betta*(y-y_0(1))*(y-y_End(1))+betta*(y-y_0(1))^2+...
                2*dUp*( (y-y_0(1))*(y-y_End(1))^2 + (y-y_End(1))*(y-y_0(1))^2 ))*YUp(j);
    YD(j)=k*y+b+alpha*(y-y_0(1))*(y-y_End(1))+betta*(y-y_0(1))^2*(y-y_End(1))+...
                dD*(y-y_0(1))^2*(y-y_End(1))^2;
    dYD(j)=(k+alpha*(y-y_0(1))+alpha*(y-y_End(1))+2*betta*(y-y_0(1))*(y-y_End(1))+betta*(y-y_0(1))^2+...
                2*dD*( (y-y_0(1))*(y-y_End(1))^2 + (y-y_End(1))*(y-y_0(1))^2 ))*YD(j);
            
    j=j+1;
end
disp('max(YUp)=');
disp(max(YUp));
disp('max(YD)=');
disp(max(YD));
disp('min(dYUp-YUp)=');
disp(min(dYUp-YUp));
disp('min(dYD-YD)=');
disp(min(dYD-YD));
%% -------------------------------------------------------------------------
dy=(y_End(1)-y_0(1))/N;
%% ����������� t_End
for d=dD:0.01:dUp
    t_End = t_0;
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
        
        f = fKan(V_y);
        g = gKan(V_y);
        
        v = (y3-f)/g;
        
        % ����������� � �������� ����������
        T = InvStateDiff(V_y);
        CA = y0;
        CB = - y1;
        
        V_x = [CA CB T];
        
        u=( v / s(V_x) - g1*k1(T)*CA^2 - g2*k2(T)*CB...
            - a0 - a1*(T-273) ) / (b0+b1*(T-273));
        
        if (y == y_0(1))
            min_u = u;
            max_u = u;
            
            min_T = T;
            max_T = T;
        else
            if (u < min_u)
                min_u = u;
            end
            if (u > max_u)
                max_u = u;
            end
            
            if (T < min_T)
                min_T = T;
            end
            if (T > max_T)
                max_T = T;
            end
        end
        
        t_End=t_End-dy/(k1(T)*y^2);
    end% y
    
    
    
    
    if (min_u>=0)&&(max_u<=1)
        d
    end
    
    delta_u = max_u - min_u;
    delta_T = max_T - min_T;
    
    if (d == dD)
        min_delta_u = delta_u;
        min_u_d = dD;
        min_delta_T = delta_T;
        min_T_d = dD;
        min_t = t_End;
        min_t_d = dD;
    else
        if (min_delta_u > delta_u)
            min_delta_u = delta_u;
            min_u_d = d;
        end
        
        if (min_delta_T > delta_T)
            min_delta_T = delta_T;
            min_T_d = d;
        end
        
        if (min_t > t_End)
            min_t = t_End;
            min_t_d = d;
        end
    end
       
end

disp('���������� �������� max_u - min_u = ');
disp(min_delta_u);
disp('d = ');
disp(min_u_d);

disp('���������� �������� max_T - min_T = ');
disp(min_delta_T);
disp('d = ');
disp(min_T_d);

disp('���������� �������� ������� t = ');
disp(min_t);
disp('d = ');
disp(min_t_d);