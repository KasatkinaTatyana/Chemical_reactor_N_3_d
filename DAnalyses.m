function [dUp dD]=DAnalyses()
global d k b y_0 y_End alpha betta
global g1 g2 b1 b0 a1 a0
d=0;
hd=0.01;
N_y=1000;
N_stop=10000;
eps = 1e-10;


dy=(y_End(1)-y_0(1))/N_y;
i=0;
Y=zeros(1,N_y+1);
dY=zeros(1,N_y+1);
while i<N_stop
    j=1;
    for y=y_0(1):dy:y_End(1)
        Y(j)=k*y+b+alpha*(y-y_0(1))*(y-y_End(1))+betta*(y-y_0(1))^2*(y-y_End(1))+...
                d*(y-y_0(1))^2*(y-y_End(1))^2;
        dY(j)=(k+alpha*(y-y_0(1))+alpha*(y-y_End(1))+2*betta*(y-y_0(1))*(y-y_End(1))+betta*(y-y_0(1))^2+...
                2*d*( (y-y_0(1))*(y-y_End(1))^2 + (y-y_End(1))*(y-y_0(1))^2 ) )*Y(j);
            
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
        
        % возвращение к исходным переменным
        T = InvStateDiff(V_y);
        CA = y0;
        CB = - y1;
        
        V_x = [CA CB T];
        
        u=( v / s(V_x) - g1*k1(T)*CA^2 - g2*k2(T)*CB...
            - a0 - a1*(T-273) ) / (b0+b1*(T-273));
        
        if (y == y_0(1))
            min_u = u;
            max_u = u;
        end
        if (u < min_u)
            min_u = u;
        end
        if (u > max_u)
            max_u = u;
        end
        j=j+1;
    end
    if (max(Y)>=-eps)||(min(dY-Y)<=eps)||(min(Y)<=-1+eps)||(max_u>1)||(min_u<0)
        break
    end
    d=d+hd;
    i=i+1;
end
dUp=d-hd;
    
d=0;
i=0;
while i<N_stop
    j=1;
    for y=y_0(1):dy:y_End(1)
        Y(j)=k*y+b+alpha*(y-y_0(1))*(y-y_End(1))+betta*(y-y_0(1))^2*(y-y_End(1))+...
                d*(y-y_0(1))^2*(y-y_End(1))^2;
        dY(j)=(k+alpha*(y-y_0(1))+alpha*(y-y_End(1))+2*betta*(y-y_0(1))*(y-y_End(1))+betta*(y-y_0(1))^2+...
                2*d*( (y-y_0(1))*(y-y_End(1))^2 + (y-y_End(1))*(y-y_0(1))^2 ) )*Y(j);
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
        
        % возвращение к исходным переменным
        T = InvStateDiff(V_y);
        CA = y0;
        CB = - y1;
        
        V_x = [CA CB T];
        
        u=( v / s(V_x) - g1*k1(T)*CA^2 - g2*k2(T)*CB...
            - a0 - a1*(T-273) ) / (b0+b1*(T-273));
        
        if (y == y_0(1))
            min_u = u;
            max_u = u;
        end
        if (u < min_u)
            min_u = u;
        end
        if (u > max_u)
            max_u = u;
        end
        j=j+1;
    end
    if (max(Y)>=-eps)||(min(dY-Y)<=eps)||(min(Y)<=-1+eps)||(max_u>1)||(min_u<0)
        break
    end
    d=d-hd;
    i=i+1;
end
dD=d+hd;