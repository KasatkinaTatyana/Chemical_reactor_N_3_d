function [dUp dD]=DAnalyses()
global d k b y_0 y_End alpha betta
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
        j=j+1;
    end
    if (max(Y)>=-eps)||(min(dY-Y)<=eps)
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
        j=j+1;
    end
    if (max(Y)>=-eps)||(min(dY-Y)<=eps)
        break
    end
    d=d-hd;
    i=i+1;
end
dD=d+hd;