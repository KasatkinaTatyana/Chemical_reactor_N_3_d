function dy=Traject(t,y)
global y_0 y_End k b alpha betta d
dy=zeros(1,1);

dy(1)=k*y+b+alpha*(y-y_0(1))*(y-y_End(1))+betta*(y-y_0(1))^2*(y-y_End(1))+...
    d*(y-y_0(1))^2*(y-y_End(1))^2;