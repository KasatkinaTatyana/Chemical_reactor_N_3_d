function F=y4Polynom(y)
global k b alpha betta y_0 y_End d

F=1./(k*y+b+alpha*(y-y_0(1)).*(y-y_End(1))+betta*((y-y_0(1)).^2).*(y-y_End(1))+...
    d*((y-y_0(1)).^2).*((y-y_End(1)).^2));