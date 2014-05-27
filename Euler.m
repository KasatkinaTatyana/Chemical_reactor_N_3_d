
function F=Euler(t, h_t)
% h_t - шаг интегрирования
global y_0 y_End k b alpha betta t_0
F=y_0(1);

if (t > t_0)
    y=y_0(1);
    for dt=t_0 : h_t : t
        y=y+h_t*(k*y+b+alpha*(y-y_0(1))*(y-y_End(1))+betta*(y-y_0(1))^2*(y-y_End(1)));
    end
    F=y;
end
