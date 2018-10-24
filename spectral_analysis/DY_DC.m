function F = dy_dc(x,Y1,Y2,c1,c2,v)
% DY_DC:  equation with the following parameters to be solved:
%
% x(1) = c_alpha = composition at phase alpha boundary
% x(2) = c_beta = composition at phase beta boundary
% x(3) = f_alpha_1 = fraction of phase alpha at composition c1

ca = x(1);
cb = x(2);
fa = x(3);

dy_dc = dAbs_dX([Y1 Y2],[c1 c2],v);

%if (((c1-ca)^2)-((cb-c1)^2)+(cb-c1)*(cb-c2)) == 0
%    keyboard
%end
    
YB = ((cb-ca).*(Y1-fa.*Y2))./(((c1-ca)^2)-((cb-c1)^2)+(cb-c1)*(cb-c2));

YA = ((cb-ca).*(Y1-(1-fa).*Y2))./((cb-c1)*(cb-ca)-(c1-ca)*(cb-c2));

F = dy_dc-YB+YA;

F = [F; ca*fa+cb*(1-fa)-c1];