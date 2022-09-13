function V = boundary_spectra_fun(x,c)
%x = [kp a_alpha a_beta c_alpha c_beta]

%F = (a-x(3)).*(x(2)-x(4)) - (b-x(4)).*(x(1)-x(3));

nc = length(c);
V = repmat(0,nc,1);

%for i = 1:nc
%    if c(i) <= x(3)
%        F(i) = x(1);
%    elseif c(i) > x(3) & c(i) < x(4)
%        F(i) = ((x(1)-x(2))/(x(3)-x(4)))*(c(i)-x(3)) + x(1);
%    elseif c(i) >= x(4)
%        F(i) = x(2);
%    end
%end

%for i = 1:nc
%    F(i) = ((x(1).*(x(5)-c(i)))./(c(i)-x(4)+(x(1).*(x(5)-c(i))))).*(x(2)-x(3)) + x(3);
%end

for i = 1:nc
    if c(i) <= x(4)
        V(i) = x(2);
    elseif c(i) > x(4) & c(i) < x(5)
        V(i) = ((x(1).*(x(5)-c(i)))./(c(i)-x(4)+(x(1).*(x(5)-c(i))))).*(x(2)-x(3)) + x(3);
    elseif c(i) >= x(5)
        V(i) = x(3);
    end
end

return