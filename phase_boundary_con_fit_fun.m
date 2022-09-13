function [V,D] = phase_boundary_con_fit_fun(x,S,C,e_spec)
% x = [Kp C_alpha C_beta]
% S = spectra matrix
% C = spectal composition vector
% e_spec = experimental or data spectra vector (each spectrum stacked on
% top of one another)

[nb,ns] = size(S);
nc = ns;

temp = [x(2) x(3)]';
low = temp < C(1);
up = temp > C(end);

if any(low)
    temp(low) = C(1);
end

if any(up)
    temp(up) = C(end);
end

x(2:3) = temp;

% check compositions with C_alpha and C_beta to insure no generation of NaN
for c = 1:nc
    if x(3) == (C(c)*(x(1)-1)+x(2))/x(1)
        if x(2) == x(3) & x(2) == C(1)
            x(3) = x(2)+0.0001;
        elseif x(2) == x(3) & x(2) == C(end)
            x(2) = x(3)-0.0001;
        elseif x(2) == x(3) & x(2) == C(c)
            x(3) = x(2)+0.0001;
        else
            x(1) = x(1)+0.0001;
        end
    else
        if x(2) == x(3) & x(2) == C(1)
            x(3) = x(2)+0.0001;
        elseif x(2) == x(3) & x(2) == C(end)
            x(2) = x(3)-0.0001;
        elseif x(2) == x(3) & x(2) == C(c)
            x(3) = x(2)+0.0001;
        end
    end
end

% interpolate spectra matrix S to find the basis spectra at x(2) = alpha
% and x(3) = beta
Ci = [x(2) x(3)];

for b = 1:nb
    A = S(b,:);
    Ai = interp1(C,A,Ci); % default = linear
    B(b,:) = Ai;
end

% generate theoretical spectra from basis B
for c = 1:nc
    t_spec(:,c) = ((x(1).*(x(3)-C(c)))./(C(c)-x(2)+(x(1).*(x(3)-C(c))))).*B(:,1) + ((C(c)-x(2))./(C(c)-x(2)+(x(1).*(x(3)-C(c))))).*B(:,2);
end
        
t_spec = reshape(t_spec,size(e_spec));

%if any(isnan(t_spec))
%    keyboard
%end

V = sum(0.5.*(e_spec - t_spec).^2);
%V = 1 - dot(e_spec,t_spec)/(norm(e_spec)*norm(t_spec));

if nargout > 1
    D = [0 0 0];
    
    for s = 1:ns
        % calculate gradient D
        kp = x(1);
        a = x(2);
        b = x(3);
        c = C(s);
        Sa = B(:,1);
        Sb = B(:,2);
        dV = [((b-c)/(c-a+kp*(b-c))).*Sa-(kp*(b-c)^2/(c-a+kp*(b-c))^2).*Sa-((c-a)*(b-c)/(c-a+kp*(b-c))^2).*Sb,...             
            (kp*(b-c)/(c-a+kp*(b-c))^2).*Sa-(1/(c-a+kp*(b-c))).*Sb+((c-a)/(c-a+kp*(b-c))^2).*Sb,...      
            (kp/(c-a+kp*(b-c))).*Sa-(kp^2*(b-c)/(c-a+kp*(b-c))^2).*Sa-((c-a)*kp/(c-a+kp*(b-c))^2).*Sb];
        D = D+sum(dV);
        %D = D+(S(:,s)'*dV);
    end
    
    %D = -D'./(norm(e_spec)*norm(t_spec));
    D = D';
 end
 
 return

% syms kp a b c Sa Sb;
% S = ((kp*(b-c))/(c-a+(kp*(b-c))))*Sa + ((c-a)/(c-a+(kp*(b-c))))*Sb;
% dS_kp =
 
% (b-c)/(c-a+kp*(b-c))*Sa-kp*(b-c)^2/(c-a+kp*(b-c))^2*Sa-(c-a)/(c-a+kp*(b-c))^2*Sb*(b-c)
 
% dS_a =
 
% kp*(b-c)/(c-a+kp*(b-c))^2*Sa-1/(c-a+kp*(b-c))*Sb+(c-a)/(c-a+kp*(b-c))^2*Sb
 
% dS_b =
 
% kp/(c-a+kp*(b-c))*Sa-kp^2*(b-c)/(c-a+kp*(b-c))^2*Sa-(c-a)/(c-a+kp*(b-c))^2*Sb*kp

% D = [ (b-c)/(c-a+kp*(b-c))*Sa-kp*(b-c)^2/(c-a+kp*(b-c))^2*Sa-(c-a)/(c-a+kp*(b-c))^2*Sb*(b-c),             
%        kp*(b-c)/(c-a+kp*(b-c))^2*Sa-1/(c-a+kp*(b-c))*Sb+(c-a)/(c-a+kp*(b-c))^2*Sb,       
%        kp/(c-a+kp*(b-c))*Sa-kp^2*(b-c)/(c-a+kp*(b-c))^2*Sa-(c-a)/(c-a+kp*(b-c))^2*Sb*kp]