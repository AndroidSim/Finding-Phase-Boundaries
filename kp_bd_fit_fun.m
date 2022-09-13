function V = kp_bd_fit_fun(kp,c,alpha,beta) 

nc = length(c);
a = find(alpha == c);
b = find(beta == c);
V1 = ones(a-1,1);
V2 = (kp.*(beta-c(a:b)))./(c(a:b)-alpha+(kp.*(beta-c(a:b))));
V3 = zeros(nc-b,1);
V = [V1;V2';V3];