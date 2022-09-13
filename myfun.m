function F = myfun(kp,c,alpha,beta) 

nc = length(c);
a = find(alpha == c);
b = find(beta == c);
F1 = ones(a-1,1);
F2 = (kp.*(beta-c(a:b)))./(c(a:b)-alpha+(kp.*(beta-c(a:b))));
F3 = zeros(nc-b,1);
F = [F1;F2';F3];