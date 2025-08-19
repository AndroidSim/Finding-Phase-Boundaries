function B = find_boundaries(spectra,X,v,f,c)
% find_boundaries:  finds the compositional phase boundaries from the
% experimental esr spectra that vary with composition
%
% *assumption*: spectra are aligned 
%
% spectra = (number points x number of spectra) matrix, 
% X = (v x number of spectra) vector
% v = the number of composition variables
% f = degrees of freedom
% c = number of components

[np,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end

ns = nc/2;

[m,n] = size(X);

if m ~= v || n ~= ns
    error('X must be a row vector of size number variables x number of spectra');
end

% phase rule
p = c-f; % p = number of phases and thus number of boundaries

neval = ns-1;
B = zeros(neval,p+2);
s = [1 2];
options = optimset('largescale','on','display','iter','diagnostics','off','MaxFunEvals',2000);%,'NonlEqnAlgorithm','lm');
x0 = [0.001; 0.999; 0.5];

for i = 1:neval
    ci = s'*2;
    S1 = spectra(:,ci(1));
    S2 = spectra(:,ci(2));
    [x,fval,exitflag,output] = fsolve(@DY_DC,x0,options,S1,S2,X(s(1)),X(s(2)),v); 
    x0 = x;
    B(i,:) = [X(s(1)) x'];
    s = s+1;
end
    
return