function [x,chisq] = diff_spec_fit(S,d)
% S = spectra = a (b x s) matrix, b = # of magnetic field point, s = # of samples
% spectra must be aligned and only contain the intensity values
% d = difference spectrum to be fitted, must be a vector

warning off
options = optimset('display','off');
[nb,ns] = size(S); % np = number of magnetic field points of each spectrum, ns = number of composition samples
%x = repmat(NaN,[ns ns]);
x = cell([ns ns]);
chisq = repmat(NaN,[ns ns]);

if all(size(d) == 1) | all(size(d) > 1) | ndims(d) > 2
    error('second argument must be a vector');
end

% B = basis set 
% p = coefficients (percent of probe in each phase)
% S(:,i) = spectrum
        
% choose basis
% B = [S(:,1) S(:,20)];

% get statistics on all possible difference spectra
%D = repmat(NaN,[nb ns*(ns-1)/2]);
%m = 1;
%for i = 1:ns-1 
%    for k = i+1:ns
%        D(:,m) = S(:,i) - S(:,k);
%        m = m+1;
%    end
%end

%avgb = mean(D,2);
% sigmab = std(D,1,2);
% sigmab(sigmab == 0) = 10^-3;
sigmab = repmat(0.5,[nb 1]);

for i = 1:ns-1 
    for k = i+1:ns       
        % choose basis
        %B = S(:,i) - S(:,k);
        B = [S(:,i) S(:,k)];
        
        % linear no constraints
        %x(i,k) = B\d;
        %chisq(i,k) = norm((x(i,k)*B - d)./sigmab)^2;
        %chisq(i,k) = 1 - ((x(i,k)*B)'*d)/(norm(x(i,k)*B)*norm(d));
        x{i,k} = B\d;
        chisq(i,k) = norm((B*x{i,k} - d)./sigmab)^2;
        %chisq(i,k) = norm([B*x{i,k} - d])^2;
        
        % linear least squares with constraints on p
        %[p,resnorm] = lsqlin(B,d,[],[],[],[],[],[],[],options);
        %x(i,k) = p;
        %chisq(i,k) = resnorm;
     
        % solution by svd
        %[U,W,V] = svd(B,0);
        %W = diag(W);
        %x(i,k) = ((U'*d)./W)*V;
        %chisq(i,k) = norm([B*x(i,k) - d])^2;
    end
end

x = x';
chisq = chisq';
warning on

return