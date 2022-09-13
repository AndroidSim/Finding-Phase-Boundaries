function [spectra,baseline] = boundary_spectra_linprog(spectra,C)
% spectra should be derivative spectra from .dat files converted from ascii
% (.asc) files

[nb,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field intensity_values]');
end

ns = nc/2;

% convert derivative spectra to absorbance spectra
spectra = deriv2abs(spectra);

% normalize absorbance spectra
spectra = normalize_spectra(spectra,'a');

% cut off baseline and match length of cut spectra
[spectra,baseline] = cut_n_match(spectra,'a');

% get baseline statistics
for s = 1:ns
    avg_bline(s) = mean(baseline{s}(:,2));
    std_bline(s) = std(baseline{s}(:,2),1);
    highest_pos(s) = max(baseline{s}(baseline{s}(:,2) > 0,2));
    lowest_neg(s) = min(baseline{s}(baseline{s}(:,2) < 0,2));
end

% interpolate spectra for constant magnetic field increments
interval = (spectra(end,1) - spectra(1,1))/length(spectra(:,1));
interval = round(interval*1000)/1000;
Bi = [spectra(1,1):interval:spectra(end,1)]';
ispectra = zeros([length(Bi) nc]);

for s = 1:ns
    B = spectra(:,s*2-1);
    A = spectra(:,s*2);
    Ai = interp1(B,A,Bi); % default is linear
    ispectra(:,s*2-1:s*2) = [Bi Ai];
end

spectra = ispectra;

% normalize spectra so that the sum of absorbances equals one
for s = 1:ns
    spectra(:,s*2) = spectra(:,s*2)./sum(spectra(:,s*2));
end

% reduce spectra to just absorbance values and subtract mean spectrum
S = spectra(:,2:2:end);
S = S - repmat(mean(S,2),1,ns);

% perform singular value decomposition to obtain eigenspectra basis
[U,W,V] = svd(S,0);
B = U(:,1:2); % B = 2 component eigenbasis

% solve for eigencoefficient matrix D by linear programming
options = optimset('display','off');

% constraints
Aeq = [sum(U(:,1)) sum(U(:,2))];
beq = 1;
A = -B;
b = zeros(size(U(:,1)));
x0 = W(1:2,1:2)*V(:,1:2)';
% f = something
   
for s = 1:ns
    [x,fval] = linprog(f,A,b,Aeq,beq);
    D(:,s) = x;
end
   
return