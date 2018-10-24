function A = spectra_subtract_fxn(x,eS,rS)
% x(1) = scaling factor
% x(2) = shift factor

% shift reference spectrum rS
iBrS = rS(:,1)+x(2);

% interpolate to give new intensity values
iIrS = interp1(rS(:,1),rS(:,2),iBrS);% linear interpolation 
rS = [iBrS iIrS];

% calculate residual spectrum: eS - scaling_factor*(shifted rS)
S(:,1) = eS(:,1);
S(:,2) = eS(:,2) - x(1).*rS(:,2);

% calculate area of residual spectrum
dH = zeros(size(S(:,1)));
dH(2:end) = diff(S(:,1));
A = sum(prod([abs(dH) S(:,2)],2)); 

return