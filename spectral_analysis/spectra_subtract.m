function S = spectra_subtract(eS,rS)
% spectra_subtract subtracts reference spectrum rS from experimental spectrum eS and returns the
% leftover spectrum S

[nb,nc] = size(eS);

if nc < 2 || rem(nc,2) ~= 0
    error('spectrum 1 must contain 2 columns: [B-field absorbance_values]');
end

[nb,nc] = size(rS);

if nc < 2 || rem(nc,2) ~= 0
    error('spectrum 2 must contain 2 columns: [B-field absorbance_values]');
end

% align spectra
aS = align_spectra([eS rS]);
eS = aS(:,1:2);
rS = aS(:,3:4);

options = optimset('display','off');
x = fminunc(@spectral_subtract_fxn,[0.5;0],options,eS,rS); 
% shift reference spectrum rS
iBrS = rS(:,1)+x(2);
% interpolate to give new intensity values
iIrS = interp1(rS(:,1),rS(:,2),iBrS);% linear interpolation 
rS = [iBrS iIrS];
% calculate residual spectrum: eS - scaling_factor*(shifted rS)
S(:,1) = eS(:,1);
S(:,2) = eS(:,2) - x(1).*rS(:,2);

% % n0 = norm(IS1);
% n0 = sum(IS1);
% for s = -ls2:ls2
%     is2 = shift_vector(IS2,s);
%     x = is2\IS1;
% %     n = norm(IS1-x.*is2);
%     n = sum(IS1-x.*is2);
%     if n <= n0
%         n0 = n;
%         IS = is2;
%         p = x;
%     end
% end
% 
% S = S1;
% S(:,2) = IS1-p.*IS;

return