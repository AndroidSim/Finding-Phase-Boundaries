function spectra = preprocess_spectra(spectra,freq,X,varargin)
% preprocess_spectra: converts raw derivative spectra in ascii format to
% absorbance spectra and then scales, normalizes, corrects and transforms, removes baseline, and interpolates the
% spectra into a form used for analysis

[np,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end

ns = nc/2;

if length(freq) ~= ns
    error('number of spectral frequencies must = number of spectra');
end

% convert to absorbance spectra
spectra = deriv2abs(spectra);

% scale so that maximum absorbance value == 100
for i = 1:ns
    spectra(:,i*2) = (spectra(:,i*2)/max(spectra(:,i*2)))*100;
end

% correct for different frequencies, i.e. variations in g and A tensors
%
% this approximately(?) amounts to shifting magnetic fields to mean magnetic field 
% calculated from the mean frequency
B0 = ((mean(freq.*10^9))*(6.6262*10^-34))/(2.0023*(9.2740*10^-28));
B0 = num2str(B0,'%10.6g');
B0 = str2num(B0);
[A,I] = max(spectra(:,2:2:end));

for i = 1:ns
    B(i) = spectra(I(i),i*2-1);
end
%[B,indices] = peak(spectra,B0);

for i = 1:ns
    spectra(:,i*2-1) = spectra(:,i*2-1)-(B(i)-B0);
end

% normalize absorbance spectra so that total area = 1
spectra = normalize_spectra(spectra,'a');

% transform magnetic field 
spectra = transform_B(spectra);

% if desired, remove baseline, align, and fill spectra 
%if nargin > 3   
%    spectra = interp_spectra(spectra,'a',X);
%end

return