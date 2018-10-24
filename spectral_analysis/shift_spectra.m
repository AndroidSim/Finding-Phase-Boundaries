function spectra = shift_spectra(spectra,varargin)
% shift_spectra: shifts the B field of the maximum intensity by a given amount
%   spectra = (number magnetic fields x 2*number of spectra)
%   shift = vector or scalar specifying the amount to shift

[nb,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end

ns = nc/2;

% from asc2dat.m:
% [c1,i1]=max(spc(:,2));
% spc(:,1)=spc(:,1)-(spc(i1,1)-B0);

if nargin > 1
    switch varargin{1}
        case 'B0'
            if isscalar(varargin{2})
                B0 = varargin{2};
            else
                error('3rd argument must be a scalar magnetic field value, B0');
            end
            [mi,ii] = max(spectra(:,2:2:end));
            for s = 1:ns
                spectra(:,s*2-1) = spectra(:,s*2-1)-(spectra(ii(s),s*2-1)-B0);
            end
        case 'inc'
            if isscalar(varargin{2}) || isvector(varargin{2})
                shift = varargin{2};
            else
                error('3rd argument must be a scalar or vector of shift incrments');
            end
            nshifts = length(shift);
            if nshifts == 1 | nshifts == ns
                for s = 1:ns
                    spectra(:,s*2-1) = spectra(:,s*2-1)+shift(s);
                end
            else
                error('the length of shift vector must be equal to the number of spectra or 1');
            end
        otherwise
            error('2nd argument must be "B0" or "inc" specifying shift method');
    end
else
    % shift about magnetic field of nearest peak to the median magnetic field
    medB = median(median(spectra(:,1:2:end)));
%     [mi,ii] = max(spectra(:,2:2:end));
%     [H,ind] = peak(spectra,medB);
    abs_spectra = deriv2abs(spectra);
    [mi,ii] = max(abs_spectra(:,2:2:end));
%     [H,ind] = peak(abs_spectra,medB);
    for s = 1:ns
        H = abs_spectra(ii(s),s*2-1);
        dB = H-medB;
        spectra(:,s*2-1) = spectra(:,s*2-1)-dB;
    end
%     dB = H-medB;
%     for s = 1:ns   
%         spectra(:,s*2-1) = spectra(:,s*2-1)-dB(s); 
%         % spectra(:,s*2-1) = spectra(:,s*2-1)-(spectra(ii(s),s*2-1)-round(avgB));
%     end
end

% correct for different frequencies
% this approximately(?) amounts to shifting magnetic fields to mean magnetic field 
% calculated from the mean frequency
% mB = ((freq_mixS.*10^9)*(6.6262*10^-34))/(2.0023*(9.2740*10^-28));
% rB = ((freq_refS.*10^9)*(6.6262*10^-34))/(2.0023*(9.2740*10^-28));
% B0 = ((mean([freq_mixS;freq_refS].*10^9))*(6.6262*10^-34))/(2.0023*(9.2740*10^-28));
% B0 = num2str(B0,'%10.6g');
% B0 = str2num(B0);
% mixS = shift_spectra(mixS,B0-mB);
% refS = shift_spectra(refS,B0-rB);
% if length(freq_refS) ~= nrs
%     error('number of frequencies of buffer spectra must = number of buffer spectra');
% end
    
return