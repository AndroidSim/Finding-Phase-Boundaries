function abs_spec = Abs_spectra(spectra,norm_flag)
% Abs_spectra calculates the absorbance spectrum from the derivative
% spectrum given in spectra.

[np,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end

normalize = false;

if nargin == 2 & isnumeric(norm_flag)
    normalize = true;
end

ns = nc/2;
abs_spec = zeros(size(spectra));
abs_spec = spectra;

for i = 1:ns
    % noise = the noise in the spectrum (should be passed in)
    % llim_noise = lower limit of noise
    llim_noise = min(spectra(1:200,i*2)); % or noise = mean(spectra(1:200,i*2)) ?
    
    % check to see if spectra are derivative or absorbance spectra
    if ~any(spectra(:,i*2) < llim_noise)
        error(sprintf('The %d-th spectrum with values in column %d is already an absorbance spectrum',i,i*2));
    end
    
    % calculate absorbance spectrum
    abs_spec1(1,i*2) = 0;
    
    for j=2:np
        abs_spec(j,i*2) = abs_spec(j-1,i*2)+spectra(j,i*2);
    end
        
    % baseline correction
    a1 = spectra(1,i*2-1);
    b1 = abs_spec(1,i*2);
    a2 = spectra(end,i*2-1);
    b2 = abs_spec(end,i*2);
    A = (b2-b1)/(a2-a1); 
    B = b2-a2*A;
    abs_spec(:,i*2) = abs_spec(:,i*2)-(A*spectra(:,i*2-1)+B);
    
    if normalize
        % normalize absorbance spectrum
        abs_spec(:,i*2) = abs_spec(:,i*2)/sum(abs_spec(:,i*2))/mean(diff(abs_spec(:,i*2-1)));
        % mean(diff(abs_spec(:,i*2-1))) == mean(abs_spec(2:end,i*2-1)-abs_spec(1:end-1,i*2-1))
    end
    
    %spca1(sm:sn,2)=spca1(sm:sn,2)/sum(spca1(sm:sn,2))/mean...
    %(spca1((sm+1):sn,1)-spca1(1:(sn-1),1));%Normalize ABS spectrum

    %a1=spc(1,1);b1=spca1(1,2);a2=spc(sn,1);b2=spca1(sn,2);
    %A=(b2-b1)/(a2-a1); B=b2-a2*A;
    %spca1(sm:sn,2)=spca1(sm:sn,2)-(A*spc(sm:sn,1)+B);
end
return