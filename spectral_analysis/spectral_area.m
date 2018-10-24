function A = spectral_area(spectra)

[np,ncol] = size(spectra);
if ncol < 2 || rem(ncol,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end
ns = ncol/2;

for s = 1:ns
%     % get baseline statistics
%     baseline(:,s) = spectra([1:100 end-100:end]',s*2);
%     avg_bline(s) = mean(baseline(:,s));
%     std_bline(s) = std(baseline(:,s),1);
%     
%     highest_pos(s) = max(baseline(:,s));
%     if isempty(highest_pos(s))
%         highest_pos(s) = 0;
%     end
%     
%     lowest_neg(s) = min(baseline(:,s));
%     if isempty(lowest_neg(s))
%         lowest_neg(s) = 0;
%     end
% 
%     % determine if spectrum is derivative and if so convert to absorbance
%     if any(spectra(:,s*2) < lowest_neg(s))
%         % spectrum is a derivative spectrum
%         spectra(:,s*2-1:s*2) = deriv2abs(spectra(:,s*2-1:s*2));
%     end
    
    % calculate area of absorbance spectrum
	dH = zeros(size(spectra(:,s*2-1)));
	dH(2:end) = diff(spectra(:,s*2-1));
	A(s,1) = sum(prod([dH spectra(:,s*2)],2));
end

return