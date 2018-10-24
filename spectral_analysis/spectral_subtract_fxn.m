function V = spectral_subtract_fxn(x,eS,rS,objfxn)
% x(1) = scaling factor
% x(2) = shift factor

% shift reference spectrum rS
iBrS = rS(:,1)+x(2);
% iBrS = rS(:,1)+s;

% interpolate to give intensity values
iIeS = interp1(eS(:,1),eS(:,2),iBrS);% linear interpolation 
iIeS(isnan(iIeS)) = 0;
eS = [iBrS iIeS];
rS(:,1) = iBrS;

% calculate residual spectrum: eS - scaling_factor*(shifted rS)
S(:,1) = eS(:,1);
% S(:,2) = eS(:,2) - x.*rS(:,2);
S(:,2) = eS(:,2) - x(1).*rS(:,2);

switch objfxn
    case 'area'
        % calculate area of residual spectrum
		dH = zeros(size(S(:,1)));
		dH(2:end) = diff(S(:,1));
		A = abs(sum(prod([dH S(:,2)],2))); 
    case 'var2der'
        % calculate variance of 2nd deriv of residual spectrum
        % method from paper:  
        %   Loethen et al. "Second-Derivative Variance Minimization Method
        %   for Automated Spectral Subtraction".  Applied Spectroscopy
        %   58(3) 2004. p.272-278.
		S = abs2deriv(S); % actually taking second derivative of absorption
        % Savitsky-Golay smoothing
        noise = max(S(1:100,2))-min(S(1:100,2));
        signal = max(S(:,2))-min(S(:,2));
        S(:,2) = smooth(S(:,2),(noise/signal)*length(S),'sgolay',2);
		V = std(S(:,2),1)^2;
    case 'skewness'
        % calculate the absolute value of the skewness of absorbance spectrum
        % from "Electron Spin Resonance: A Comprehensive Treatise on
        % Experimental Techniques" Charles P. Poole, Jr.
        S = deriv2abs(S);
        S(S(:,2) < 0,2) = 0;
		H = S(600:900,1);
		y = S(600:900,2);
        [maxy,imax] = max(y);
		H0 = H(imax);
        dH = zeros(size(H));
		dH(2:end) = diff(H);
		A = sum(prod([dH y],2));
% 		B = (sum(prod([(H-H0) y],2)))/(sum(y));
% 		H0 = H0+B;
		H3 = (1/A)*sum(prod([dH (H-H0).^3 y],2));
		H2 = (1/A)*sum(prod([dH (H-H0).^2 y],2));
		V = H3/(H2^(3/2));
        V = abs(V);
    otherwise
        error('invalid objective fxn');
end

return