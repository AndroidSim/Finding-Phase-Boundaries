function varargout = spectral_subtract(eS,rS,varargin)
% spectra_subtract subtracts reference spectrum rS from experimental spectrum eS and returns the
% leftover spectrum S

[nb,nc] = size(eS);

if nc < 2 || rem(nc,2) ~= 0
    error('spectrum 1 must contain 2 columns: [B-field intensity_values]');
end

[nb,nc] = size(rS);

if nc < 2 || rem(nc,2) ~= 0
    error('spectrum 2 must contain 2 columns: [B-field intensity_values]');
end

if isempty(varargin)
    isauto = true;
    isgrid = false;
    objfxn = 'var2der';
else
    if ischar(varargin{1})
        if strcmp(varargin{1},'grid')
            isauto = true;
            isgrid = true;
            if isempty(varargin{2})
                objfxn = 'var2der';
            else
                if any(strcmp(varargin{2},{'var2der','skewness'}))
                    objfxn = varargin{2};
                else
                    error('objective fxn must be the strings var2der or skewness');
                end
            end
        else
            isgrid = false;
            if isempty(varargin{1})
                objfxn = 'var2der';
            else
                if any(strcmp(varargin{1},{'var2der','skewness'}))
                    isauto = true;
                    objfxn = varargin{1};
                else
                    error('objective fxn must be the strings var2der or skewness');
                end
            end
        end
    elseif isnumeric(varargin{1})    
        isauto = false;
        if any(size(varargin{1}) == 1) %isvector
            p = varargin{1};
        else
            error('var arg in must be a vector = [scaling shifting]');
        end
    else
        error('1st var arg in must be a string or numeric');
    end
end

% align spectra
aS = align_spectra([eS rS]);
eS = aS(:,1:2);
rS = aS(:,3:4);

if isauto
    if isgrid
		fval0 = feval(@spectral_subtract_fxn,[0 0],eS,rS,objfxn);
		% start grid search
		for shift = -5:0.1:5       
		    for scale = 0:0.01:1
		        fval = feval(@spectral_subtract_fxn,[scale shift],eS,rS,objfxn);
                if fval <= fval0
                    fval0 = fval;
                    p = [scale shift];
                end
            end
	    end
    else
        options = optimset('display','off','TolX',10^-8,'TolFun',10^-8);
        p = fminsearch(@spectral_subtract_fxn,[0.5;0],options,eS,rS,objfxn);
    end
    
end

% calculate subtracted difference spectra with best parameters
% shift reference spectrum rS
iBrS = rS(:,1)+p(2);
% interpolate to give intensity values
iIeS = interp1(eS(:,1),eS(:,2),iBrS);% linear interpolation 
iIeS(isnan(iIeS)) = 0;
eS = [iBrS iIeS];
rS(:,1) = iBrS;

% calculate residual spectrum: eS - scaling_factor*(shifted rS)
S(:,1) = eS(:,1);
S(:,2) = eS(:,2) - p(1).*rS(:,2);
bS = rS;
bS(:,2) = p(1).*rS(:,2);

if isauto
    varargout{1} = S;
    varargout{2} = bS;
    varargout{3} = p;
else
    varargout{1} = S;
    varargout{2} = bS;
    varargout{3} = std(S(:,2),1)^2;
end

return