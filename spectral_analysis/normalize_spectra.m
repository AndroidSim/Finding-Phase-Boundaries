function spectra = normalize_spectra(spectra,aord,varargin)
% normalizes each spectrum in spectra
%
% aord = 'a' for absorbance or 'd' for derivative

[np,nc] = size(spectra);
if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end
ns = nc/2;

if isempty(varargin)
    method = 'area';
    scale = 1;
else
    if isempty(varargin{1})
        method = 'area';
    else
        if ischar(varargin{1})
            if any(strcmp(varargin{1},{'area';'sum'}))
                method = varargin{1};
            else
                error('1st variable argument is either "area" or "sum"');
            end
        else
            error('1st variable argument is either "area" or "sum"');
        end
    end
    if isempty(varargin{2})
        scale = 1;
    else
        if isscalar(varargin{2})
            scale = varargin{2};
        else
            error('2nd variable argument must be a scalar');
        end
    end
end

if aord == 'a'
    switch method
        case 'area'
            for i = 1:ns 
                H = spectra(:,i*2-1);
                I = spectra(:,i*2);
        %         I(I < 0) = 0;

                % calculate area A
                dH = zeros(size(H));
                dH(2:end) = diff(H);
                A = sum(prod([dH I],2)); 

                spectra(:,i*2) = (I./A).*scale;
            end  
        case 'sum'
            for i = 1:ns 
                H = spectra(:,i*2-1);
                I = spectra(:,i*2);
        %         I(I < 0) = 0;

                % calculate sum A
                A = sum(I); 

                spectra(:,i*2) = (I./A).*scale;
            end  
        otherwise
            error('invalid normalization method');
    end
end
        
if aord == 'd'
    % convert to absorbance spectra
    spectra = deriv2abs(spectra);
    % normalize absorbance spectra
    spectra = normalize_spectra(spectra,'a',varargin{:});
    % convert back to derivative spectra
    spectra = abs2deriv(spectra);
end

return