function Si = interp_spectra(S,C,varargin)
% interp_spectra: interpolates to a certain number of points 
% depending on the interval of interpolation

if all(size(S) > 1) & ndims(S) == 2 % if S is a matrix
    [nb,ncol] = size(S);

    if ncol < 2 || rem(ncol,2) ~= 0
        error('each spectrum is two columns: [B-field intensity_values]');
    end
    
    ns = ncol/2;
else
    error('first argument must be a matrix of the spectra');
end

if isempty(C)
    nc = 0;
    % assume interpolation of magnetic field only
    if isempty(varargin)
        Si = S;
        return
    else
        if ischar(varargin{1}) | isempty(varargin{1}) | isnumeric(varargin{1})
            varargin = {varargin{2:end}};
        end
    
        if isempty(varargin{1}) | isequal(varargin{1},0)
            interpolateB = false;
            nbi = nb;
        elseif isnumeric(varargin{1}) & ~isempty(varargin{1})
            Bi = varargin{1};
            nbi = length(Bi);
        
            if isempty(varargin{2})
                interpolateB = true;
                Bimethod = 'pchip'; % default
            elseif ischar(varargin{2})
                interpolateB = true;
                Bimethod = varargin{2};
            else
                error('if want to interpolate magnetic field values, then please specify the method string or leave empty');
            end
        end
        
        if isempty(varargin{3}) | isequal(varargin{3},0)
            interpolateC = false;
            nci = nc;
        elseif isnumeric(varargin{3}) & ~isempty(varargin{3})
            error('if want to interpolate composition values, then please specify the compositions');
%             Ci = varargin{3};
%             nci = length(Ci);
%             
%             if isempty(varargin{4})
%                 interpolateC = true;
%                 Cimethod = 'linear'; % default
%             elseif ischar(varargin{4})
%                 interpolateC = true;
%                 Cimethod = varargin{4};
%             else
%                 error('if want to interpolate composition values, then please specify the method string or leave empty');
%             end
        end      
    end
elseif  any(size(C) == 1) & any(size(C) > 1) % if C is a vector
    nc = length(C); % nc = number of compositions, should equal ns
    nd = 1;
    [nrow,ncol] = size(C);
    
    if nrow == 1 % if row vector
        C = C'; % convert to column vector 
    end
    
    if nc ~= ns
        error('number of compositions must equal number of spectra');
    end
    
    curve = true;
    lenC = C(end)-C(1);
    
    for c = 1:nc
        pC(c) = (C(c)-C(1))./lenC;
    end
    % pC = diff(C)./lenC;
    
    if isempty(varargin)
        Si = S;
        return
    else
        if ischar(varargin{1}) | isempty(varargin{1}) | isnumeric(varargin{1})
            varargin = {varargin{2:end}};
        end
    
        if isempty(varargin{1}) | isequal(varargin{1},0)
            interpolateB = false;
            nbi = nb;
        elseif isnumeric(varargin{1}) & ~isempty(varargin{1})
            Bi = varargin{1};
            nbi = length(Bi);
        
            if isempty(varargin{2})
                interpolateB = true;
                Bimethod = 'pchip'; % default
            elseif ischar(varargin{2})
                interpolateB = true;
                Bimethod = varargin{2};
            else
                error('if want to interpolate magnetic field values, then please specify the method string or leave empty');
            end
        end
        
        if isempty(varargin{3}) | isequal(varargin{3},0)
            interpolateC = false;
            nci = nc;
        elseif isnumeric(varargin{3}) & ~isempty(varargin{3})
            Ci = varargin{3};
            nci = length(Ci);
            
            if isempty(varargin{4})
                interpolateC = true;
                Cimethod = 'linear'; % default
            elseif ischar(varargin{4})
                interpolateC = true;
                Cimethod = varargin{4};
            else
                error('if want to interpolate composition values, then please specify the method string or leave empty');
            end
        end      
    end
elseif all(size(C) > 1) & ndims(C) == 2 % if C is a matrix 
    [nc,nd] = size(C);
    
    if nc ~= ns
        error('number of compositions must equal number of spectra');
    end
    
    if nd == 2
        tern = false;
    elseif nd == 3
        tern = true;
    else
        error('composition dimension > than which currently can be handled');
    end
    
    if isempty(varargin)
        Si = S;
        return
    else
        switch varargin{1}
            case 'curve'
                curve = true;
                [lenC,pC] = bdy_fxn(C,'linear');
                pC = pC(:,1);
            case 'scattered'
                curve = false;
            otherwise
                error('if compositional dims > 1,1st argument of varargin must be either "curve" or "scattered" specifying how the compositions of the spectra are arranged'); 
        end
    
        if isempty(varargin{2}) | isequal(varargin{2},0)
            interpolateB = false;  
            nbi = nb;
        elseif isnumeric(varargin{2}) & ~isempty(varargin{2})
            Bi = varargin{2};
            nbi = length(Bi);
        
            if isempty(varargin{3})
                interpolateB = true;
                Bimethod = 'pchip'; % default
            elseif ischar(varargin{3})
                interpolateB = true;
                Bimethod = varargin{3};
            else
                error('if want to interpolate magnetic field values, then please specify the method string or leave empty');
            end
        end
        
        if isempty(varargin{4}) | isequal(varargin{4},0)
            interpolateC = false;
            nci = nc;
        elseif isnumeric(varargin{4}) & ~isempty(varargin{4})
            Ci = varargin{4};
            [nci,ndi] = size(Ci);
            
            if size(Ci,2) ~= size(C,2)
                error('the interpolated comps must have same dims as spectral comps');
            end
                
            if isempty(varargin{5})
                interpolateC = true;
                Cimethod = 'linear'; % default
            elseif ischar(varargin{5})
                interpolateC = true;
                Cimethod = varargin{5};
            else
                error('if want to interpolate composition values, then please specify the method string or leave empty');
            end    
        end
    end
end

Si = S;

if interpolateC
    % S = Si;
    % [nb,ncol] = size(S);
    S = align_spectra(S);
    nb = size(S,1);
    Si = zeros([nb nci*2]);

    if curve
        % Si = zeros([nb nci*2]);
        pCi = bdypt2b(C,Ci);
        I = S(:,2:2:end);
        Ii = interp1(pC,I',pCi,Cimethod);
        Si(:,1:2:end) = S(:,1:2:nci*2);
        Si(:,2:2:end) = Ii';
    else
        % Si = zeros([nb nci*2]);
    
        if tern
            C = tern2cart(C,1);
            Ci = tern2cart(Ci,1);
        end
        
        [Ci1,Ci2] = meshgrid(Ci(:,1),Ci(:,2));
        
        warning off
        for b = 1:nb
            Ii(b,:) = diag(griddata(C(:,1),C(:,2),S(b,2:2:end)',Ci1,Ci2,Cimethod))'; % 'nearest' and 'v4' does not seem to return NaNs
        end
        warning on
        
        Si(:,1:2:end) = repmat(S(:,1),1,nci); % S(:,1:2:nci*2);
        Si(:,2:2:end) = Ii;
    end
end

if interpolateB
    % interval = (spectra(end,1) - spectra(1,1))/length(spectra(:,1));
    % interval = round(interval*1000)/1000;
    % Bi = [spectra(1,1):interval:spectra(end,1)]';
    % Si = zeros([nbi nci]);
    S = Si;
    [nb,ncol] = size(Si);
    nsi = ncol/2;
    clear Si;

    for s = 1:nsi
        B = S(:,s*2-1);
    
        if any(diff(B) == 0)
            i = find(diff(B) == 0);
        
            for k = 1:length(i)
                inc = (B(i(k)+2)-B(i(k)-1))/3;
                B(i(k)) = B(i(k)-1)+inc;
                B(i(k)+1) = B(i(k))+inc;
            end
        end
    
        if any(diff(B) < 0)
            i = find(diff(B) < 0);
        
            for k = 1:length(i)
                low = 1;
                high = 2;
                m = (B(i(k)+high)-B(i(k)-low))/(high+low);
            
                while m <= 0
                    low = low-1;
                    high = high+1;
                    m = (B(i(k)+high)-B(i(k)-low))/(high+low);
                end
            
                B(i(k)) = B(i(k)-1)+m;
                B(i(k)+1) = B(i(k))+m;
            end
        end
    
        I = S(:,s*2);
        Ii = interp1(B,I,Bi,Bimethod); 
        Si(:,s*2-1:s*2) = [Bi Ii];
    end

    clear I B Ii;
end

return