function varargout = spectral_prop(spectra,aord,varargin)
% spectra = either one or more spectra with magnetic field values and
% absorbance values
%
%***spectra must be normalized***
%
% aord = 'a' if spectra are absorbance spectra, and 'b' if derivative
% spectra
%
% varargin = prop_type = {'moment', n}, 'skewness', 'kurtosis',
% 'shape_mat', 'lwp_mat', {'characteristic_fxn', order}

[np,nc] = size(spectra);

if nc < 2 || rem(nc,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end

ns = nc/2;
H0 = 3326;

switch varargin{1}
    case 'area'     
        if aord == 'a'
            for i = 1:ns 
                H_field = spectra(:,i*2-1);
                specP = spectra(:,i*2);
                specP(specP < 0) = 0;
    
                % calculate area A
                dH = zeros(size(H_field));
                dH(2:end) = diff(H_field);
                A(i) = sum(prod([dH specP],2)); 
            end
        end
        
        if aord == 'd'
            for i = 1:ns 
                H_field = spectra(:,i*2-1);
                specP = spectra(:,i*2);
    
                % calculate area A
                dH = zeros(size(H_field));
                dH(2:end) = diff(H_field);
                
                sA = 0;
                for j = 1:np
                    sA = sA + prod([dH(j)^2 sum(specP(1:j))],2);
                end
                
                A(i) = sA;
                
                % or 
                % A = sum(prod([abs(dH) dH0 specP],2));
            end
        end
        
        varargout{1} = A';
    case 'moment'
        n = varargin{2};
        
        if aord == 'a'
            for i = 1:ns 
                H_field = spectra(:,i*2-1);
                specP = spectra(:,i*2);
    
                % calculate area A
                dH = zeros(size(H_field));
                dH(2:end) = diff(H_field);
                A = sum(prod([dH specP],2));
            
                % find field H of central peak
                dH0 = H_field-H0;
                H = (sum(prod([dH0 specP],2)))/(sum(specP));
            
                % calculate moment of nth order
                dHH = H_field-H;
                mono(i) = (sum(prod([dH specP dH0.^n],2)))/A; 
            end
        end
        
        if aord == 'd'
            for i = 1:ns 
                H_field = spectra(:,i*2-1);
                specP = spectra(:,i*2);
    
                % calculate area A
                dH = zeros(size(H_field));
                dH(2:end) = diff(H_field);
                
                A = 0;
                for j = 1:np
                    A = A + prod([dH.^2 sum(specP(1:j))],2);
                end
                
                % or 
                % A = sum(prod([abs(dH) dH0 specP],2));
            
                % find field H of central peak
                dH0 = H_field-H0;
                H = 0.5*(sum(prod([dH0.^2 specP],2)))/(sum(prod([dH0 specP],2)));
            
                % calculate moment of nth order
                dHH = H_field-H;
                
                mno = 0;
                for j = 1:np
                    mno = mono + prod([dH.^2 sum(prod([specP(1:j) dH0(1:j).^n],2))],2);
                end
                
                mono(i) = mno/A; 
                
                % or
                % mono = (sum(prod([abs(dH) dHO.^(n+1) specP],2)))/((n+1)*A);
            end
        end
        
        varargout{1} = mono;
    case 'skewness'
        numer = spectral_prop(spectra,aord,'moment',3);
        denom = (spectral_prop(spectra,aord,'moment',2)).^(3/2);
        skew = numer./denom;
        varargout{1} = skew;
    case 'kurtosis'
        numer = spectral_prop(spectra,aord,'moment',4);
        denom = (spectral_prop(spectra,aord,'moment',2)).^2;
        kurt = numer./denom;
        varargout{1} = kurt;
    case 'characteristic_fxn'
        order = varargin{2};
        
        if order == 0
            char_fxn = ones([1 ns]);
        else
            char_fxn = ones([order+1 ns]);
            
            for j = 1:order
                char_fxn(j+1,:) = spectral_prop(spectra,aord,'moment',j);
            end
        end
        
        varargout{1} = char_fxn;
    case 'shape_mat'
        % shape_mat is a matrix of pairwise distances from one point to
        % every other point  sqrt(sum(A.^2,2))
        
        shape_mat = zeros([np np ns]);
        
        for i = 1:ns
            for j = 1:np
                delta_x = spectra(j,i*2-1)-spectra(:,i*2-1);
                delta_y = spectra(j,i*2)-spectra(:,i*2);
                shape_mat(j,:,i) = sqrt(sum([delta_x delta_y].^2,2));
            end
        end
        
        varargout{1} = shape_mat;
    otherwise
        error('no spectral property was given');
end

return        