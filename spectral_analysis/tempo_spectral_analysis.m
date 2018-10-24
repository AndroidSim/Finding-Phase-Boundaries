function varargout = tempo_spectral_analysis(S,refS,varargin)
% this function subtracts the water component from mixture spectra of tempo
% in lipid membrane suspensions and calculates the concentration of tempo
% using the reference spectra.
%
% S = spectra matrix of lipid suspension containing components from the
%        membrane and buffer phase (i.e. spectra of pellet).
% freq_S = frequencies recorded from frequency counter for S.
% pmass = mass of pellet of each S.
% bufS = spectra matrix of tempo in buffer (i.e. spectra of supernatant).
% freq_bufS = frequencies recorded from frequency counter for bufS.
% smass = mass of supernatant of each bufS.
% refS = spectrum of 1 mM tempo in buffer.
% freq_refS = frequency recorded from frequency counter for refS.

[nmp,ncol] = size(S);
if ncol < 2 || rem(ncol,2) ~= 0
    error('each spectrum must contain 2 columns: [B-field absorbance_values]');
end
nms = ncol/2;

% if length(freq_S) ~= nms
%     error('number of frequencies of mixture spectra must = number of mix spectra');
% end

if isempty(refS)
    subtract = false;
else
	[nrp,ncol] = size(refS);
	if ncol < 2 || rem(ncol,2) ~= 0
        error('each spectrum must contain 2 columns: [B-field absorbance_values]');
	end
	nrs = ncol/2;
end

nvarg = length(varargin);
if isempty(varargin)
    normalize = true;
    method = 'spectralp';
else
    if islogical(varargin{1})
        normalize = varargin{1};
        if nvarg > 1
            if any(strcmp(varargin{2},{'spectralp','concentration'}))
                method = varargin{2};
                if strcmp(method,'concentration')
                    comps = varargin{3};
                    suspconc = varargin{4};
                    suspvol = varargin{5};
                    Msolid = varargin{6};
                    Mfluid = varargin{7};
                end
            else
                error('method is invalid')
            end
        else
            method = 'spectralp';
        end
    else
        error('1st var arg must be true or false for normalization');
    end
end
                
% shift middle peak to median field
Bfields = S(:,1:2:end);
[maxmid,imm] = max(S(400:600,2:2:end));
medianB = median(S(:,1));
for s = 1:nms
    S(:,s*2-1:s*2) = shift_spectra(S(:,s*2-1:s*2),medianB-Bfields(399+imm(s),s));
end
if subtract
    [maxmidr,immr] = max(refS(400:600,2:2:end));
    refS = shift_spectra(refS,medianB-refS(immr,1));
end

% align spectra
if subtract
    temp = [S refS];
    temp = align_spectra(temp);
    S = temp(:,1:nms*2);
    refS = temp(:,end-1:end);
else
    S = align_spectra(S);
end

% smooth spectra (and/or interpolate?)
for s = 1:nms
    S(:,s*2) = smooth(S(:,s*2-1),S(:,s*2),3,'moving');
%     S(:,s*2) = smooth(S(:,s*2-1),S(:,s*2),20,'sgolay',2);
end

if subtract
    refS(:,2) = smooth(refS(:,1),refS(:,2),3,'moving');
%     refS(:,2) = smooth(refS(:,1),refS(:,2),20,'sgolay',2);
end
% interpolate
% [np,ncol] = size(S);
% for s = 1:nms
%     temp(:,s*2-1) = linspace(S(1,s*2-1),S(end,s*2-1),np*2)';
%     temp(:,s*2) = interp1(S(:,s*2-1),S(:,s*2),linspace(S(1,s*2-1),S(end,s*2-1),np*2)','pchip');% spline interpolation
%     temp(isnan(temp(:,s*2)),s*2) = 0;
% end
% S = temp;
% clear temp;
% if subtract
%     refS(:,2) = interp1(refS(:,1),refS(:,2),'spline');
%     refS(isnan(refS(:,2)),2) = 0;
% end

if subtract
    % subtract buffer spectra from mixture spectra
    for m = 1:nms
        [sS(:,m*2-1:m*2),bS(:,m*2-1:m*2),p(:,m)] = spectral_subtract(S(:,m*2-1:m*2),refS);
    end
    
    switch method
        case 'spectralp'
            M = max(sS(600:900,2:2:end));
            W = max(bS(600:900,2:2:end));
            varargout{1} = sS;
            varargout{2} = bS;
            varargout{3} = refS;
            varargout{4} = M;
            varargout{5} = W;
            varargout{6} = M./(M+W);
        case 'concentration' 
            % calculate area of absorption spectrum of subtracted spectrum and background
			% spectrum
			sS = deriv2abs(sS);
			bS = deriv2abs(bS);
			rS = deriv2abs(refS);
			% AsS = spectral_prop(sS,'a','area');
			% AbS = spectral_prop(bS,'a','area');
			% ArS = spectral_prop(rS,'a','area');
			for m = 1:nms
                AsS(m,1) = max(sS(200:400,m*2)) + max(sS(400:600,m*2)) + max(sS(600:800,m*2));
                AbS(m,1) = max(bS(200:400,m*2)) + max(bS(400:600,m*2)) + max(bS(600:800,m*2));
			end
			ArS = max(rS(200:400,2)) + max(rS(400:600,2)) + max(rS(600:800,2));
			sS = abs2deriv(sS);
			bS = abs2deriv(bS);
			rS = abs2deriv(rS);
            
            varargout{1} = sS;
            varargout{2} = bS;
            varargout{3} = rS;
            varargout{4} = AsS;
            varargout{5} = AbS;
            varargout{6} = ArS;
        otherwise
            error('invalid method');
	end
else 
    switch method
        case 'spectralp'
            % convert spectra to absorbance, normalize, convert back to derivative
            if normalize
                S = deriv2abs(S);
                S = normalize_spectra(S,'a');
                S = abs2deriv(S);
            end
            % find membrane peak and water(buffer) valley
            [B,M,indices] = peak(S,3332); %3332 3490 3330
            [B,W,indices] = valley(S,3336);% 3335 3493 3339
		%     M = max(S(660:720,2:2:end));
		%     W = min(S(720:780,2:2:end));
            varargout{1} = M';
            varargout{2} = abs(W');
            varargout{3} = M'./(M'+abs(W'));
            varargout{4} = S;
        case 'concentration' 
            % calculate area of absorption spectrum 
			absS = deriv2abs(S);
			A = spectral_area(S);
            mass = suspconc.*suspvol;
			amt = comps2amt(comps,mass,Msolid,Mfluid);
            Xpm = A./amt;
            
            varargout{1} = A;
            varargout{2} = absS;
            varargout{3} = mass;
            varargout{4} = Xpm;
        otherwise
            error('invalid method');
	end
end

return