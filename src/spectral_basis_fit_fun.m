function V = spectral_basis_fit_fun(x,c,method,varargin)

switch method
    case 'grid'
        % varargin = alpha,beta,B;
        alpha = varargin{1};
        beta = varargin{2};
        B = varargin{3};
        
        V = ((x.*(beta-c))./(c-alpha+(x.*(beta-c)))).*B(:,1) + ((c-alpha)./(c-alpha+(x.*(beta-c)))).*B(:,2);
    case 'cc'
        S = varargin{1};
        C = varargin{2};
        % interpolate spectra matrix S to find the basis spectra at x(2) = alpha
        % and x(3) = beta
        [nb,nc] = size(S);
        Ci = [x(2) x(3)];

        for b = 1:nb
            %C = c;
            A = S(b,:);
            Ai = interp1(C,A,Ci); % default = linear
            B(b,:) = Ai;
        end
        
        %for s = 1:nc
        %    D(:,s) = ((x(1).*(x(3)-c(s)))./(c(s)-x(2)+(x(1).*(x(3)-c(s))))).*B(:,1) + ((c(s)-x(2))./(c(s)-x(2)+(x(1).*(x(3)-c(s))))).*B(:,2);
        %end
        
        %F = reshape(D,nb*nc,1); 
        
        V = ((x(1).*(x(3)-c))./(c-x(2)+(x(1).*(x(3)-c)))).*B(:,1) + ((c-x(2))./(c-x(2)+(x(1).*(x(3)-c)))).*B(:,2);
    otherwise
        error('invalid method in fit function');
end