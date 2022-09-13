function [fit_out] = spectral_basis_fit(S,C,method)
% S = spectra
% C = compositions of the spectra
% method = basis search method, either 'grid' or 'cc'
% 'cc' = 'continuous composition'

warning off
nc = length(C);
[nb,ns] = size(S);

if nc ~= ns
    error('number of compositions must equal number of spectra');
end

% B = basis set 
% S(:,i) = spectrum
  
switch method
    case 'grid'
        options = optimset('display','off');
        kp = repmat(NaN,[nc nc]);
        stdev_kp = repmat(NaN,[nc nc]);
        chisq = repmat(NaN,[nc nc]);
        
        % grid search for basis
        for i = 2:nc-3
            alpha = C(i);
    
            for k = i+2:nc-1
                beta = C(k);
        
                % choose basis
                B = [S(:,i) S(:,k)];
        
                % only fitting parameter is kp
                for s = 1:ns
                    [x(s),resnorm(s)] = lsqcurvefit(@spectral_basis_fit_fun,1,C(s),S(:,s),0.1,10,options,method,alpha,beta,B); 
                end
        
                kp(i,k) = mean(x);
                stdev_kp(i,k) = std(x);
                chisq(i,k) = (1/(k-i+1))*sum(resnorm); %*stdev_kp(i,k); 
            end
        end
        
        fit_out = struct('kp',kp,'stdevkp',stdev_kp,'chisq',chisq);
    case 'cc'
        % continuous composition search
        % fitting parameters = x = [kp alpha beta]
        options = optimset('display','iter');
        %D = reshape(S,nb*ns,1); 
        %[x,resnorm] = lsqcurvefit(@spectral_basis_fit_fun,[1 0.0001 0.9999],C,D,[0.1 0 0],[10 1 1],options,method,S);  
        %fit_out = struct('kp',x(1),'calpha',x(2),'cbeta',x(3),'chisq',resnorm);
        
        for s = 1:ns
            [x(:,s),resnorm(s)] = lsqcurvefit(@spectral_basis_fit_fun,[1 0.0001 0.9999]',C(s),S(:,s),[0.1 0 0]',[10 1 1]',options,method,S,C);
        end
        
        fit_out = struct('kp',x(1,:),'calpha',x(2,:),'cbeta',x(3,:),'chisq',resnorm);
    otherwise
        error('invalid method');
end

warning on
return