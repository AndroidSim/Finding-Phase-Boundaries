function fit_output = phase_boundary_fit(spectra,C,p,basis_method,varargin)
% spectra should be derivative spectra from .dat files converted from ascii
% (.asc) files
%
% spectra = derivative or absorbance spectra in two columns: 
%   [magnetic field values  intensity values]
%
% C = vector or matrix of compositions of each spectrum
%
% p = the number of phases of the region to find the boundaries for
%
% basis_method = the method used to form the basis for the spectra
%   either "svd" for singular value decomposition or
%          "spectral" for the spectra themselves

% checking the input arguments
disp(sprintf('checking input arguments and spectrum pre-processing:\n'));

if all(size(spectra) > 1) & ndims(spectra) == 2 % if spectra is a matrix
    [nb,ncol] = size(spectra);

    if ncol < 2 || rem(ncol,2) ~= 0
        error('each spectrum is two columns: [B-field intensity_values]');
    end
    
    ns = ncol/2;
else
    error('first argument must be a matrix of the spectra');
end

if  any(size(C) == 1) & any(size(C) > 1) % if C is a vector
    nc = length(C); % nc = number of compositions, should equal ns
    nd = 1;
    [nrow,ncol] = size(C);
    
    if nrow == 1
        C = C';
    end
    
    if nc ~= ns
        error('number of compositions must equal number of spectra');
    end
elseif all(size(C) > 1) & ndims(C) == 2 % if C is a matrix 
    [nc,nd] = size(C);
    
    if nd == 2
        tern = false;
        cart = true;
    elseif nd == 3
        tern = true;
        cart = false;
    end
else % if C is not a vector, ie all(size(C) == 1) | all(size(C) > 1) | ndims(C) > 2
    error('second argument must be a vector specifying compositions of spectra');
end

if all(size(p) == 1) & ndims(p) == 2 % if p is a scalar
    disp(sprintf('check: finding the boundaries of a %d-phase region.\n',p));
else
    error('third argument must be a scalar specifying the number of phases');
end

% check to see if spectra are derivative or absorbance spectra
for s = 1:ns
    approx_bline = [spectra(1:50,s*2);spectra(end-50:end,s*2)]; % approx_bline = approximate baseline
    llim_noise = min(approx_bline); % llim_noise = lower limit of noise and noise = mean absolute value of baseline
    
    if any(spectra(400:600,s*2) < llim_noise)
        absorbance(s) = false;
    else
        absorbance(s) = true;
    end
end
if all(absorbance == true)
    aord = 'a';
elseif all(absorbance == false)
    aord = 'd';
else
    indt = find(absorbance == true);
    indf = find(absorbance == false);
    disp(sprintf('The %d-th spectrum is an absorbance spectrum\n',indt));
    disp(sprintf('The %d-th spectrum is a derivative spectrum\n',indf));
    error('all spectra must be of the same type');
end

% perform the fitting procedure to determine the phase boundaries
switch basis_method
    case 'svd'
        % the following code is modified from a previous m-file called
        % 'boundary_spectra_lsqlin.m'
        
        disp(sprintf('beginning the svd phase boundary fitting procedure\n'));
        disp(sprintf('performing spectra pre-processing:\n'));
        
        % align all spectra
        spectra = align_spectra(spectra);
        if aord == 'd'
            % convert derivative spectra to absorbance spectra
            spectra = deriv2abs(spectra);
            aord = 'a';
            spectra = normalize_spectra(spectra,aord);
            if any(round(spectral_area(spectra).*100)./100 ~= 1)
                % normalize absorbance spectra
                spectra = normalize_spectra(spectra,aord);
            end
        else
            if any(round(spectral_area(spectra).*100)./100 ~= 1)
                % normalize absorbance spectra
                spectra = normalize_spectra(spectra,aord);
            end
        end
        
        % get baseline statistics
        for s = 1:ns
            Sbaseline(:,s) = spectra([1:200]',s*2);
            std_Sbline(s) = std(Sbaseline(:,s));
        end 
%         stdS = std(S')';
        
        % set baseline (negative absorbances) to zero
        % and   
        % renormalize spectra so that the sum of absorbances equals one
        for s = 1:ns
            spectra(find(spectra(:,s*2) <= std_Sbline(s)),s*2) = std_Sbline(s);
%             spectra(find(spectra(:,s*2) <= std_Sbline(s)),s*2) = 10^-4;
%             spectra(find(spectra(:,s*2) < 0),s*2) = 0;
%             ind = find(spectra(:,s*2) < 0);
%             spectra(ind,s*2) = abs(spectra(ind,s*2));
            spectra(:,s*2) = spectra(:,s*2)./sum(spectra(:,s*2));
        end
        
        % reduce spectra to just absorbance values 
        S = spectra(:,2:2:end);
        
%         weight_fxn = ones(size(stdS));
%         weight_fxn = stdbS./sum(stdbS);
%         w = zeros(size(weight_fxn));
%         w(weight_fxn >= 10*mean(std_bSbline)) = 1;
%         w(weight_fxn >= 100*mean(std_bSbline)) = 2;
%         w(weight_fxn >= 1000*mean(std_bSbline)) = 3;
%         w(weight_fxn >= 10000*mean(std_bSbline)) = 4;
        
        % remove baseline or unimportant magnetic fields
%         S = S([200:400 550:850]',:);
%         I = S(200:end,:);
%         S = S(200:end-200,:);
%         S = S(200:1000,:);
%         S = S([stdS >= mean(std_Sbline)],:);
        
        % subtract mean spectrum
%         S = S-repmat(mean(S,2),1,ns);

        % perform singular value decomposition to obtain eigenspectra basis
        disp(sprintf('performing singular value decomposition to obtain eigenspectra basis:\n'));
        [U,W,V] = svd(S,0);
        B = U(:,1:p); % B = p component eigenbasis

        % solve for eigenspectra coefficient matrix D by linear least squares
        disp(sprintf('solving for eigenspectra coefficient matrix by linear least squares:\n'));
        options = optimset('display','off','largescale','off');
        warning off;

        % constraints
%         Aeq = [sum(U(:,1)) sum(U(:,2))];
%         Aeq = sum(U(:,1:p));
        Aeq = sum(B);
        beq = sum(S);
%         beq = ones(1,ns);
        A = -B;
        b = zeros(size(U(:,1)));
        D = W(1:p,1:p)*V(:,1:p)';

        x0 = D;
        for s = 1:ns
            disp(sprintf('\t fitting spectrum %d, %f percent done',s,((s-1)/ns)*100));
%             [x,resnorm,residual] = lsqlin(B,S(:,s),A,b,Aeq,beq,[],[],x0(:,s),options); % ...[],[],x0,options
            [x,resnorm,residual] = lsqlin(B,S(:,s),A,b,Aeq,beq(s),[],[],x0(:,s),options);
%             [x1,resnorm,residual] = lsqlin(B,S(:,s),A,b,Aeq,beq(s),[],[],x0(:,s),options);
%             [x2,resnorm,residual] = lsqlin(B,S(:,s),A,b,[],[],[],[],x0(:,s),options);
%             [x3,resnorm,residual] = lsqlin(B,S(:,s),[],[],Aeq,beq(s),[],[],x0(:,s),options);
%             x = lsqlin(B,S(:,s),[],[],Aeq,beq);
%             x = lsqlin(B,S(:,s),A,b);
            % x = B\S(:,s);
            D_con(:,s) = x;
%             D_con1(:,s) = x1;
%             D_con2(:,s) = x2;
%             D_con3(:,s) = x3;
        end
 
        Sbar = B*D_con;
        Sstar = B*D;
        D_con = D_con';
        D = D';
        matrix = (B'*B);
        covar = inv(matrix);
                    
        % solve for boundary coordinates a_alpha,b_alpha,a_beta,and b_beta using D
        % D = [a b]; x = [kp a_alpha a_beta c_alpha c_beta];

        disp(sprintf('solving for boundary coordinates using coefficient matrix:\n'));
        options = optimset('display','iter','largescale','on','Jacobian','off','tolfun',10^-10);
        maxa = max(D_con(:,1));
        mina = min(D_con(:,1));
%         x0 = [1 D_con(1,1) D_con(end,1) C(1)+0.0001 C(end)-0.0001];
        x0 = [1 mina maxa C(1) C(end)];
        lb = [0.1 mina mina C(1) C(1)];
        ub = [10 maxa maxa C(end) C(end)];
%         x0 = [1 C(1) C(end)];
%         lb = [0.1 C(1) C(1)];
%         ub = [10 C(end) C(end)];
        [x,resnorm] = lsqcurvefit(@phase_boundary_fit_fun,x0,C,D_con(:,1),lb,ub,options,basis_method,D_con(:,1));
%         options = optimset('display','iter','largescale','off','tolfun',10^-6);
%         [x,resnorm] = fmincon(@phase_boundary_fit_fun,x0,[],[],[],[],lb,ub,[],options,basis_method,D_con(:,1),C);
        T = feval(@phase_boundary_fit_fun,x,C,basis_method,D_con(:,1));

        disp(sprintf('done!\n'));
        fit_output = struct('Kp',x(1),'C_alpha',x(4),'C_beta',x(5),'resnorm',resnorm,'eig_coeff_uncon',D,'eig_coeff_con',D_con,'theory_coeff',T,...
            'U',U,'V',V,'W',W,'B',B,'Sbar',Sbar,'Sstar',Sstar,'S',S,'spectra',spectra,'stdbline',std_Sbline);
%         fit_output = struct('Kp',x(1),'C_alpha',x(2),'C_beta',x(3),'resnorm',resnorm,'eig_coeff_uncon',D,'eig_coeff_con',D_con,'theory_coeff',T,...
%             'U',U,'V',V,'W',W,'B',B,'Sbar',Sbar,'S',S);
        
        disp(sprintf('plotting constrained coefficients of 1st eigenspectra (experiment) and fitted curve (theory) as a function of composition\n'));
        figure;
        plot(C,[D(:,1) D_con(:,1) T]);
        warning on;
    case 'spectral'
        if isempty(varargin)
            search_method = 'cc';
        else  
            if ischar(varargin{1})
                if any(strcmp(varargin{1},{'grid';'cc'}))
                    search_method = varargin{1};
                else
                    disp('using default cc search method for spectral basis method');
                    search_method = 'cc';
                end
            else
                error('fifth argument must be a string');
            end
        end
        
        disp(sprintf('beginning the spectral phase boundary fitting procedure\n'));
        disp(sprintf('performing spectra pre-processing:\n'));
        
        % cut off baseline and match length of cut spectra
        [spectra,baseline] = cut_n_match(spectra,aord);
        
        % get baseline statistics
        avg_bline = repmat(NaN,ns,1);
        std_bline = repmat(NaN,ns,1);
        highest_pos = repmat(NaN,ns,1);
        lowest_neg = repmat(NaN,ns,1);

        for s = 1:ns
            if isempty(baseline{s}), continue, end
            avg_bline(s) = mean(baseline{s}(:,2));
            std_bline(s) = std(baseline{s}(:,2),1);
    
            temp = max(baseline{s}(baseline{s}(:,2) > 0,2));
            if ~isempty(temp)
                highest_pos(s) = temp;
            end
    
            temp = min(baseline{s}(baseline{s}(:,2) < 0,2));
            if ~isempty(temp)
                lowest_neg(s) = temp;
            end
        end

        % interpolate spectra for constant magnetic field increments
        interval = (spectra(end,1) - spectra(1,1))/length(spectra(:,1));
        interval = round(interval*1000)/1000;
        Bi = [spectra(1,1):interval:spectra(end,1)]';
        ispectra = zeros([length(Bi) nc]);

        for s = 1:ns
            B = spectra(:,s*2-1);
    
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
    
            A = spectra(:,s*2);
            Ai = interp1(B,A,Bi); % default is linear
            ispectra(:,s*2-1:s*2) = [Bi Ai];
        end

        clear A B Ai Bi;
        spectra = ispectra;

        % normalize spectra so that the sum of absorbances equals one
        %for s = 1:ns
        %    spectra(:,s*2) = spectra(:,s*2)./sum(spectra(:,s*2));
        %end
        
        % reduce spectra to just intensity values 
        S = spectra(:,2:2:end);
        
        switch search_method
            case 'grid'
                disp(sprintf('beginning grid search method:\n'));
                % from kp_bd_fit.m
                options = optimset('display','off');
                warning off;
                Kp = cell([ns ns]);
                stdev_Kp = cell([ns ns]);
                chisq = cell([ns ns]);
                coeff_theory = cell([ns ns]);

                % start grid search
                for sa = 2:ns-2
                    alpha = C(sa);
    
                    for sb = sa+1:ns-1
                        beta = C(sb);
        
                        % choose spectral basis
                        B = [S(:,sa) S(:,sb)];
        
                        % linear least squares with constraints on
                        % parameters
        
                        disp(sprintf('solving for basis coefficients and fitting spectra for basis %d - %d:\n',sa,sb));
                        % constraints
                        Aeq = [1 1];
                        beq = 1;
                        A = [-1 0;0 -1];
                        b = [0;0];
        
                        for s = 1:ns
                            disp(sprintf('\t spectrum %d, %f percent done\n',s,((s-1)/ns)*100));
                            coeff(:,s) = lsqlin(B,S(:,s),A,b,Aeq,beq,[],[],[],options);% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options)
                            [x(s),resnorm(s)] = lsqcurvefit(@phase_boundary_fit_fun,1,C(s),S(:,s),0.1,10,options,basis_method,search_method,'spectra',alpha,beta,B); 
                        end
        
                        disp(sprintf('fitting basis coefficients:\n'))
                        coeff = coeff';
                        [kp,resnormp] = lsqcurvefit(@phase_boundary_fit_fun,1,C,coeff(:,1),0.1,10,options,basis_method,search_method,'coeff',alpha,beta); 
                        
                        Kp{sa,sb} = [kp mean(x)];
                        stdev_Kp{sa,sb} = [0 std(x)];
                        chisq{sa,sb} = [resnormp (1/(sb-sa+1))*sum(resnorm)] ; %*stdev_kp(i,k);
                        coeff_theory{sa,sb} = feval(@phase_boundary_fit_fun,kp,C,basis_method,search_method,'coeff',alpha,beta);
                        clear coeff;              
                    end
                end
        
                disp(sprintf('done!\n'));
                fit_output = struct('Kp',{Kp},'stdev_kp',{stdev_Kp},'chisq',{chisq},'theory_coeff',{coeff_theory});
                warning on;
            case 'cc'
                disp(sprintf('beginning continuous composition search method:\n'));
                % continuous composition search
                % fitting parameters = x = [kp alpha beta]
                lsq_nlcon = varargin{2};
                options = optimset('display','iter');
                warning off;
                Kp = cell([1 2]);
                C_alpha = cell([1 2]);
                C_beta = cell([1 2]);
                chisq = cell([1 2]);
                
                disp(sprintf('fitting all spectra, parameters = C_phase_boundaries:\n'));
                % fit the coefficients of the basis spectra as a function
                % of the parameters C_alpha, and C_beta
                [numr,numc] = size(S);
                D_data = reshape(S,numr*numc,1);
                
                try
                    x0 = [C(1)+0.0001 C(end)-0.0001];
                    params = lsqcurvefit(@phase_boundary_fit_fun,x0,C,D_data,[C(1) C(1)],[C(end) C(end)],options,basis_method,search_method,'coeff',S);
                catch
                    % disp(sprintf('fitting procedure crashed with following parameters: %f\n',params));
                    rethrow(lasterror);
                end
                
                disp(sprintf('recomputing theoretical spectra based of basis spectra at phase boundaries:\n'));
                D_theory = feval(@phase_boundary_fit_fun,params,C,basis_method,search_method,'coeff',S);
                D_theory = reshape(D_theory,numr,numc);
                Ci = [params(1) params(2)];

                for b = 1:numr
                    A = D_theory(b,:);
                    Ai = interp1(C,A,Ci); % default = linear
                    B(b,:) = Ai;
                end
                
                disp(sprintf('solving for basis coefficients by linear lsq using theoretic spectra and basis spectra:\n'));
                ops = optimset('display','off');
                % constraints
                Aeq = [1 1];
                beq = 1;
                A = [-1 0;0 -1];
                b = [0;0];
                
                for s = 1:ns
                    coeffs(:,s) = lsqlin(B,D_theory(:,s),A,b,Aeq,beq,[],[],[],ops);% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options)
                end
                
                coeffs = coeffs';
                
                disp(sprintf('fitting basis coefficients for Kp:\n'));
                
                try
                    [kp,resnormparams] = lsqcurvefit(@phase_boundary_fit_fun,1,C,coeffs(:,1),0.1,10,options,basis_method,'grid','coeff',params(1),params(2));
                catch
                    % disp(sprintf('fitting procedure crashed with following parameters: %f\n',kp));
                    rethrow(lasterror);
                end
                    
                % fit the actual data spectra as a function of the
                % parameters Kp, C_alpha, and C_beta and the basis spectra
                disp(sprintf('fitting actual data spectra, parameters = [Kp C_phase_boundaries]:\n'));

                try
                    switch lsq_nlcon
                        case 'no'
                            options = optimset('display','iter','GradObj','on','largescale','on');
                            x0 = [1 C(1)+0.0001 C(end)-0.0001]';
                            [x,resnorm,residual,exitflag] = lsqcurvefit(@phase_boundary_fit_fun,x0,C,D_data,[0.1 C(1) C(1)]',[10 C(end) C(end)]',options,basis_method,search_method,'spectra',S);
                        case 'yes'
                            options = optimset('display','iter','GradObj','on','GradConstr','on','largescale','on','MaxFunEvals',500);
                            x0 = [1 C(1)+0.0001 C(end)-0.0001]';
                            % constraints
                            A = [0 1 -1];
                            b = 0;
                            [x,resnorm,exitflag] = fmincon(@phase_boundary_con_fit_fun,x0,A,b,[],[],[0.1 C(1) C(1)]',[10 C(end) C(end)]',[],options,S,C,D_data);%@phase_boundary_nlcon
                    end
                catch
                    disp(sprintf('lsqcurvefit for spectra fitting crashed at spectrum %d (composition = %f)\n',s,C(s)));
                    % disp(sprintf('and with the following parameters: %f  %f  %f\n',x));
                    rethrow(lasterror);
                end
        
                Kp = {kp,x(1)};
                C_alpha = {params(1),x(2)};
                C_beta = {params(2),x(3)};
                chisq = {resnormparams,resnorm};
                
                disp(sprintf('done!\n'));
                fit_output = struct('Kp',{Kp},'C_alpha',{C_alpha},'C_beta',{C_beta},'chisq',{chisq});
                warning on;
            otherwise
                error('invalid search method');
        end
    otherwise
        error('fourth argument must be string designating the basis method');
end

return