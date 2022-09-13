function fit_output = phase_boundary_fitting(C,data,phase_rule,data_specs,varargin)
% input arguments:
% data = FRET vector or derivative/absorbance spectra matrixd divided into
%        a number of two columns: [magnetic field values  intensity values]
% C = vector or matrix of compositions of each spectrum
% phase_rule = [number_components number_phases]
% data_specs = structure with fields:
%   type = "FRET" or "spectra"
%   basis = "svd" or "spectra" (only with data_specs.type = "spectra")
% varargin = structure of fitting conditions with some fields:
%   fit_method = "lsqnonlin","lsqcurvefit","simplex","fmincon"
%   parameter position ("range","located") and values
%   x0 = initial parameter vector
%
% output arguments:
% fit_output = structure of the results of the fitting

% checking the input arguments
disp(sprintf('checking input arguments and performing data pre-processing:\n'));

if ~isvector(data) && ~isscalar(data) % all(size(spectra) > 1) && ndims(spectra) == 2 
    % if data is a matrix
    [nrow,ncol] = size(data);
    switch data_type
        case 'FRET'
            ndata = nrow;
            ntraj = ncol;
        case 'spectra'
            nBfield = nrow;
            if rem(ncol,2) == 0
                nS = ncol/2;
            else
                error('each spectrum is two columns: [B-field intensity_values]');
            end 
        otherwise
            error('invalid data type');
    end
elseif isvector(data)
    % if data is a vector
    [nrow,ncol] = size(data);
    if nrow == 1
        data = data';
        nrow = ncol;
    end
    switch data_type
        case 'FRET'
            ndata = nrow;
            ntraj = 1;
        case 'spectra'
            error('each spectrum is two columns: [B-field intensity_values]');
%           nBfield = nrow;
%           nS = 1; 
        otherwise
            error('invalid data type');
    end
else
    error('first argument must be a vector or matrix of data');
end

if ~isvector(C) && ~isscalar(C) % all(size(C) > 1) && ndims(C) == 2 
    % if C is a matrix
    [nrow,ncol] = size(C);
    nC = nrow;
    nCdim = ncol;  
    Ctraj = linspace(0,1,nC)';
elseif isvector(C) % any(size(C) == 1) && any(size(C) > 1) 
    % if C is a vector
    [nrow,ncol] = size(C);
    if nrow == 1
        C = C';
        nrow = ncol;
    end 
    nC = nrow;
    nCdim = 1;
    Ctraj = linspace(0,1,nC)';
else % if C is a scalar, ie all(size(C) == 1) 
    error('second argument must be a matrix or vector of compositions');
end

%declare and initialize phase system structure
phase_system = struct('comps',C,'nCdim',nCdim,'Ctraj',Ctraj,'data',data);

if isscalar(phase_rule) % && isinteger(phase_rule)  
    % if phase_rule is a scalar assumed to specify number of phases
%     disp(sprintf('check: finding the boundaries of a %d-phase region.\n',p));
    phase_system.p = phase_rule;
    phase_system.c = nCdim;
    phase_system.f = phase_system.c-phase_system.p;
elseif isvector(phase_rule) 
    phase_system.c = phase_rule(1);
    phase_system.p = phase_rule(2);
%     if isinteger(phase_rule(1))
%         phase_system.c = phase_rule(1);
%     else
%         error('number of chemical components must be an integer');
%     end
%     if isinteger(phase_rule(2))
%         phase_system.p = phase_rule(2);
%     else
%         error('number of phases must be an integer');
%     end
    phase_system.f = phase_system.c-phase_system.p;
elseif isstruct(phase_rule)
    if isfield(phase_rule,'ncomps')
        phase_system.c = phase_rule.ncomps;
    else
        error('phase rule structure must have a field "ncomps" containing the number of chemical components');
    end
    if isfield(phase_rule,'nphases')
        phase_system.p = phase_rule.nphases;
    else
        error('phase rule structure must have a field "nphases" containing the number of phases');
    end
    phase_system.f = phase_system.c-phase_system.p;
else  
    error('phase_rule must be a scalar = number of phases, vector = [#comps #phases], or a structure with fields "ncomps" and "nphases"');
end

if phase_system.p < 2
    error('number of phase must be >= 2');
elseif phase_system.c < 1
    error('number of chemical components must be >= 1');
elseif phase_system.c < nCdim
    error('number of chemical components must be >= nCdim');
else
    switch phase_system.p
        case 2
            switch phase_system.c
                case 1
                    phase_system.name = 'binary';
                    phase_system.Ccoord = '1X';
                case 2
                    switch nCdim
                        case 1
                            phase_system.name = 'binary';
                            phase_system.Ccoord = '1X';
                        otherwise % nCdim == 2
                            if all(round(sum(C,2)) == 1)
                                phase_system.name = 'binary';
                                phase_system.Ccoord = '2X';
                            else 
                                error('binary compositions must sum to 1');
                            end
                    end  
                case 3
                    switch nCdim
                        case 1
                            phase_system.name = 'ternary';
                            phase_system.Ccoord = '1X';
                        case 2
                            if all(round(sum(C,2)) == 1)
                                error('C for binary system ~= 3 components');
                            else 
                                phase_system.name = 'ternary2p';
                                phase_system.Ccoord = '2Xcart';
                                C = cart2tern(C); 
                                error('cannot determine for 2-phase regions in ternary composition space, for now');
                            end
                        otherwise % nCdim == 3
                            if all(round(sum(C,2)) == 1)
                                phase_system.name = 'ternary2p';
                                phase_system.Ccoord = '3X';
                            else  
                                error('ternary compositions must sum to 1');
                            end
                    end
                otherwise
                    % ncomp > 3, more than ternary
                    error('number of chemical components must be <= 3, for now');
            end
        case 3
            switch ncomps
                case 1
                    % ternary trajectory through 3-phase region
                    ternary3 = true;
                    error('cannot determine for 3-phase regions in ternary composition space, for now');
                case 2
                    if all(round(sum(C,2)) == 1)
                        % binary
                        error('degrees of freedom(f) = -1, a violation of the phase rule');
                    else % ternary cartesian
                        ternary3 = true;
                        C = cart2tern(C);
                        error('cannot determine for 3-phase regions in ternary composition space, for now');
                    end
                case 3
                    % ternary
                    ternary3 = true;
                    error('cannot determine for 3-phase regions in ternary composition space, for now');
                otherwise
                    % ncomps > 3, more than ternary
                    error('number of chemical components must be <= 3, for now');
             end
        otherwise
            error('cannot determine boundaries when number of phases > 3');
    end
end

if isstruct(data_specs)
    if isfield(data_specs,'type')
        if strcmp(data_specs.type,'FRET')
            phase_system.data_type = data_specs.type;
        elseif strcmp(data_specs.type,'spectra')
            phase_system.data_type = data_specs.type;
            if isfield(data_specs,'basis')
                if any(strcmp(data_specs.basis,{'svd','spectra'}))
                    phase_system.data_basis = data_specs.basis;
                else
                    error('data basis must be "svd" or "spectra"');
                end
            else
                error('data_specs must contain a field "basis" if data type = "spectra"');
            end
        else
            error('data type must be "FRET" or "spectra"');
        end
    else
        error('data_specs must contain a field "type"');
    end
    
else
    error('data_specs must be a structure');
end

if isempty(varargin)
    interactive = true;
    reply = input(sprintf('"lsqnonlin", "fmincon", "lsqcurvefit", or "simplex" fit method?\n'),'s');
    while isempty(reply) || ~any(strcmp(reply,{'lsqnonlin','fmincon','lsqcurvefit','simplex'}))
        reply = input(sprintf('"lsqnonlin", "fmincon", "lsqcurvefit", or "simplex" fit method?\n'),'s');
    end
    fit_conditions.fit_method = reply;  
elseif isstruct(varargin{1})
    interactive = false;
    fit_structure = varargin{1};
    if isfield(fit_structure,'fit_method')
        if isempty(fit_structure.fit_method)
            % default fit method
            fit_conditions.fit_method = 'lsqnonlin';
        elseif any(strcmp(fit_structure.fit_method,{'lsqnonlin','fmincon','lsqcurvefit','simplex'}))
            fit_conditions.fit_method = fit_structure.fit_method; 
        else
            error('fit method must be either "lsqnonlin","fmincon","lsqcurvefit",or "simplex"');
        end   
    else
        fit_conditions.fit_method = 'lsqnonlin';
%         error('fit_conditions structure must have a field "fit_method" containing the fitting method');
    end
    if isfield(fit_structure,'ta')
        if isscalar(fit_structure.ta)
            fit_conditions.taposition = 'located';
            fit_conditions.taparameter = fit_structure.ta;
        elseif isvector(fit_structure.ta)
            fit_conditions.taposition = 'range';
            temp = fit_structure.ta;
            if temp(1) > temp(end)
                fit_conditions.taparameter = [temp(end) temp(1)];
            elseif temp(1) < temp(end)
                fit_conditions.taparameter = [temp(1) temp(end)];
            else
                error('ua range is same point');
            end
        elseif isempty(fit_structure.ta)
            fit_conditions.taposition = 'range'; % or 'unknown' ?
            fit_conditions.taparameter = [0 1];
        else
            error('field "ta" must be a scalar, a vector specifying a range, or empty');
        end
    else
        error('fit structure must have a field "ta" containing the parameter of the alpha phase boundary');
    end
    if isfield(fit_structure,'tb')
        if isscalar(fit_structure.tb)
            fit_conditions.tbposition = 'located';
            fit_conditions.tbparameter = fit_structure.tb;
        elseif isvector(fit_structure.tb)
            fit_conditions.tbposition = 'range';
            temp = fit_structure.tb;
            if temp(1) > temp(end)
                fit_conditions.tbparameter = [temp(end) temp(1)];
            elseif temp(1) < temp(end)
                fit_conditions.tbparameter = [temp(1) temp(end)];
            else
                error('tb range is same point');
            end
        elseif isempty(fit_structure.tb)
            fit_conditions.tbposition = 'range'; % or 'unknown' ?
            fit_conditions.tbparameter = [0 1];
        else
            error('field "tb" must be a scalar, a vector specifying a range, or empty');
        end
    else
        error('fit structure must have a field "tb" containing the parameter of the beta phase boundary');
    end
else
    error('varargin must either be empty or a structure specifying the fit conditions');
end

switch phase_system.data_type
    case 'FRET'
        if nC ~= ndata
            error('number of compositions must equal number of FRET data pts');
        end
        switch phase_system.name
            case 'binary'
                disp(sprintf('solving for phase boundary coordinates:\n'));
                options = optimset('Display','iter','Largescale','on','Jacobian','off','MaxFunEvals',1000,...
                    'MaxIter',100,'TolFun',10^-8,'TolX',10^-8);
                % x = [KpA KpD Ctraj_alpha Ctraj_beta];
                switch fit_conditions.taposition
                    case 'located'
                        switch fit_conditions.tbposition
                            case 'located'
                                ta = fit_conditions.taparameter;
                                tb = fit_conditions.tbparameter;
                                disp(sprintf('starting search and fitting...'));
                                tic;
                                switch fit_conditions.fit_method
                                    case 'lsqnonlin'
%                                         lb = [m-1e-8 -Inf]';
%                                         ub = [m+1e-8 Inf]';
                                        lb = [0.001 0.001 ta-1e-8 tb-1e-8]';
                                        ub = [1 1 ta+1e-8 tb+1e-8]';
                                        x0 = [0.06 0.8 ta tb]';
                                        [x,resnorm,residual,exitflag,output] =...
                                            lsqnonlin(@phase_boundary_fitting_fxn,x0,lb,ub,options,Ctraj,'FRET',C,data);
                                    case 'lsqcurvefit'
                                        lb = [0.001 0.001 ta-1e-8 tb-1e-8]';
                                        ub = [1000 1000 ta+1e-8 tb+1e-8]';
                                        x0 = [1 1 ta tb]';
                                        [x,resnorm,residual,exitflag,output] =...
                                            lsqcurvefit(@phase_boundary_fitting_fxn,x0,Ctraj,data,lb,ub,options,'FRET',C,data);
                                    case 'fmincon'
                                        lb = [0.001 0.001 ta tb]';
                                        ub = [1000 1000 ta tb]';
                                        x0 = [0.01 0.01 ta tb]';
                                        A = [];
                                        b = [];
                                        nonlcon = [];
                                        [x,resnorm,exitflag,output] =...
                                            fmincon(@phase_boundary_fitting_fxn,x0,A,b,Aeq,beq,lb,ub,nonlcon,options,varargin);
                                    case 'simplex'
                                        lb = [0.001 0.001 ta tb]';
                                        ub = [1 1 ta tb]';
                                        x0 = [0.001 1 ta tb]';
                                        A = [];
                                        b = [];
                                        nonlcon = [];
                                        [x,resnorm,exitflag,output] =...
                                            fminsearchcon(@phase_boundary_fitting_fxn,x0,lb,ub,A,b,nonlcon,options,Ctraj,'FRET',C,data);
                                    otherwise
                                        error('invalid fit method');
                                end
                                
%                                 x0 = [0.1+(10-0.1).*rand(1) 0.1+(10-0.1).*rand(1) 0.5.*rand(1) 0.5+(1-0.5).*rand(1)]';
%                                 x0 = [1 1 0.5.*rand(1) 0.5+(1-0.5).*rand(1)]';
                                
                                Ca = interp1(Ctraj,C,x(3));
                                Cb = interp1(Ctraj,C,x(4));
                                [fitmeasure,T,residual] = phase_boundary_fitting_fxn(x,Ctraj,'FRET',C,data);
                                fit_output = struct('KpA',x(1),'KpD',x(2),'ta',x(3),'tb',x(4),'C_alpha',Ca,'C_beta',...
                                    Cb,'resnorm',resnorm,'residual',residual,'exitflag',exitflag,'output',output,...
                                    't',Ctraj,'C',C,'data',data,'fit',T);
                                
                                disp(sprintf('\nsearch done...analyzing and displaying results'));
                                disp(sprintf('exitflag = %d',exitflag));
                                display = [resnorm;x];
                                disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                disp(sprintf('plotting data and fit as a function of composition:\n'));
                                figure,plot(Ctraj,data,'.k','markersize',15),hold on,plot(Ctraj,T,'-r'),hold off;
                                disp(sprintf('done!'));
                                toc;
                                return % fitting finished
                            case 'range'
                                switch fit_conditions.fit_method
                                    case 'lsqnonlin'
                                    case 'lsqcurvefit'
                                    case 'fmincon'
                                    case 'simplex'
                                    otherwise
                                        error('invalid fit method');
                                end
                            otherwise
                                error('invalid tb position');
                        end
                    case 'range'
                        switch fit_conditions.tbposition
                            case 'located'
                                switch fit_conditions.fit_method
                                    case 'lsqnonlin'
                                    case 'lsqcurvefit'
                                    case 'fmincon'
                                    case 'simplex'
                                    otherwise
                                        error('invalid fit method');
                                end
                            case 'range'
                                switch fit_conditions.fit_method
                                    case 'lsqnonlin'
                                    case 'lsqcurvefit'
                                    case 'fmincon'
                                    case 'simplex'
                                    otherwise
                                        error('invalid fit method');
                                end
                            otherwise
                                error('invalid tb position');
                        end
                    otherwise
                        error('invalid ta position');
                end
            case 'ternary2p'
                error('no code');
            case 'ternary3p'
                error('no code');
            otherwise
                error('invalid phase boundary description, for now');
        end    
    case 'spectra'
        if nC ~= nS
            error('number of compositions must equal number of spectra');
        end
        % check to see if spectra are derivative or absorbance spectra
        for s = 1:nS
            Bfield = data(:,s*2-1);
            approx_bline = data([1:50 end-50:end]',s*2); 
            % approx_bline = approximate baseline
            llim_noise = min(approx_bline); 
            % llim_noise = lower limit of baseline
            % noise = stdev of baseline
            if any(data(nBfield/2-100:nBfield/2+100,s*2) < llim_noise)
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
        switch phase_system.data_basis
            case 'svd'
                % the following code is modified from a previous m-file called
                % 'boundary_spectra_lsqlin.m'

                disp(sprintf('beginning the svd phase boundary fitting procedure\n'));
                disp(sprintf('performing spectra pre-processing:\n'));

                % align all spectra
                spectra = align_spectra(data);
                if aord == 'd'
                    % convert derivative spectra to absorbance spectra
                    spectra = deriv2abs(spectra);
                    aord = 'a';
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

                % set negative absorbances to zero
                % and   
                % renormalize spectra so that the sum of absorbances equals one
                for s = 1:ns
                    spectra(find(spectra(:,s*2) < 0),s*2) = 0;
                    spectra(:,s*2) = spectra(:,s*2)./sum(spectra(:,s*2));
                end

                % reduce spectra to just absorbance values 
                S = spectra(:,2:2:end);

                % get baseline statistics
                for s = 1:ns
                    Sbaseline(:,s) = S([1:200 end-200:end]',s);
                    std_Sbline(s) = std(Sbaseline(:,s));
                end 
                stdS = std(S')';

                % remove baseline or unimportant magnetic fields
                S = S(200:end,:);
        %         S = S(200:end-200,:);
        %         S = S([stdS >= mean(std_Sbline)],:);

                % perform singular value decomposition to obtain eigenspectra basis
                disp(sprintf('performing singular value decomposition to obtain eigenspectra basis:\n'));
                [U,W,V] = svd(S,0);
                B = U(:,1:phase_system.p); % B = p component eigenbasis

                % solve for eigenspectra coefficient matrix D by linear least squares
                disp(sprintf('solving for eigenspectra coefficient matrix by linear least squares:\n'));
                options = optimset('display','off','largescale','off');
                warning off;

                % constraints
        %         Aeq = [sum(U(:,1)) sum(U(:,2))];
        %         Aeq = sum(U(:,1:p));
                Aeq = sum(B);
                beq = sum(S);
                A = -B;
                b = zeros(size(U(:,1)));
                D = W(1:p,1:p)*V(:,1:p)';

        %         x0 = D;
                for s = 1:ns
                    disp(sprintf('\t fitting spectrum %d, %f percent done',s,((s-1)/ns)*100));
                    % [x,resnorm,residual] = lsqlin(B,S(:,s),A,b,Aeq,beq,[],[],x0(:,s)); % ...[],[],x0,options
                    [x,resnorm,residual] = lsqlin(B,S(:,s),A,b,Aeq,beq(s),[],[],[],options);
        %             [x,resnorm,residual] = lsqlin(B,S(:,s),A,b,Aeq,beq(s),[],[],x0(:,s),options);
        %             [x,resnorm,residual] = lsqlin(B,S(:,s),A,b,[],[],[],[],x0(:,s),options);
        %             x = lsqlin(B,S(:,s),[],[],Aeq,beq);
        %             x = lsqlin(B,S(:,s),A,b);
                    % x = B\S(:,s);
                    D_con(:,s) = x;
                end

                Sbar = B*D_con;
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
                    'U',U,'V',V,'W',W,'B',B,'Sbar',Sbar,'S',S);
        %         fit_output = struct('Kp',x(1),'C_alpha',x(2),'C_beta',x(3),'resnorm',resnorm,'eig_coeff_uncon',D,'eig_coeff_con',D_con,'theory_coeff',T,...
        %             'U',U,'V',V,'W',W,'B',B,'Sbar',Sbar,'S',S);

                disp(sprintf('plotting constrained coefficients of 1st eigenspectra (experiment) and fitted curve (theory) as a function of composition\n'));
                figure;
                plot(C,[D(:,1) D_con(:,1) T]);
                warning on;
                switch phase_system.name
                    case 'binary'
                        disp(sprintf('solving for phase boundary coordinates:\n'));
                        options = optimset('Display','iter','Largescale','on','Jacobian','off','MaxFunEvals',1000,...
                            'MaxIter',100,'TolFun',10^-8,'TolX',10^-8);
                        % x = [KpA KpD Ctraj_alpha Ctraj_beta];
                        switch fit_conditions.taposition
                            case 'located'
                                switch fit_conditions.tbposition
                                    case 'located'
                                        ta = fit_conditions.taparameter;
                                        tb = fit_conditions.tbparameter;
                                        disp(sprintf('starting search and fitting...'));
                                        tic;
                                        switch fit_conditions.fit_method
                                            case 'lsqnonlin'
        %                                         lb = [m-1e-8 -Inf]';
        %                                         ub = [m+1e-8 Inf]';
                                                lb = [0.001 0.001 ta-1e-8 tb-1e-8]';
                                                ub = [1 1 ta+1e-8 tb+1e-8]';
                                                x0 = [0.06 0.8 ta tb]';
                                                [x,resnorm,residual,exitflag,output] =...
                                                    lsqnonlin(@phase_boundary_fitting_fxn,x0,lb,ub,options,Ctraj,'FRET',C,data);
                                            case 'lsqcurvefit'
                                                lb = [0.001 0.001 ta-1e-8 tb-1e-8]';
                                                ub = [1000 1000 ta+1e-8 tb+1e-8]';
                                                x0 = [1 1 ta tb]';
                                                [x,resnorm,residual,exitflag,output] =...
                                                    lsqcurvefit(@phase_boundary_fitting_fxn,x0,Ctraj,data,lb,ub,options,'FRET',C,data);
                                            case 'fmincon'
                                                lb = [0.001 0.001 ta tb]';
                                                ub = [1000 1000 ta tb]';
                                                x0 = [0.01 0.01 ta tb]';
                                                A = [];
                                                b = [];
                                                nonlcon = [];
                                                [x,resnorm,exitflag,output] =...
                                                    fmincon(@phase_boundary_fitting_fxn,x0,A,b,Aeq,beq,lb,ub,nonlcon,options,varargin);
                                            case 'simplex'
                                                lb = [0.001 0.001 ta tb]';
                                                ub = [1 1 ta tb]';
                                                x0 = [0.001 1 ta tb]';
                                                A = [];
                                                b = [];
                                                nonlcon = [];
                                                [x,resnorm,exitflag,output] =...
                                                    fminsearchcon(@phase_boundary_fitting_fxn,x0,lb,ub,A,b,nonlcon,options,Ctraj,'FRET',C,data);
                                            otherwise
                                                error('invalid fit method');
                                        end

        %                                 x0 = [0.1+(10-0.1).*rand(1) 0.1+(10-0.1).*rand(1) 0.5.*rand(1) 0.5+(1-0.5).*rand(1)]';
        %                                 x0 = [1 1 0.5.*rand(1) 0.5+(1-0.5).*rand(1)]';

                                        Ca = interp1(Ctraj,C,x(3));
                                        Cb = interp1(Ctraj,C,x(4));
                                        [fitmeasure,T,residual] = phase_boundary_fitting_fxn(x,Ctraj,'FRET',C,data);
                                        fit_output = struct('KpA',x(1),'KpD',x(2),'ta',x(3),'tb',x(4),'C_alpha',Ca,'C_beta',...
                                            Cb,'resnorm',resnorm,'residual',residual,'exitflag',exitflag,'output',output,...
                                            't',Ctraj,'C',C,'data',data,'fit',T);

                                        disp(sprintf('\nsearch done...analyzing and displaying results'));
                                        disp(sprintf('exitflag = %d',exitflag));
                                        display = [resnorm;x];
                                        disp(sprintf('best fit measure = %f, with parameters [%f %f %f %f]\n',display));
                                        disp(sprintf('plotting data and fit as a function of composition:\n'));
                                        figure,plot(Ctraj,data,'.k','markersize',15),hold on,plot(Ctraj,T,'-r'),hold off;
                                        disp(sprintf('done!'));
                                        toc;
                                        return % fitting finished
                                    case 'range'
                                        switch fit_conditions.fit_method
                                            case 'lsqnonlin'
                                            case 'lsqcurvefit'
                                            case 'fmincon'
                                            case 'simplex'
                                            otherwise
                                                error('invalid fit method');
                                        end
                                    otherwise
                                        error('invalid tb position');
                                end
                            case 'range'
                                switch fit_conditions.tbposition
                                    case 'located'
                                        switch fit_conditions.fit_method
                                            case 'lsqnonlin'
                                            case 'lsqcurvefit'
                                            case 'fmincon'
                                            case 'simplex'
                                            otherwise
                                                error('invalid fit method');
                                        end
                                    case 'range'
                                        switch fit_conditions.fit_method
                                            case 'lsqnonlin'
                                            case 'lsqcurvefit'
                                            case 'fmincon'
                                            case 'simplex'
                                            otherwise
                                                error('invalid fit method');
                                        end
                                    otherwise
                                        error('invalid tb position');
                                end
                            otherwise
                                error('invalid ta position');
                        end
                    case 'ternary2p'
                        error('no code');
                    case 'ternary3p'
                        error('no code');
                    otherwise
                        error('invalid phase boundary description, for now');
                end
            case 'spectra'
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
    otherwise
        error('invalid data type');
end
return