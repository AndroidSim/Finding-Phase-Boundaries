function varargout = phase_boundary_fitting_fxn(x,u,data_type,varargin)

if nargin > 3
    switch data_type
        case 'FRET'
            C = varargin{1};
            F = varargin{2};
        case 'spectra'
            basis_method = varargin{1};
            switch basis_method
                case 'svd'
                    C = varargin{2};
                case 'spectra'
                    fit_data = varargin{2};
                    switch fit_data
                        case 'coeff'
                            C = varargin{3};
                            coeff = varargin{4};
                        case 'spectra'
                            C = varargin{3};
                            S = varargin{4};
                        otherwise 
                            error('invalid data to be fit for spectral basis method');
                    end
                otherwise
                    error('invalid basis method');
            end
        otherwise
            error('invalid data type');
    end
end

switch data_type
    case 'FRET'
        % x = [KpA KpD ua ub];
        KpA = x(1);
        KpD = x(2);
        ua = x(3);
        ub = x(4);
        [nC,nCdim] = size(C);
        Fa = interp1(u,F,ua);
        Fb = interp1(u,F,ub);
        Ca = interp1(u,C,ua);
        Cb = interp1(u,C,ub);
        V = repmat(0,nC,1); 
        V(u <= ua) = F(u <= ua);
%         V(u <= ua) = Fa;
%         V(u <= ua) = NaN;
        phib = (C(u > ua & u < ub)-Ca)./(Cb-Ca);
        phia = 1-phib;
        V(u > ua & u < ub) = Fa*(phia./((phia + KpD.*phib).*(phia + KpA.*phib))) + Fb*(KpD.*KpA.*phib./((phia + KpD.*phib).*(phia + KpA.*phib)));
        V(u >= ub) = F(u >= ub);
%         V(u >= ub) = Fb;
%         V(u >= ub) = NaN;
%         for c = 1:nC
%             if u(c) <= ua
%                 V(c) = F(c);
%             elseif u(c) > ua & u(c) < ub
%                 phib = norm(C(c)-Ca)/norm(Cb-Ca);
%                 phia = 1-phib;
%                 V(c) = Fa*(phia/((phia + KpD*phib)*(phia + KpA*phib))) + Fb*(KpD*KpA*phib/((phia + KpD*phib)*(phia + KpA*phib)));
%             elseif u(c) >= ub
%                 V(c) = F(c);
%             end
%         end 
        ndatafit = nnz(u > ua & u < ub);
        temp = V-F;
        residual = temp(~isnan(temp));
        if nargout > 1
%             varargout{1} = (V-F)./sqrt(ndatafit);
            varargout{1} = mean(residual.^2);
            varargout{2} = V;
            varargout{3} = residual;
        else
%             varargout{1} = residual./sqrt(ndatafit);
%             varargout{1} = residual;
%             varargout{1} = V;
            varargout{1} = sum(residual.^2);
        end
    case 'spectra'
        switch basis_method
            case 'svd'
                % x = [kp a_alpha a_beta c_alpha c_beta]
                nc = length(C);
                V = repmat(0,nc,1);

                for c = 1:nc
                    if C(c) <= x(4)
                        V(c) = x(2);
                    elseif C(c) > x(4) & C(c) < x(5)
                        V(c) = ((x(1).*(x(5)-C(c)))./(C(c)-x(4)+(x(1).*(x(5)-C(c))))).*(x(2)-x(3)) + x(3);
                    elseif C(c) >= x(5)
                        V(c) = x(3);
                    end
                end

                if nargout > 1
                    kp = x(1);
                    a = x(4);
                    b = x(5);
                    d1 = x(2);
                    d2 = x(3);
                    D = zeros(nc,5);

                    for ce = 1:nc
                        c = C(ce);

                        if c <= x(4)
                            D(ce,:) = [0 1 0 0 0];
                        elseif c > x(4) & c < x(5)
                            D(ce,:) = [(b-c)/(c-a+kp*(b-c))*(d1-d2)-kp*(b-c)^2/(c-a+kp*(b-c))^2*(d1-d2),...
                                        kp*(b-c)/(c-a+kp*(b-c)),...
                                        -kp*(b-c)/(c-a+kp*(b-c))+1,...
                                        kp*(b-c)/(c-a+kp*(b-c))^2*(d1-d2),...    
                                        kp/(c-a+kp*(b-c))*(d1-d2)-kp^2*(b-c)/(c-a+kp*(b-c))^2*(d1-d2)];
                        elseif c >= x(5)
                            D(ce,:) = [0 0 1 0 0];
                        end
                    end
                end
            case 'spectral'
                switch search_method
                    case 'grid'
                        switch data_type
                            case 'coeff'
                                % varargin = search_method,data_type,alpha,beta;
                                % x = kp;
                                alpha = varargin{3};
                                beta = varargin{4};
                                nc = length(C);
                                V = repmat(0,nc,1);

                                for c = 1:nc
                                    if C(c) <= alpha
                                        V(c) = 1;
                                    elseif C(c) > alpha & C(c) < beta
                                        V(c) = (x.*(beta-C(c)))./(C(c)-alpha+(x.*(beta-C(c))));
                                    elseif C(c) >= beta
                                        V(c) = 0;
                                    end
                                end
                            case 'spectra'
                                % varargin = search_method,data_type,alpha,beta,B;
                                % x = kp;
                                alpha = varargin{3};
                                beta = varargin{4};
                                B = varargin{5};
                                V = ((x.*(beta-C))./(C-alpha+(x.*(beta-C)))).*B(:,1) + ((C-alpha)./(C-alpha+(x.*(beta-C)))).*B(:,2);
                            otherwise
                                error('invalid data type for spectral search method in fit function');
                        end
                    case 'cc'
                        switch data_type
                            case 'coeff'
                                % varargin = search_method,data_type,S
                                % x = [alpha beta]
                                S = varargin{3};
                                % interpolate spectra matrix S to find the basis spectra at x(1) = alpha
                                % and x(2) = beta
                                [nb,ns] = size(S);
                                Ci = [x(1) x(2)];

                                for b = 1:nb
                                    A = S(b,:);
                                    Ai = interp1(C,A,Ci); % default = linear
                                    B(b,:) = Ai;
                                end

                                ops = optimset('display','off');
                                % constraints
                                Aeq = [1 1];
                                beq = 1;
                                A = [-1 0;0 -1];
                                b = [0;0];

                                for s = 1:ns
                                    p(:,s) = lsqlin(B,S(:,s),A,b,Aeq,beq,[],[],[],ops);% x = lsqlin(C,d,A,b,Aeq,beq,lb,ub,x0,options)
                                end

                                % generate spectra from coefficients of basis B
                                for s = 1:ns
                                    T(:,s) = p(1,s).*B(:,1) + p(2,s).*B(:,2);
                                end

                                V = reshape(T,nb*ns,1);
                            case 'spectra'
                                % varargin = search_method,data_type,S
                                % x = [kp alpha beta]

                                S = varargin{3};
                                [nb,ns] = size(S);
                                nc = ns;

                                temp = [x(2) x(3)]';
                                low = temp < C(1);
                                up = temp > C(end);

                                if any(low)
                                    temp(low) = C(1);
                                end

                                if any(up)
                                    temp(up) = C(end);
                                end

                                x(2:3) = temp;

                                % check compositions with C_alpha and C_beta to insure no generation of NaN
                                for c = 1:nc
                                    if x(3) == (C(c)*(x(1)-1)+x(2))/x(1)
                                        if x(2) == x(3) & x(2) == C(1)
                                            x(3) = x(2)+0.0001;
                                        elseif x(2) == x(3) & x(2) == C(end)
                                            x(2) = x(3)-0.0001;
                                        elseif x(2) == x(3) & x(2) == C(c)
                                            x(3) = x(2)+0.0001;
                                        else
                                            x(1) = x(1)+0.0001;
                                        end
                                    else
                                        if x(2) == x(3) & x(2) == C(1)
                                            x(3) = x(2)+0.0001;
                                        elseif x(2) == x(3) & x(2) == C(end)
                                            x(2) = x(3)-0.0001;
                                        elseif x(2) == x(3) & x(2) == C(c)
                                            x(3) = x(2)+0.0001;
                                        end
                                    end
                                end

                                % interpolate spectra matrix S to find the basis spectra at x(2) = alpha
                                % and x(3) = beta
                                Ci = [x(2) x(3)];

                                for b = 1:nb
                                    A = S(b,:);
                                    Ai = interp1(C,A,Ci); % default = linear
                                    B(b,:) = Ai;
                                end

                                % generate theoretical spectra from basis B
                                for c = 1:nc
                                    t_spec(:,c) = ((x(1).*(x(3)-C(c)))./(C(c)-x(2)+(x(1).*(x(3)-C(c))))).*B(:,1) + ((C(c)-x(2))./(C(c)-x(2)+(x(1).*(x(3)-C(c))))).*B(:,2);
                                end

                                V = reshape(t_spec,nb*ns,1);

                                if nargout > 1
                                    D = [0 0 0];

                                    for s = 1:ns
                                        % calculate gradient D
                                        kp = x(1);
                                        a = x(2);
                                        b = x(3);
                                        c = C(s);
                                        Sa = B(:,1);
                                        Sb = B(:,2);
                                        dV = [((b-c)/(c-a+kp*(b-c))).*Sa-(kp*(b-c)^2/(c-a+kp*(b-c))^2).*Sa-((c-a)*(b-c)/(c-a+kp*(b-c))^2).*Sb,...             
                                              (kp*(b-c)/(c-a+kp*(b-c))^2).*Sa-(1/(c-a+kp*(b-c))).*Sb+((c-a)/(c-a+kp*(b-c))^2).*Sb,...      
                                              (kp/(c-a+kp*(b-c))).*Sa-(kp^2*(b-c)/(c-a+kp*(b-c))^2).*Sa-((c-a)*kp/(c-a+kp*(b-c))^2).*Sb];
                                        D = D+sum(dV);
                                    end

                                    D = D';
                                end
                            otherwise
                                error('invalid data type for spectral search method in fit function');
                        end
            otherwise
                error('invalid search method for spectral basis method in fit function');
        end
            otherwise
                error('invalid basis method in fit function');
        end
    otherwise
        error('invalid data type');
end

return