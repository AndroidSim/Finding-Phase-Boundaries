function [V,D] = phase_boundary_fit_fun(x,C,basis_method,varargin)

if nargin > 3
    switch basis_method
        case 'svd'
            data = varargin{1};
%             C = varargin{2};
        case 'spectral'
            if any(strcmp(varargin{1},{'grid';'cc'}))
                search_method = varargin{1};
            else
                error('search method for spectral basis method must be "grid" or "cc"');
            end

            if nargin > 4
                if any(strcmp(varargin{2},{'coeff';'spectra'}))
                    data_type = varargin{2};
                else
                    error('data type for spectral basis method must be "coeff" or "spectra"');
                end
            end
        otherwise
            error('invalid basis_method in phase_boundary_fit_fun');
    end
end
    
switch basis_method
    case 'svd'
        % x = [kp a_alpha a_beta c_alpha c_beta]
        nc = length(C);
        V = repmat(0,nc,1);
        
%         data_alpha = interp1(C,data,x(2));
%         data_beta = interp1(C,data,x(3));
%         for c = 1:nc
%             if C(c) <= x(2)
%                 V(c) = data_alpha;
%             elseif C(c) > x(2) & C(c) < x(3)
%                 V(c) = ((x(1).*(x(3)-C(c)))./(C(c)-x(2)+(x(1).*(x(3)-C(c))))).*(data_alpha-data_beta) + data_beta;
%             elseif C(c) >= x(3)
%                 V(c) = data_beta;
%             end
%         end
        
        for c = 1:nc
            if C(c) <= x(4)
                V(c) = x(2);
            elseif C(c) > x(4) & C(c) < x(5)
                fb = (C(c)-x(4))./(x(5)-x(4));
                fa = 1-fb;
                fpb = (x(1).*fb)./(fa + x(1).*fb);
                fpa = 1-fpb;
                V(c) = x(2).*fpa + x(3).*fpb;
%                 V(c) = ((x(1).*(x(5)-C(c)))./(C(c)-x(4)+(x(1).*(x(5)-C(c))))).*(x(2)-x(3)) + x(3);
            elseif C(c) >= x(5)
                V(c) = x(3);
            end
        end
        
%         V = norm(V-data)/std(data);
        
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

return