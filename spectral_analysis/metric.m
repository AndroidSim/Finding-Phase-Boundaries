function m = metric(spec_A,spec_B,varargin)
% metric calculates a distance or metric between two spectra, where the
% type of metric is given as in 1st element in varargin
%
% metric type can be any of the following strings:
%   see "pdist" matlab fxn
%   'c_of_d'
%   'strain_energy'
%   'mismatch_mat'
%   'eign_basis'

[np,nc] = size(spec_A);

if nc < 2 || rem(nc,2) ~= 0
    error('1st spectrum must contain 2 columns: [B-field absorbance_values]');
end

[np,nc] = size(spec_B);

if nc < 2 || rem(nc,2) ~= 0
    error('2nd spectrum must contain 2 columns: [B-field absorbance_values]');
end

type = varargin{1};

if ~ischar(type)
    error('the input argument type must be a string specifying metric type');
end

SA = spec_A;
SB = spec_B;
B0 = 3326;

switch type
    case 'euclidean'
        a_spec = align_spectra([SA SB],[B0 0]);
        ai = find(~any(isnan(a_spec),2));
        data = [SA(ai,2) SB(ai,2)]';
        m = pdist(data,type);
    case 'seuclidean'
        a_spec = align_spectra([SA SB],[B0 0]);
        ai = find(~any(isnan(a_spec),2));
        data = [SA(ai,2) SB(ai,2)]';
        m = pdist(data,type);
    case 'mahalanobis'
        a_spec = align_spectra([SA SB],[B0 0]);
        ai = find(~any(isnan(a_spec),2));
        data = [SA(ai,2) SB(ai,2)]';
        m = pdist(data,type);
    case 'cityblock'
        a_spec = align_spectra([SA SB],[B0 0]);
        ai = find(~any(isnan(a_spec),2));
        data = [SA(ai,2) SB(ai,2)]';
        m = pdist(data,type);
    case 'minkowski'
        p = varargin{2};
        a_spec = align_spectra([SA SB],[B0 0]);
        ai = find(~any(isnan(a_spec),2));
        data = [SA(ai,2) SB(ai,2)]';
        m = pdist(data,type,p);
    case 'cosine'
        a_spec = align_spectra([SA SB],[B0 0]);
        ai = find(~any(isnan(a_spec),2));
        data = [SA(ai,2) SB(ai,2)]';
        m = pdist(data,type);
    case 'correlation'
        a_spec = align_spectra([SA SB],[B0 0]);
        ai = find(~any(isnan(a_spec),2));
        data = [SA(ai,2) SB(ai,2)]';
        m = pdist(data,type);
    case 'hamming'
        a_spec = align_spectra([SA SB],[B0 0]);
        ai = find(~any(isnan(a_spec),2));
        data = [SA(ai,2) SB(ai,2)]';
        m = pdist(data,type);
    case 'jaccard'
        a_spec = align_spectra([SA SB],[B0 0]);
        ai = find(~any(isnan(a_spec),2));
        data = [SA(ai,2) SB(ai,2)]';
        m = pdist(data,type);
    case 'c_of_d'
        a_spec = align_spectra([SA SB],[B0 0]);
        ai = find(~any(isnan(a_spec),2));
        B_field = a_spec(ai,1);
        data = [SA(ai,2) SB(ai,2)];
        deltas = abs(diff(data,1,2));
        m = (sum(B_field.*deltas))/sum(deltas);
    case 'strain_energy'
    case 'mismatch_mat'
        % the following procedure and metric was adopted from the following
        % paper:
        %   Bilge Gunsel and A. Murat Tekalp.  "Shape Similarity Matching
        %   for Query-by-Example"  Pattern Recognition.  vol. 31 no. 7 pp. 931-944, 1998
        
        sigma = varargin{2};
        
        shape_mat_A = spectral_prop(SA,'d','shape_mat');
        shape_mat_B = spectral_prop(SB,'d','shape_mat');
        %HA = exp((shape_mat_A.^2)/2*sigma^2);
        %HB = exp((shape_mat_B.^2)/2*sigma^2);
        HA = shape_mat_A;
        HB = shape_mat_B;
        [VA,DA] = eig(HA);
        [VB,DB] = eig(HB);
        
        % create mismatch matrix Z between A and B
        for i = 1:np
            for j = 1:np
                Z(i,j) = (norm(VA(:,i)-VB(:,j))).^2;
            end
        end
        
        % find matching feature or boundary points
        K=0;
        Kmax = np;
        
        for i = 1:np
            mini = min(Z(i,:));
            jmini = find(mini == Z(i,:));
            minjmini = min(Z(:,jmini));
            iminjmini = find(minjmini == Z(:,jmini));
            
            if iminjmini == i || minjmini == mini
                K = K+1;
                % mp == matched points, 1st column = index of spec_A and
                % 2nd column = index of spec_B
                mp(K,:) = [i jmini];
                s(K) = mini;
            end
        end
        
        % calculate similarity metric
        if K == Kmax
            m = (norm(s).^2)/K;
        else
            m = ((norm(s).^2)+2*(Kmax-K))/Kmax;
        end
    case 'eign_basis'
        % the following procedure and metric i created following the
        % procedure and metric for 'mismatch_mat'
        
        %sigma = varargin{2};
        
        shape_mat_A = spectral_prop(SA,'d','shape_mat');
        shape_mat_B = spectral_prop(SB,'d','shape_mat');
        %HA = exp((shape_mat_A.^2)/2*sigma^2);
        %HB = exp((shape_mat_B.^2)/2*sigma^2);
        HA = shape_mat_A;
        HB = shape_mat_B;
        [VA,DA] = eig(HA);
        [VB,DB] = eig(HB);
        
        % create mismatch matrix Z between A and B
        Z = VA'*VB;
        
        % find matching feature or boundary points
        K=0;
        
        for i = 1:np
            maxi = max(Z(i,:));
            jmaxi = find(maxi == Z(i,:));
            maxjmaxi = max(Z(:,jmaxi));
            imaxjmaxi = find(maxjmaxi == Z(:,jmaxi));
            
            if imaxjmaxi == i || maxjmaxi == maxi
                K = K+1;
                % mp == matched points, 1st column = index of spec_A and
                % 2nd column = index of spec_B
                mp(K,:) = [i jmaxi];
                s(K) = maxi;
            end
        end
        
        % calculate similarity metric
        m = norm(s)/length(s); 
        d = 1-diag(Z);
        m = sum(d)/length(d);
    otherwise
        error('no metric type was given');
end
return        