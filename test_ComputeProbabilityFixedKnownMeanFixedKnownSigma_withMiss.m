function [Prphi1_S,Prphi2_S]=test_ComputeProbabilityFixedKnownMeanFixedKnownSigma_withMiss(Vectors,HypParam,LabelFcnValues,Obs,varargin)


Mu1 = HypParam.templates(1,:)';
Mu2 = HypParam.templates(2,:)';
p = size(HypParam.variances,1);
Sigma1 = HypParam.variances(1:p/2,:);
Sigma2 = HypParam.variances(p/2+1:end,:);
%d = size(Mu1,2);
NumPartitions = size(LabelFcnValues.phsi_i1,1);

% Computing Conditional Probabilities for every possible label function
Prphi1_S=[];
Prphi2_S=[];
for j=1:NumPartitions
    % Partitioning the set of points using the label function phsi_i1
    Labels1 = LabelFcnValues.phsi_i1(j,:)'+1;
    
    % Class 1
    x1 = Vectors(Labels1==1,:);     % n1 by d matrix of datapoints in cluster 1
%    n1 = size(x1,1);
    idx_obs1 = Obs(Labels1==1,:);    % observed features indicator matrix of cluster 1
    % find patterns of missing data in 1st cluster
    [g1,~,ig1] = unique(idx_obs1,'rows','stable');
    ng1 = full(diag(sparse(ig1,ig1,1)));
    G1 = length(ng1);
%    mg1 = cell(G1,1);
%    Sg1 = cell(G1,1);
    logdet1_S1 = 0;
    logdet2_S1 = 0;
    trc1_S1 = 0;
    trc2_S1 = 0;
    quad1_S1 = 0;
    quad2_S1 = 0;
    for gg=1:G1
        idx = g1(gg,:)==1;
        temp = x1(ig1==gg,idx);
        if ng1(gg)>1
            mg = mean(temp)';
            Sg = (ng1(gg)-1)*cov(temp);
        elseif ng1(gg)==1
            mg = temp';
            Sg = zeros(length(temp));
        else
            error('Problem');
        end
        logdet1_S1 = logdet1_S1 + ng1(gg)*log(max(det(2*pi*Sigma1(idx,idx)),realmin));
        logdet2_S1 = logdet2_S1 + ng1(gg)*log(max(det(2*pi*Sigma2(idx,idx)),realmin));
        trc1_S1 = trc1_S1 + trace(Sg/Sigma1(idx,idx));
        trc2_S1 = trc2_S1 + trace(Sg/Sigma2(idx,idx));
        temp = mg-Mu1(idx);
        quad1_S1 = quad1_S1 + ng1(gg)*temp'*inv(Sigma1(idx,idx))*temp;
        temp = mg-Mu2(idx);
        quad2_S1 = quad2_S1 + ng1(gg)*temp'*inv(Sigma2(idx,idx))*temp;
    end
    
    % Class 2
    x2 = Vectors(Labels1==2,:);     % n2 by d matrix of datapoints in cluster 2
%    n2 = size(x2,1);
    idx_obs2 = Obs(Labels1==2,:);    % observed features indicator matrix of cluster 2
    % find patterns of missing data in 1st cluster
    [g2,~,ig2] = unique(idx_obs2,'rows','stable');
    ng2 = full(diag(sparse(ig2,ig2,1)));
    G2 = length(ng2);
%    mg2 = cell(G2,1);
%    Sg2 = cell(G2,1);
    logdet1_S2 = 0;
    logdet2_S2 = 0;
    trc1_S2 = 0;
    trc2_S2 = 0;
    quad1_S2 = 0;
    quad2_S2 = 0;
    for gg=1:G2
        idx = g2(gg,:)==1;
        temp = x2(ig2==gg,idx);
        if ng2(gg)>1
            mg = mean(temp)';
            Sg = (ng2(gg)-1)*cov(temp);
        elseif ng2(gg)==1
            mg = temp';
            Sg = zeros(length(temp));
        else
            error('Problem');
        end
        logdet1_S2 = logdet1_S2 + ng2(gg)*log(max(det(2*pi*Sigma1(idx,idx)),realmin));
        logdet2_S2 = logdet2_S2 + ng2(gg)*log(max(det(2*pi*Sigma2(idx,idx)),realmin));
        trc1_S2 = trc1_S2 + trace(Sg/Sigma1(idx,idx));
        trc2_S2 = trc2_S2 + trace(Sg/Sigma2(idx,idx));
        temp = mg-Mu1(idx);
        quad1_S2 = quad1_S2 + ng2(gg)*temp'*inv(Sigma1(idx,idx))*temp;
        temp = mg-Mu2(idx);
        quad2_S2 = quad1_S2 + ng2(gg)*temp'*inv(Sigma2(idx,idx))*temp;
    end
    
    % log Probabilities for Label Function 1
    Prphi1_S(j,1) = -0.5*(logdet1_S1+trc1_S1+quad1_S1...
        + logdet2_S2+trc2_S2+quad2_S2);
    
    % log Probabilities for Label Function 2
    Prphi2_S(j,1) = -0.5*(logdet2_S1+trc2_S1+quad2_S1...
        + logdet1_S2+trc1_S2+quad1_S2);
end


if nargin==5
    try
        aux=varargin{1};
    end
    if isempty(aux)==1
        varargin{1}='Normalize';
    else
        if iscell(aux)==1
            varargin=varargin{1};
        else
            varargin{1}=aux;
        end
    end
end
if nargin==4
    varargin{1}='Normalize';
end
% Normalize probabilities
if strcmp(varargin{1},'Normalize')==1
mp = max([Prphi1_S;Prphi2_S]);
Prphi1_S = exp(Prphi1_S-mp);
Prphi2_S = exp(Prphi2_S-mp);
Den = sum(Prphi1_S+Prphi2_S);
Prphi1_S = Prphi1_S/Den;
Prphi2_S = Prphi2_S/Den;
else
Prphi1_S = exp(Prphi1_S);
Prphi2_S = exp(Prphi2_S);   
end
end
