function [Prphi1_S,Prphi2_S]=test_ComputeProbabilityNormalMeanInverseWishartSigma_withMiss(Vectors,HypParam,LabelFcnValues,Obs,varargin)


m1 = HypParam.templates(1,:)';
m2 = HypParam.templates(2,:)';
[p q]=size(HypParam.S);
S1 = HypParam.S(1:p/2,:);
S2 = HypParam.S(p/2+1:end,:);
v1 = HypParam.v(1,1);
v2 = HypParam.v(1,2);
k1 = HypParam.k(1,1);
k2 = HypParam.k(1,2);
d = size(m1,2);
J = 1000;                      % number of draws for MC approximation

NumPartitions = size(LabelFcnValues.phsi_i1,1);

% Computing Conditional Probabilities for every possible label function
Prphi1_S = zeros(NumPartitions,1);
Prphi2_S = zeros(NumPartitions,1);
for j=1:NumPartitions
  if j==1
   
   NumPartitions

  end

    % Partitioning the set of points using the label function phsi_i1
    Labels1=LabelFcnValues.phsi_i1(j,:)'+1;
    
    % Class 1
    x1 = Vectors(Labels1==1,:);     % n1 by d matrix of datapoints in cluster 1
%    n1 = size(x1,1);
    idx_obs1 = Obs(Labels1==1,:);    % observed features indicator matrix of cluster 1
    % find patterns of missing data in 1st cluster
    [g1,~,ig1] = unique(idx_obs1,'rows','stable');
    ng1 = full(diag(sparse(ig1,ig1,1)));
    G1 = length(ng1);
    mg1 = cell(G1,1);
    Sg1 = cell(G1,1);
    for gg=1:G1
        idx = g1(gg,:)==1;
        temp = x1(ig1==gg,idx);
        if ng1(gg)>1
            mg1{gg} = mean(temp)';
            Sg1{gg} = (ng1(gg)-1)*cov(temp);
        elseif ng1(gg)==1
            mg1{gg} = temp';
            Sg1{gg} = zeros(length(temp));
        else
            error('Problem');
        end
    end
    
    % Class 2
    x2 = Vectors(Labels1==2,:);     % n2 by d matrix of datapoints in cluster 2
%    n2 = size(x2,1);
    idx_obs2 = Obs(Labels1==2,:);    % observed features indicator matrix of cluster 2
    % find patterns of missing data in 1st cluster
    [g2,~,ig2] = unique(idx_obs2,'rows','stable');
    ng2 = full(diag(sparse(ig2,ig2,1)));
    G2 = length(ng2);
    mg2 = cell(G2,1);
    Sg2 = cell(G2,1);
    for gg=1:G2
        idx = g2(gg,:)==1;
        temp = x2(ig2==gg,idx);
        if ng2(gg)>1
            mg2{gg} = mean(temp)';
            Sg2{gg} = (ng2(gg)-1)*cov(temp);
        elseif ng2(gg)==1
            mg2{gg} = temp';
            Sg2{gg} = zeros(length(temp));
        else
            error('Problem');
        end
    end
    
    prob1_S1 = zeros(J,1);
    prob2_S1 = zeros(J,1);
    prob1_S2 = zeros(J,1);
    prob2_S2 = zeros(J,1);
    for jj=1:J
        Sigma1 = iwishrnd(S1,k1);
        Sigma2 = iwishrnd(S2,k2);
        prob1_S1(jj) = calc_int(G1,g1,mg1,Sg1,ng1,v1,m1,Sigma1);
        prob2_S1(jj) = calc_int(G1,g1,mg1,Sg1,ng1,v2,m2,Sigma2);
        prob1_S2(jj) = calc_int(G2,g2,mg2,Sg2,ng2,v1,m1,Sigma1);
        prob2_S2(jj) = calc_int(G2,g2,mg2,Sg2,ng2,v2,m2,Sigma2);
    end
    Prphi1_S(j) = log_sumExp(prob1_S1)+log_sumExp(prob2_S2)-2*log(J);
    Prphi2_S(j) = log_sumExp(prob2_S1)+log_sumExp(prob1_S2)-2*log(J);

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