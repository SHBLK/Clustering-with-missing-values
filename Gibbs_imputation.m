function [y] = Gibbs_imputation(x,obs,iter)
% Gibbs sampling imputation for multivariate normal distributed matrix


if ~exist('iter','var')
    iter = 200;
end
[n,d] = size(x);
y = x;
% Initialization
mu = (sum(y.*obs)./sum(obs))';
for i=1:d
    y(obs(:,i)==0,i) = mu(i);
end
Sigma = cov(y)+1e-3*eye(d);

idx_obs = obs>0;
% Priors
e = zeros(d,1);
Binv = eye(d);
nu = d+5;
S = eye(d);

idx = find(sum(obs,2)<d)';
for i=1:iter
    %% Imputation
    for j=idx
        temp = inv(Sigma(idx_obs(j,:),idx_obs(j,:)));
        mui = mu(~idx_obs(j,:)) + Sigma(~idx_obs(j,:),idx_obs(j,:))*temp*... 
                (y(j,idx_obs(j,:))'-mu(idx_obs(j,:)));
        sigmai = Sigma(~idx_obs(j,:),~idx_obs(j,:))-Sigma(~idx_obs(j,:),idx_obs(j,:))*temp*...
            Sigma(idx_obs(j,:),~idx_obs(j,:));
        try
            y(j,~idx_obs(j,:)) = mvnrnd(mui',sigmai);
        catch
            
            sigmai = (sigmai + sigmai')./2;
            y(j,~idx_obs(j,:)) = mvnrnd(mui',sigmai);
        end
    end
    %% parameter inference
    temp = y' - mu(:,ones(1,n));
    Sigma = iwishrnd(S+temp*temp', nu+n);
    sigmainv = inv(Sigma);
    A = inv(Binv+n*sigmainv);
    mu = mvnrnd((A*(sigmainv*sum(y)'+Binv*e))', A)';
end

end

