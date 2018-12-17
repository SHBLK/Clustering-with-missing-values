function [feature_vectors,S_,m_,v_,k_] = test_ConstructPriors(data_,calibSize)

% Function to calibrate priors from calibSize last columns of data.
% Also, discards features used in calibration and returns the feature
% vectors.

compSize=2;
remove_outliers = 1;%removes perc% largest values for sample moment calculations
perc=0.1;

[~,p_]=size(data_);
D=p_ - calibSize;
feature_vectors = data_(:,1:D);

data_calib = data_(:,(D + 1):end);




if remove_outliers==1
    Sigma_hat = cov(data_calib);
    sorted_mu_E=sort(mean(data_calib),'ascend');
    sorted_diag_Sigma_E=sort(diag(Sigma_hat),'ascend');
    E_new = floor((1-perc).*calibSize);
    mu_E = sorted_mu_E(1:E_new);
    diag_Sigma_E = sorted_diag_Sigma_E(1:E_new);
    E_mu1_hat = mean(mu_E);
    E_Sigma11_hat = mean(diag_Sigma_E);
    E_Sigma12_hat=sum(sum(Sigma_hat - diag(diag(Sigma_hat))))./(calibSize*(calibSize-1));
    Var_mu1_hat = (sum((mu_E-E_mu1_hat).^2))./(E_new-1);
    Var_Sigma11_hat = (sum((diag_Sigma_E-E_Sigma11_hat).^2))./(E_new-1);
else
    E_mu1_hat = mean(mean(data_calib));
    Sigma_hat = cov(data_calib);
    E_Sigma11_hat = mean(diag(Sigma_hat));
    E_Sigma12_hat=sum(sum(Sigma_hat - diag(diag(Sigma_hat))))./(calibSize*(calibSize-1));
    Var_mu1_hat = (sum((mean(data_calib)-E_mu1_hat).^2))./(calibSize-1);
    Var_Sigma11_hat = (sum((diag(Sigma_hat)-E_Sigma11_hat).^2))./(calibSize-1);
end



sigma2_ = (2.*E_Sigma11_hat).*( (((E_Sigma11_hat).^2)./(Var_Sigma11_hat)) + 1 );
rho_ = E_Sigma12_hat./E_Sigma11_hat;
m_ = ones(compSize,D).*E_mu1_hat;
v_ = ones(1,compSize).*(E_Sigma11_hat./Var_mu1_hat);
k_ = ones(1,compSize).*( ( (2.*((E_Sigma11_hat).^2))./(Var_Sigma11_hat) ) + D + 3);

SS= sigma2_ .* ( ((ones(D,D)-diag(ones(1,D))).*rho_)+  ( diag(ones(1,D)) ) );

S_=zeros(D,D,compSize);
for i=1:compSize
   S_(:,:,i) = SS; 
end

end

