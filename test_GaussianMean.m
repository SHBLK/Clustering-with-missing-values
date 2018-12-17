function RandMu=test_GaussianMean(Centers,Sigma)
% Developed by Marco Benalcázar - September 25/2012
%
% test_GaussianMean  generates random means drawn from a Normal Distribution
%
% RandMU = test_GaussianMean(Centers,Sigma)
%
% INPUTS:
%       
%   Centers          - LxD. L is the number of classes. D is the dimention
%                      of the working space. Each row of Centers contains 
%                      the mean for class Li.
%   Sigma            - DxDxL. Sigma is an hypermatrix where each layer
%                      is a DxD matrix containing the covariance matrix 
%                      for class Li.
%
% OUTPUTS:
%
%   RandMu           - LxD. Each row of RandMU(i,:) contains the random mean
%                      drawn from a Gaussian with mean Centers(i,:) and 
%                      Covariance Sigma(:,:,i)
RandMu = mvnrnd(Centers,Sigma);
return
