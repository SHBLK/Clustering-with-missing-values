function [FeatureVectors Labels]=test_GaussianSamples(Mu,Sigma,NumberSamples)
% Developed by Marco Benalcázar - September 25/2012
%
% test_GaussianSamples  generates random samples (feature vectors and
%                       labels) drawn from a mixture of 2 Gaussians
%
% [FeatureVectors Labels] = test_GaussianSamples(Mu,Sigma,NumberSamples)
%
% INPUTS:
%       
%   Mu               - LxD. L is the number of classes. D is the dimention
%                      of the working space. Each row of Mu contains 
%                      the mean for class Li.
%   Sigma            - DxDxL. Sigma is an hypermatrix where each layer
%                      is a DxD matrix containing the covariance matrix 
%                      for class Li.
%   NumberSamples    - 1x2. NumberSamples is a row vector where the first
%                      column, Numbersamples(1,1), is the number of points
%                      to be drawn from the Gaussian 1 and the second
%                      column, Numbersamples(1,2), is the number of points
%                      from the Gaussian 2.
%
% OUTPUTS:
%
%   FeatureVectors   - [Numbersamples(1,1)+Numbersamples(1,2)]xD. Each row 
%                      of FeatureVectors(i,:) contains a random feature 
%                      vector drawn from the 2 Gaussians Model.  
%   Labels           - 1x[Numbersamples(1,1)+Numbersamples(1,2)]. Labels is 
%                      a row vector where each column, Labels(1,i), is the
%                      true label for the random vector
%                      FeatureVectors(i,:).
Sigma_NewFormat=[];
for i=1:size(Sigma,3)
    Sigma_NewFormat{i}=Sigma(:,:,i);
end

[FeatureVectors Labels] = bb_cl_simdata(Mu,Sigma_NewFormat,NumberSamples,1,'ga');
return