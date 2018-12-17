function RandSigma=test_InverseWishartCovariance(k,S)
% test_InverseWishartCovariance generates two covaraince matrices using an
% Inverse Wishart distribution
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%
RandSigma(:,:,1) = iwishrnd(S(:,:,1),k(1,1));
RandSigma(:,:,2) = iwishrnd(S(:,:,2),k(1,2));
return