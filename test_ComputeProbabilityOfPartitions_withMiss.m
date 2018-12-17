function [Prphi1_S,Prphi2_S]=test_ComputeProbabilityOfPartitions_withMiss(TypeModel,Vectors,HypParam,Partitions,MissMatr,varargin)
% test_ComputeProbabilityOfPartitions_withMiss computes the probability of each
% vector of labels of the list 'Partitions' using the model from which the
% hyperparameteres 'HypParam' come
%

%
% See also test_ComputeProbabilityFixedKnownMeanFixedKnownSigma, 
% test_ComputeProbabilityNormalMeanFixedKnownSigma, and
% test_ComputeProbabilityNormalMeanInverseWishartSigma

%
if strcmp(TypeModel,'EXP_FixedKnownMean_FixedKnownCovariance')
    [Prphi1_S,Prphi2_S]=test_ComputeProbabilityFixedKnownMeanFixedKnownSigma_withMiss(Vectors,HypParam,Partitions,MissMatr,varargin);
elseif strcmp(TypeModel,'EXP_NormalMean_FixedKnownCovariance')
    [Prphi1_S,Prphi2_S]=test_ComputeProbabilityNormalMeanFixedKnownSigma_withMiss(Vectors,HypParam,Partitions,MissMatr,varargin);
elseif strcmp(TypeModel,'EXP_NormalMean_InverseWishartCovariance')
    [Prphi1_S,Prphi2_S]=test_ComputeProbabilityNormalMeanInverseWishartSigma_withMiss(Vectors,HypParam,Partitions,MissMatr,varargin);
end
return