function [Prphi1_S,Prphi2_S,normI]=test_ComputeProbabilityOfPartitions_withMiss_Real(TypeModel,Vectors,HypParam,Partitions,MissMatr,flag_normVal,varargin)
% test_ComputeProbabilityOfPartitions computes the probability of each
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
    [Prphi1_S,Prphi2_S,normI]=test_ComputeProbabilityNormalMeanInverseWishartSigma_withMiss_R(Vectors,HypParam,Partitions,MissMatr,flag_normVal,varargin);
end
return