function val=test_ComputeBoundX1(ProbabiliesPartitionsSortedDescendent,ValidLabelFcnValuesSortedDescendent_phsi_i1,M,n,varargin)
% test_ComputeBoundX2 computes the value of the bound to compare with the
% maximum probability P1, in order to determine that the optimal
% partition for a given point set is the partition with the probability
% P1 (Theorem 1 of the paper 'Bayes Labeling and Bayes Clustering Operators 
% for Known Random Labeled Point Processes')
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%
% See also test_ComputeBoundXi, test_ComputeBoundX2
IdxOfLabelFcnsHammingDist1=test_FindLabelFcnsHammingDistance(ValidLabelFcnValuesSortedDescendent_phsi_i1,ValidLabelFcnValuesSortedDescendent_phsi_i1(1,:),1,varargin);  
val=((M*n)/(M*n+2))*(1-sum(ProbabiliesPartitionsSortedDescendent(IdxOfLabelFcnsHammingDist1)));
return