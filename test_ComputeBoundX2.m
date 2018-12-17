function val=test_ComputeBoundX2(ProbabiliesPartitionsSortedDescendent,ValidLabelFcnValuesSortedDescendent_phsi_i1,M,n,varargin)
% test_ComputeBoundX2 computes the value of the bound to compare with the
% maximum probability P1, in order to look for the optimal partition 
% using only a set of candidate partitions composed of those partitions with
% maximum Hamming distance iqual to or less than 2 from the partition P1 
% instead of all the partitions (Theorem 1 of the paper 'Bayes Labeling and 
% Bayes Clustering Operators for Known Random Labeled Point Processes')
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%
% See also test_ComputeBoundXi, test_ComputeBoundX1
suma=0;
for j=1:2
    IdxOfLabelFcnsHammingDistj=test_FindLabelFcnsHammingDistance(ValidLabelFcnValuesSortedDescendent_phsi_i1,ValidLabelFcnValuesSortedDescendent_phsi_i1(1,:),j,varargin);
    suma=suma+(M*n+2-2*j+1)*sum(ProbabiliesPartitionsSortedDescendent(IdxOfLabelFcnsHammingDistj));
    clear IdxOfLabelFcnsHammingDistj;
end
val=(M*n-suma)/(M*n+2+1);
return