function val=test_ComputeBoundXi(ProbabiliesPartitionsSortedDescendent,ValidLabelFcnValuesSortedDescendent_phsi_i1,M,n,X,varargin)
% test_ComputeBoundXi computes the values of the bounds stated by Theorem 1
% of the paper 'Bayes Labeling and Bayes Clustering Operators for Known 
% Random Labeled Point Processes'
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%
% See also test_ComputeBoundX2, test_ComputeBoundX1

suma=0;
for j=1:X
    IdxOfLabelFcnsHammingDistj=test_FindLabelFcnsHammingDistance(ValidLabelFcnValuesSortedDescendent_phsi_i1,ValidLabelFcnValuesSortedDescendent_phsi_i1(1,:),j,varargin);
    suma=suma+(M*n+X-2*j+1)*sum(ProbabiliesPartitionsSortedDescendent(IdxOfLabelFcnsHammingDistj));
    clear IdxOfLabelFcnsHammingDistj;
end
val=(M*n-suma)/(M*n+X+1);
return