function NumPartitions=test_ComputeNumberPartitions(n1,n2,NumClustersPerPartition)
% test_ComputeNumberPartitions finds the number of partitions that assign
% the correct number of points for each cluster
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%
n=n1+n2;
for i=NumClustersPerPartition
    if n1==n2
        NumLabelFcns=factorial(n)/(factorial(n-n1)*factorial(n1));
    else
        NumCombinations0s=factorial(n)/(factorial(n-n1)*factorial(n1));
        NumCombinations1s=factorial(n)/(factorial(n-n2)*factorial(n2));
        NumLabelFcns=NumCombinations0s+NumCombinations1s;
    end
end
NumPartitions=NumLabelFcns/2;
return