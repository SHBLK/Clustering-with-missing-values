function TableOut = test_FindValidLabels(TableIn,n1,n2)
% Developed by Marco Benalcázar - September 27/2012
%
% test_FindValidLabels      finds all the valid label functions between all the
%                       possible label functions for a given dataset 
%                       S = {x1,...,xn}, of N points, to be partitioned in 
%                       two clusters with equal number of points for each
%                       cluster.
%  
%
% TableOut = test_FindValidLabels(TableIn)
%
% INPUTS:
%       
%   TableIn           - MxN. M = 2^N is the number of all possible label
%                       functions for a dataset of N points. TableIn
%                       contains all the possible label functions.
%
% OUTPUTS:
%
%   TableOut           - PxN. P is the number of all valid possible label
%                       functions, where each valid label function partitions 
%                       the dataset S, of N points, in two clusters, where
%                       each cluster has N/2 points. TableOut contains all 
%                       valid possible label functions.
%
%
% SEE ALSO: test_TruthTable, test_LabelFcn, test_ComputeW.

Numb1sPerRow=sum(TableIn,2);
TableOut=TableIn(((Numb1sPerRow==n1)|(Numb1sPerRow==n2)),:);
return