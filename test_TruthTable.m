function Table = test_TruthTable(n)
% Developed by Marco Benalcázar - September 27/2012
%
% test_TruthTable       finds all the possible label functions for a given 
%                       dataset S = {x1,...,xn}, of N points, to be partitioned in 
%                       two clusters with equal number of points for each
%                       cluster.
%  
%
% Table = test_TruthTable(n)
%
% INPUTS:
%
%      n              - Scalar number: n is the number of points of the
%                       dataset S to be partitioned in two clusters, with
%                       equal number of points for each cluster
%
% OUTPUTS:
%
%   Table           - MxN. M = 2^N is the number of all possible label
%                       functions for a dataset of N points. Table
%                       contains all the possible label functions.
%
%
% SEE ALSO: test_ValidLabels, test_LabelFcn, test_ComputeW.
%
Table=fliplr(rem(floor([0:((2^n)-1)].'* pow2(0:-1:-n+1)),2));
return
