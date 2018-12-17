function [LabelFcn1 LabelFcn2]=test_SplitLabelFcns(ValidLabels)
% Developed by Marco Benalcázar - September 27/2012
%
% test_SplitLabelFcns   Finds the two label functions that produce the ith
%                       partition of the dataset S = {x1,...,xn}, of N 
%                       points.
%  
%
% [LabelFcn1 LabelFcn2] = test_SplitLabelFcns(ValidLabels)
%
% INPUTS:
%
%   ValidLabels        - PxN. P is the number of all valid possible label
%                       functions, where each valid label function partitions 
%                       the dataset S, of N points, in two clusters. 
%                       ValidLabels contains all valid possible label 
%                       functions.
% OUTPUTS:
%
%   LabelFcn1           - (P/2)xN. Each row of LabelFcn1 contains the vector 
%                       of labels for each possible partition of S.
%   LabelFcn2           - (P/2)xN. Each row of LabelFcn2 is equal to  
%                       1 - LabelFcn2.
%
%
% SEE ALSO: test_TruthTable, test_ValidLabels, test_ComputeW.
%
[m n]=size(ValidLabels);
LabelFcn1=ValidLabels(1:m/2,:);
LabelFcn2=1-LabelFcn1;
return