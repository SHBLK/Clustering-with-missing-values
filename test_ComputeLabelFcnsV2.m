function [data0,data1]=test_ComputeLabelFcnsV2(n1,n2,GenericExpFolderName,ProjectName)
% Developed by Marco Benalcázar - September 27/2012
%
% test_ComputeLabelFcnsV2  Computes all the possible label functions in 
%                          order to split a given set of points
%                          S={x1,...,xn} into two clusters of n1 and n2
%                          points each
%
%  test_ComputeLabelFcnsV2(n1,n2,GenericExpFolderName,ProjectName)
%
% INPUTS:
%
%      GenericExpFolderName   - String vector. GenericExpFolderName
%                             contains the pathname to save the data used
%                             to generate the dataset. If the provided
%                             pathname already exists,
%                             test_SimNormalMeanFixedKnownSigma keeps the
%                             existing files.
%      ModelNumber            - String vector: ModelNumber contains the tag
%                             used to identify the different samples drawn
%                             with test_SimNormalMeanFixedKnownSigma.
%    
%
% OUTPUTS:
%
%  LabelFunctions     - .def file with the valid label functions for
%                        the clusters of each possible partition.
%                        In the .def file, the first row of phsi_i1 
%                        is equal to 1 minus the first row of phsi_i2,
%                        and so on. Each column represent the label for the
%                        each point from S (e.g., x1,...,xn).
%                        
%
% NOTE: This function uses some functions form the Toolbox bbcl developed by
% Marcel Brun, so make sure you have that Toolbox before run this function.
Tstart=tic;
n=n1+n2;

% Computing a Truth table with 2^n combinations
Table = test_TruthTable(n);

% Saving All Possible Label Functions
[Labels1 Labels2] = test_SplitLabelFcns(Table);
data0=struct('FileInformation','All Possible Label Functions',...
            'phsi',Labels1);

% Finding the Valid Label Functions
ValidLabels = test_FindValidLabels(Table,n1,n2);
[Cluster1 Cluster2] = test_SplitLabelFcns(ValidLabels);

% Saving Valid Label functions
data1=struct('FileInformation','Valid Label Functions',...
            'phsi_i1',Cluster1,...
            'phsi_i2',Cluster2);
fprintf('Computing label functions. Elapsed time: %d s\r\r',toc(Tstart));
return