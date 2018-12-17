function [Indicator,ErrorPartitions,PredictedLabels,HammingDistance_PredictedPmaxPartitions]=test_ComputeErrorClustererV2(AllLabelFcnValues,ValidLabelFcnValues,n1,n2,Prphi1_S,Prphi2_S,varargin)
% test_ComputeErrorClustererV2 finds the optimal partition for a given point
% set applying exhaustive search (Bayes partition)
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%

Tstart1=tic;
% Computing and sorting Probabilities for possible (valid) partitions
ProbabiliesPartitions=Prphi1_S+Prphi2_S;
[ProbabiliesPartitionsSortedDescendent,IDX]=sort(ProbabiliesPartitions,'descend');
fprintf('Sorting probabilities. Elapsed time: %d s\r',toc(Tstart1));

if nargin==6
    varargin{1}='NoFullMatrixW';
end

NumPartitions=test_ComputeNumberPartitions(8,8,2);
if strcmp(varargin{1},'FullMatrixW') && size(ProbabiliesPartitionsSortedDescendent,1)>=NumPartitions
    varargin{1}='NoFullMatrixW';
    fprintf('\rWARNING: Method changed from "FullMatrixW" to "NoFullMatrixW" because matrices are out of the RAM memory size\r');
end
Method=varargin{1};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: In the current code, Columns of matrix W are all the candidate    %
% partitions, whereas rows are all the possible partitions                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Tstart2=tic;
%Sequential search of the best vector of labels
n=n1+n2;
if mod(n,2)==0
    M=0.5;
else
    M=(n-1)/(2*n);
end

if ProbabiliesPartitionsSortedDescendent(1) > (M*n/(M*n+1))
    Indicator=0;
    % Computing Matrix W for the most likely label function
    Rows=(1:size(ValidLabelFcnValues.phsi_i1,1))';
    Columns=1;
    W_likely=test_ComputeMatrixWv2(Rows,Columns,ValidLabelFcnValues.phsi_i1(IDX,:),ValidLabelFcnValues.phsi_i1(IDX,:),ValidLabelFcnValues.phsi_i2(IDX,:));
    
    % Error of the betst partition
    ErrorPartitions = W_likely'*ProbabiliesPartitionsSortedDescendent;
    
    % Finding the optimal vector of labels
    PredictedLabels=ValidLabelFcnValues.phsi_i1(IDX(1),:)'+1;
    fprintf('The Best Label Function --> Highest Probability. Elapsed time %d s\r',toc(Tstart2));
elseif ProbabiliesPartitionsSortedDescendent(1) > test_ComputeBoundX1(ProbabiliesPartitionsSortedDescendent,ValidLabelFcnValues.phsi_i1(IDX,:),M,n)
    Indicator=1;
    [ValidLabelFcnValuesSortedDescendent_phsi_i1,...
        ValidLabelFcnValuesSortedDescendent_phsi_i2,...
        LabelFcnValues,...
        ProbabiliesPartitionsSortedDescendent]=test_SortProbabilities(IDX,ValidLabelFcnValues,AllLabelFcnValues,ProbabiliesPartitionsSortedDescendent,n1,n2);
    % Finding vector of labels that have a hamming distance of 1 from the
    % vector of labels that corresponds to tha highest probability
    HammingDistance=1;
    Columns=test_FindLabelFcnsHammingDistance(LabelFcnValues,ValidLabelFcnValuesSortedDescendent_phsi_i1(1,:),HammingDistance,'LesserOrEqual');
    [ErrorPartitions,IdxOfBestLabelFcn]=test_FindOptimalLabelFcn(ProbabiliesPartitionsSortedDescendent,LabelFcnValues(Columns,:),ValidLabelFcnValuesSortedDescendent_phsi_i1,ValidLabelFcnValuesSortedDescendent_phsi_i2,M,Method);
    IdxOfBestPartition=Columns(IdxOfBestLabelFcn);
    PredictedLabels=LabelFcnValues(IdxOfBestPartition,:)'+1;
    fprintf('The Best Label Function --> Set of labels with at most 1 mismatch. Elapsed time %d s\r',toc(Tstart2));
elseif ProbabiliesPartitionsSortedDescendent(1) > test_ComputeBoundX2(ProbabiliesPartitionsSortedDescendent,ValidLabelFcnValues.phsi_i1(IDX,:),M,n)
    Indicator=2;
    [ValidLabelFcnValuesSortedDescendent_phsi_i1,...
        ValidLabelFcnValuesSortedDescendent_phsi_i2,...
        LabelFcnValues,...
        ProbabiliesPartitionsSortedDescendent]=test_SortProbabilities(IDX,ValidLabelFcnValues,AllLabelFcnValues,ProbabiliesPartitionsSortedDescendent,n1,n2);
    % Finding vector of labels that have a hamming distance of 2 from the
    % vector of labels that corresponds to tha highest probability
    HammingDistance=2;
    Columns=test_FindLabelFcnsHammingDistance(LabelFcnValues,ValidLabelFcnValuesSortedDescendent_phsi_i1(1,:),HammingDistance,'LesserOrEqual');
    [ErrorPartitions,IdxOfBestLabelFcn]=test_FindOptimalLabelFcn(ProbabiliesPartitionsSortedDescendent,LabelFcnValues(Columns,:),ValidLabelFcnValuesSortedDescendent_phsi_i1,ValidLabelFcnValuesSortedDescendent_phsi_i2,M,Method);
    IdxOfBestPartition=Columns(IdxOfBestLabelFcn);
    PredictedLabels=LabelFcnValues(IdxOfBestPartition,:)'+1;
    fprintf('The Best Label Function --> Set of labels with at most 2 mismatches. Elapsed time %d s\r',toc(Tstart2));
else
    Indicator=3;
    while Indicator <= M*n
        if ProbabiliesPartitionsSortedDescendent(1) > test_ComputeBoundXi(ProbabiliesPartitionsSortedDescendent,ValidLabelFcnValues.phsi_i1(IDX,:),M,n,Indicator) && Indicator < M*n
            [ValidLabelFcnValuesSortedDescendent_phsi_i1,...
                ValidLabelFcnValuesSortedDescendent_phsi_i2,...
                LabelFcnValues,...
                ProbabiliesPartitionsSortedDescendent]=test_SortProbabilities(IDX,ValidLabelFcnValues,AllLabelFcnValues,ProbabiliesPartitionsSortedDescendent,n1,n2);
            % Finding vector of labels that have a hamming distance of 2 from the
            % vector of labels that corresponds to tha highest probability
            HammingDistance=Indicator;
            Columns=test_FindLabelFcnsHammingDistance(LabelFcnValues,ValidLabelFcnValuesSortedDescendent_phsi_i1(1,:),HammingDistance,'LesserOrEqual');
            [ErrorPartitions,IdxOfBestLabelFcn]=test_FindOptimalLabelFcn(ProbabiliesPartitionsSortedDescendent,LabelFcnValues(Columns,:),ValidLabelFcnValuesSortedDescendent_phsi_i1,ValidLabelFcnValuesSortedDescendent_phsi_i2,M,Method);
            IdxOfBestPartition=Columns(IdxOfBestLabelFcn);
            PredictedLabels=LabelFcnValues(IdxOfBestPartition,:)'+1;
            fprintf('The Best Label Function --> Set of labels with at most %d mismatches. Elapsed time %d s\r',Indicator,toc(Tstart2));
            break;
        elseif Indicator == M*n
            [ValidLabelFcnValuesSortedDescendent_phsi_i1,...
                ValidLabelFcnValuesSortedDescendent_phsi_i2,...
                LabelFcnValues,...
                ProbabiliesPartitionsSortedDescendent]=test_SortProbabilities(IDX,ValidLabelFcnValues,AllLabelFcnValues,ProbabiliesPartitionsSortedDescendent,n1,n2);
            [ErrorPartitions,IdxOfBestLabelFcn]=test_FindOptimalLabelFcn(ProbabiliesPartitionsSortedDescendent,LabelFcnValues,ValidLabelFcnValuesSortedDescendent_phsi_i1,ValidLabelFcnValuesSortedDescendent_phsi_i2,M,Method);
            PredictedLabels=LabelFcnValues(IdxOfBestLabelFcn,:)'+1;
            fprintf('The Best Label Function --> Set of labels with at most %d mismatches (Full matrix W). Elapsed time %d s\r',Indicator,toc(Tstart2));
            break;
        else
            Indicator=Indicator+1;
        end
    end
end
LabelFcnsIn=ValidLabelFcnValues.phsi_i1(IDX(1),:);
LabelFcnTargetMatrix=PredictedLabels'-1;
HammingDistance_PredictedPmaxPartitions=min(sum(LabelFcnsIn~=LabelFcnTargetMatrix,2),sum(LabelFcnsIn~=(1-LabelFcnTargetMatrix),2));
return