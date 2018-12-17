function [ErrorPartitions,BestIndices]=test_FindOptimalLabelFcn(ProbabiliesPartitionsSortedDescendent,LabelFcns,ValidLabelFcnsSortedDescendent_phsi_i1,ValidLabelFcnsSortedDescendent_phsi_i2,M,Method)
% test_FindOptimalLabelFcn finds the optimal partition for a given point
% set applying exhaustive search (Bayes partition)
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%

step=50;
count=0;
if size(ProbabiliesPartitionsSortedDescendent,1)<=step || strcmp(Method,'FullMatrixW')
    Columns=(1:size(LabelFcns,1))';
    Rows=(1:size(ProbabiliesPartitionsSortedDescendent,1))';
    W_likely=test_ComputeMatrixWv2(Rows,Columns,LabelFcns,ValidLabelFcnsSortedDescendent_phsi_i1,ValidLabelFcnsSortedDescendent_phsi_i2);
    ErrorPartitions=W_likely'*ProbabiliesPartitionsSortedDescendent;
    [~,BestIndices]=min(ErrorPartitions);
else
    likely_index = 0;
    index = 0;
    likely_probability = 0;
    surviving_index = [0 0];
    i = 0;
    lower_bound = 0;
    while length(surviving_index)~=1
        index = index+step;
        if size(ProbabiliesPartitionsSortedDescendent,1)<index(end)
            try
                likely_index = likely_index(end)+(1:(size(ProbabiliesPartitionsSortedDescendent,1)-likely_index(end)));
            catch
                [ErrorPartitions,surviving_index]=min(lower_bound);
                BestIndices=BestIndices(surviving_index);
            end
        else
            likely_index = likely_index(end) + (1:step);
        end
        likely_probability = likely_probability + sum(ProbabiliesPartitionsSortedDescendent(likely_index));
        constant = M*(1 - likely_probability);
        P_likely = ProbabiliesPartitionsSortedDescendent(likely_index);
        if i == 0
            Rows=likely_index';
            Columns=(1:size(LabelFcns,1))';
            W_likely=test_ComputeMatrixWv2(Rows,Columns,LabelFcns,ValidLabelFcnsSortedDescendent_phsi_i1,ValidLabelFcnsSortedDescendent_phsi_i2);
            lower_bound = W_likely'*P_likely;
            % Error Partitions has the error values of those likely partions
            ErrorPartitions=lower_bound;
            BestIndices=(1:size(LabelFcns,1))';
        else
            Columns=BestIndices;
            Rows=likely_index;
            W_likely=test_ComputeMatrixWv2(Rows,Columns,LabelFcns,ValidLabelFcnsSortedDescendent_phsi_i1,ValidLabelFcnsSortedDescendent_phsi_i2);
            lower_bound = lower_bound(surviving_index) + W_likely'*P_likely;
        end
        surviving_index = find(lower_bound < (min(lower_bound) + constant));
        if isempty(surviving_index)
            [ErrorPartitions,surviving_index]=min(lower_bound);
            BestIndices=BestIndices(surviving_index);
            break;
        end
        OldSizeBestIndices=length(BestIndices);
        BestIndices=BestIndices(surviving_index);
        NewSizeBestIndices=length(BestIndices);
        i = i + 1;
    end
end
return