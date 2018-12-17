function TheoreticalError=test_ComputeTheoreticalClusteringErrorValue(Prphi1_S,Prphi2_S,PredictedLabels,ValidLabelFcnValues,M)
% test_ComputeTheoreticalClusteringErrorValue computes the exact error value for
% the partition 'PredictedLabels' of the set S using the probabilities 
% 'Prphi1_S','Prphi2_S' computed using the model from which the set S was drawn.
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%
Columns=size(PredictedLabels,1);
Rows=(1:size(ValidLabelFcnValues.phsi_i1,1))';
for i=1:Columns
    if PredictedLabels(i,1)==1
        PredictedLabels(i,:)=1-PredictedLabels(i,:);
    end
end

% Computing mismatch between the predicted partition and all the possible
% partitions
W=test_ComputeMatrixWv2(Rows,Columns,PredictedLabels,ValidLabelFcnValues.phsi_i1,ValidLabelFcnValues.phsi_i2);
TheoreticalError=100*(W'*(Prphi1_S+Prphi2_S));
return