function TheoreticalError=test_ComputeTheoreticalClusteringError(PredictedLabels,ValidLabelFcnValues,M,GenericExpFolderName,ModelNumber)
% test_ComputeTheoreticalClusteringError computes the exact error for
% the partition 'PredictedLabels' of a point set S using the probabilities 
% 'Prphi1_S','Prphi2_S' computed using the model from which the set was drawn.
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%
% See also test_ComputeTheoreticalClusteringErrorValue
%

% Loading probabilities
PathName=[GenericExpFolderName '/SAMPLES/Sample_' ModelNumber];
% Filenames for files that contain probabilities for the optimal operator
FileName1='Probabilities1_LabelFcns_Optimal.def';
ModelProbabilities1=fullfile(PathName,FileName1);
FileName2='Probabilities2_LabelFcns_Optimal.def';
ModelProbabilities2=fullfile(PathName,FileName2);
% Saving Probabilities
Prphi1_S=load(ModelProbabilities1,'Prphi1_S');
Prphi2_S=load(ModelProbabilities2,'Prphi2_S');

TheoreticalError=test_ComputeTheoreticalClusteringErrorValue(Prphi1_S,Prphi2_S,PredictedLabels,ValidLabelFcnValues,M);
return