function TheoreticalError=test_ComputeTheoreticalClusteringErrorLargeSamples(PredictedLabels,M,MaxHammingDistance,GenericExpFolderName,ModelNumber)
% test_ComputeTheoreticalClusteringErrorLargeSamples computes the exact error for
% the partition 'PredictedLabels' of a point set S using the probabilities 
% 'Prphi1_S', 'Prphi2_S' computed using the model from which the set was drawn.
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

% Loading all possible label function values
FileName6=['PossibleLabelFcnValues_' num2str(MaxHammingDistance) '.mat'];
ListOfPossibleLabelFunctionValues=fullfile(PathName,FileName6);
load(ListOfPossibleLabelFunctionValues,'CandidateLabelFcns');

% Filenames for files that contain probabilities for the optimal operator
FileName1=['Probabilities1_LabelFcns_SubOptimalPseed_HamDist_' num2str(MaxHammingDistance) '.def'];
ModelProbabilities1=fullfile(PathName,FileName1);
FileName2=['Probabilities2_LabelFcns_SubOptimalPseed_HamDist_' num2str(MaxHammingDistance) '.def'];
ModelProbabilities2=fullfile(PathName,FileName2);
% Saving Probabilities
Prphi1_S=load(ModelProbabilities1,'Pr_phi1');
Prphi2_S=load(ModelProbabilities2,'Pr_phi2');

TheoreticalError=test_ComputeTheoreticalClusteringErrorValue(Prphi1_S,Prphi2_S,PredictedLabels,CandidateLabelFcns,M);
return