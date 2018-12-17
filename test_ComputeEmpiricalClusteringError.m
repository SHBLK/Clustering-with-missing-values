function Error=test_ComputeEmpiricalClusteringError(GenericExpFolderName,ModelNumber,Mask)
% test_ComputeEmpiricalClusteringError computes the empirical error of
% clustering a point set by comparying the labels assigned by the tested
% algorithm with the true ones obtained from the random labeled point
% process from which the point set tested come
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%

PathName=[GenericExpFolderName '/SAMPLES/Sample_' ModelNumber];

% Loading actual labels
FileName1='Sample_Labels.def';
ModelNameLabelsSamples=fullfile(PathName,FileName1);
ActualLabels=load(ModelNameLabelsSamples,'Labels');

% Loading predicted labels
FileName2=[Mask '.def'];
ModelOptimalLabels=fullfile(PathName,FileName2);
PredictedLabels=load(ModelOptimalLabels,'PredictedLabels');

Error=100*(1/length(ActualLabels))*min(sum(ActualLabels~=PredictedLabels),sum(ActualLabels~=(3-PredictedLabels)));
return