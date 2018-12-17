function PredictedLabels=test_Hierarchical_withMissing(GenericExpFolderName,ModelNumber,TypeDist,Method)
% test_Hierarchical clusters point sets using the Hierarchical
% algorithm (implements the single and complete linkage of this algorithm)
%

%

% Initial Values
NumClustersPerPartition=2;

PathName=[GenericExpFolderName '/SAMPLES/Sample_' ModelNumber];

% Loading feature vectors
FileName1='Sample_Data.def';
ModelNameDataSamples=fullfile(PathName,FileName1);
Vectors=load(ModelNameDataSamples,'Vectors');
NumSamples=size(Vectors,1);

FileName_='Sample_imputed.def';
ModelNameImputed=fullfile(PathName,FileName_);
Imputed=load(ModelNameImputed,'Imputed');

% Aplying Hierarchical Clustering
[PredictedLabels,CEN,LINK] = bb_cl_hierarchical(Imputed,NumClustersPerPartition,TypeDist,Method);
PredictedLabels=PredictedLabels';

% Saving Predicted lables
FileName2=['Predicted_Labels_Hierarchical_' Method '.def'];
ModelPredictedLabels=fullfile(PathName,FileName2);
save(ModelPredictedLabels,'PredictedLabels','-ascii');
return