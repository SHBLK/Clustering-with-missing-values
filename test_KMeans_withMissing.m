function PredictedLabels=test_KMeans_withMissing(GenericExpFolderName,ModelNumber)
% test_KMeans clusters point sets using the K-means algorithm
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

% Aplying K Means Clustering
[PredictedLabels,Centers] = kmeans(Imputed, NumClustersPerPartition);

% Saving Predicted lables
FileName2='Predicted_Labels_KM.def';
ModelPredictedLabels=fullfile(PathName,FileName2);
save(ModelPredictedLabels,'PredictedLabels','-ascii');
return