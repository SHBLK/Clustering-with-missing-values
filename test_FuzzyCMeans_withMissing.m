function PredictedLabels=test_FuzzyCMeans_withMissing(GenericExpFolderName,ModelNumber)
% test_FuzzyCMeans clusters point sets using the Fuzzy C-means algorithm
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



% Aplying Fuzzy C Means Clustering
options= [NaN,NaN,NaN,0];
[center,U] = fcm(Imputed, NumClustersPerPartition, options);
maxU = max(U);
PredictedLabels=double(U(1,:) == maxU)'+1;

% Saving Predicted lables
FileName2='Predicted_Labels_FCM.def';
ModelPredictedLabels=fullfile(PathName,FileName2);
save(ModelPredictedLabels,'PredictedLabels','-ascii');
return