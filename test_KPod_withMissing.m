function PredictedLabels=test_KPod_withMissing(GenericExpFolderName,ModelNumber)
% test_KPod clusters point sets using the K-POD algorithm
%
%

% Initial Values
NumClustersPerPartition=2;
max_iter=50;

PathName=[GenericExpFolderName '/SAMPLES/Sample_' ModelNumber];

% Loading feature vectors
FileName1='Sample_Data.def';
ModelNameDataSamples=fullfile(PathName,FileName1);
Vectors=load(ModelNameDataSamples,'Vectors');
NumSamples=size(Vectors,1);

FileName_='Sample_imputed.def';
ModelNameImputed=fullfile(PathName,FileName_);
Imputed=load(ModelNameImputed,'Imputed');

FileName_='Sample_missing.def';
ModelNameImputed=fullfile(PathName,FileName_);
Missing=load(ModelNameImputed,'Missing');

Latest=Vectors.*Missing;
mu = (sum(Vectors.*Missing)./sum(Missing))';
for i=1:size(Vectors,2)
    Latest(Missing(:,i)==0,i) = mu(i);
end

% Aplying KPod Means Clustering
[PredictedLabels,Centers] = kmeans(Latest, NumClustersPerPartition);

for j=2:max_iter
   A = zeros(NumSamples,NumClustersPerPartition);
   for ii=1:NumClustersPerPartition
       A(PredictedLabels==ii,ii)=1;
   end
   Latest = (Vectors.*Missing) + ((A*Centers).*(1-Missing));
   [PredictedLabels,Centers] = kmeans(Latest, NumClustersPerPartition); 
end

% Saving Predicted lables
FileName2='Predicted_Labels_KP.def';
ModelPredictedLabels=fullfile(PathName,FileName2);
save(ModelPredictedLabels,'PredictedLabels','-ascii');
return