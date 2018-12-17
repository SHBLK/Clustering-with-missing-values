function PredictedLabels=test_FuzzyCMeansOPS_withMissing(GenericExpFolderName,ModelNumber)
% test_FuzzyCMeansOPS clusters point sets using the Fuzzy C-means with Optimal Completion Strategy algorithm
%
%

% Initial Values
NumClustersPerPartition=2;
max_iter=100;
thresh_=1e-5;
m_=2;

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

p_=size(Vectors,2);

Complete=Vectors(sum(Missing,2)==p_,:);
options= [NaN,NaN,NaN,0];
[center,U] = fcm(Complete, NumClustersPerPartition, options);

D_=zeros(NumClustersPerPartition,NumSamples);

for i=1:NumClustersPerPartition
   D_(i,:) = ((p_./sum(Missing,2)).*(sum(((Vectors-repmat(center(i,:),NumSamples,1)).^2).*Missing,2) ) )'; 
end


[~,Label_temp] = min(D_);
A = zeros(NumSamples,NumClustersPerPartition);
for ii=1:NumClustersPerPartition
    A(Label_temp'==ii,ii)=1;
end
Latest = (Vectors.*Missing) + ((A*center).*(1-Missing));
% Aplying Fuzzy C Means Clustering
DD=zeros(NumClustersPerPartition,NumSamples);
center_prev=center;
for j=2:max_iter

    for i=1:NumClustersPerPartition
        DD(i,:) = (sum(((Latest-repmat(center(i,:),NumSamples,1)).^2),2) )';
    end
    U = (DD.^(1/(1-m_)));
    U = U ./ repmat(sum(U,1),NumClustersPerPartition,1);
    center= ( (U.^m_)*Latest )./( repmat(sum(U.^m_,2),1,p_)  );
    
    if max(max(abs(center-center_prev)))<thresh_
       break; 
    end
    center_prev=center;

    temp_impute = ( (U'.^m_)*center )./( repmat(sum(U'.^m_,2),1,p_)  );
    Latest = (Vectors.*Missing) + (temp_impute.*(1-Missing));
end

maxU = max(U);
PredictedLabels=double(U(1,:) == maxU)'+1;

% Saving Predicted lables
FileName2='Predicted_Labels_FCM_OPS.def';
ModelPredictedLabels=fullfile(PathName,FileName2);
save(ModelPredictedLabels,'PredictedLabels','-ascii');
return