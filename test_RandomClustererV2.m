function PredictedLabels=test_RandomClustererV2(GenericExpFolderName,ModelNumber,n)
% test_RandomClustererV2 clusters a small point set randomly
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%
% See also test_RandomClustererV2LargeSamples

% Initial Values
NumClustersPerPartition=2;

PathName=[GenericExpFolderName '/SAMPLES/Sample_' ModelNumber];

% Choosing a label function randomly
RandIdx=randi((2^n)-1,1);
BinRandomNumber=dec2bin(RandIdx,n);
PredictedLabels=double(BinRandomNumber==49);
if PredictedLabels(1,1)==1
    PredictedLabels=1-PredictedLabels;
end
PredictedLabels=PredictedLabels'+1;

% Saving Predicted lables
FileName2='Predicted_Labels_RandomOp.def';
ModelPredictedLabels=fullfile(PathName,FileName2);
save(ModelPredictedLabels,'PredictedLabels','-ascii');
return