function PredictedLabels=test_RandomClustererV2LargeSamples(GenericExpFolderName,ModelNumber,n)
% test_RandomClustererV2LargeSamples clusters large point sets randomly
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%
% See also test_RandomClustererV2

% Initial Values
NumClustersPerPartition=2;

PathName=[GenericExpFolderName '/SAMPLES/Sample_' ModelNumber];

% Choosing a label function randomly
if n<=53
    RandIdx=randi((2^n)-1,1);
    BinRandomNumber=dec2bin(RandIdx,n);
    PredictedLabels=double(BinRandomNumber==49);
else
    NumParts=ceil(n/53);
    PredictedLabels=[];
    for i=1:NumParts
        if i==NumParts
            RandIdx=randi((2^(n-53*(i-1)))-1,1);
            BinRandomNumber=dec2bin(RandIdx,(n-53*(i-1)));
        else
            RandIdx=randi((2^53)-1,1);
            BinRandomNumber=dec2bin(RandIdx,53);
        end
        PredictedLabels=[PredictedLabels double(BinRandomNumber==49)];
    end
end
if PredictedLabels(1,1)==1
    PredictedLabels=1-PredictedLabels;
end
PredictedLabels=PredictedLabels'+1;

% Saving Predicted lables
FileName2='Predicted_Labels_RandomOp.def';
ModelPredictedLabels=fullfile(PathName,FileName2);
save(ModelPredictedLabels,'PredictedLabels','-ascii');
return