function test_RealDataSampling_withMissing(real_data0,real_data1,calibrationSize,NumberSamples,...
                                               GenericExpFolderName,MissingProb,...
                                               ModelNumber)

% test_RealDataSampling_withMissing    generates feature vectors and labels
%                                      and constructs data-based priors for
%                                      a 2 component real data
%
%  test_RealDataSampling_withMissing(real_data0,real_data1,calibrationSize,NumberSamples,...
%                                               GenericExpFolderName,MissingProb,...
%                                               ModelNumber)
%
% INPUTS:
%
%      real_data0             - Real data for component 1 (each row is an observation).
%      real_data1             - Real data for component 2.
%      calibrationSize        - Number of features (last columns) to discard and use for prior
%                               calibration
%      NumberSamples          - 1xL: Each column contains the number of
%                             points for class Li.
%      GenericExpFolderName   - String vector. GenericExpFolderName
%                             contains the pathname to save the data used
%                             to generate the dataset. If the provided
%                             pathname already exists,
%                             test_RealDataSampling_withMissing keeps the
%                             existing files.
%      ModelNumber            - String vector: ModelNumber contains the tag
%                             used to identify the different samples drawn
%                             with test_RealDataSampling_withMissing.
%      
%
% OUTPUTS:
%
%      - ModelHyperparameters: .def file with the hyperparameters
%      
%      - Data_Sample:          .def file with the drawn feature vectors
%      - Labels_Sample:        .def file with the drawn true labels
%      - Missing pattern:      .def file containing the missing pattern
%      - Imputed data:         .def file containing feature vectors after
%                              imputation
%      
%
% 
%

model='gau gau';
type='mixture';

[Vectors_complete,Labels] = test_RealSamples(real_data0,real_data1,NumberSamples);
[Vectors,S,m,v,k] = test_ConstructPriors(Vectors_complete,calibrationSize);

% Saving Hyperparameters
PathName=[GenericExpFolderName '/SAMPLES/Sample_' ModelNumber];
FileName1='Model_Hyperparameters.def';
ModelNameHypParam=fullfile(PathName,FileName1);
mkdir(PathName);
data1=struct('ModelInformation',['Model with ' num2str(size(m,2)) '-dimensional space'],...
            'classes',size(m,1),'dim',size(m,2),...
            'folder',ModelNameHypParam(1,1:end-4),...
            'templates',m,...
            'S',[S(:,:,1);S(:,:,2)],...
            'quantities',NumberSamples,...
            'v',v,...
            'k',k,...
            'model',model,'type',type,'missingRate',num2str(MissingProb)); 
bb_cl_saveparam(ModelNameHypParam,data1);            


Labels=Labels';

Missing = test_GenerateMissingPattern(sum(NumberSamples),size(m,2),MissingProb);

Imputed=Gibbs_imputation(Vectors,Missing);

% Saving feature vectors in a .txt file
FileName3='Sample_Data.def';
ModelNameDataSamples=fullfile(PathName,FileName3);
save(ModelNameDataSamples,'Vectors','-ascii');

% Saving true labels in a .txt file
FileName4='Sample_Labels.def';
ModelNameLabelsSamples=fullfile(PathName,FileName4);
save(ModelNameLabelsSamples,'Labels','-ascii');

% Saving missing pattern in a .txt file
FileName5='Sample_missing.def';
ModelNameMissing=fullfile(PathName,FileName5);
save(ModelNameMissing,'Missing','-ascii');

% Saving imputed data in a .txt file
FileName6='Sample_imputed.def';
ModelNameImputed=fullfile(PathName,FileName6);
save(ModelNameImputed,'Imputed','-ascii');
return