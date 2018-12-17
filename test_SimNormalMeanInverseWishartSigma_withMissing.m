function test_SimNormalMeanInverseWishartSigma_withMissing(m,v,k,S,NumberSamples,...
                                               GenericExpFolderName,MissingProb,...
                                               ModelNumber)

%
% test_SimNormalMeanFixedKnownSigma    generates feature vectors, labels,
%                                      missing pattern, and imputed data
%                                      from a 2 Gaussians Model with normal
%                                      means and known covariances
%
%  test_SimNormalMeanFixedKnownSigma(m,v,Sigma,NumberSamples,...
%                                    GenericExpFolderName,...
%                                    ModelNumber)
%
% INPUTS:
%
%      m, v, k, S             - Gaussian inverse-Wishart distribution
%                             parameters
%      NumberSamples          - 1xL: Each column contains the number of
%                             points for class Li.
%      GenericExpFolderName   - String vector. GenericExpFolderName
%                             contains the pathname to save the data used
%                             to generate the dataset. If the provided
%                             pathname already exists,
%                             test_SimNormalMeanFixedKnownSigma keeps the
%                             existing files.
%      ModelNumber            - String vector: ModelNumber contains the tag
%                             used to identify the different samples drawn
%                             with test_SimNormalMeanFixedKnownSigma.
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
% Defining Hyperparameters and file names

% Generating Random Covariance Matrices drawn from a Inverse Wishart Dist. 
RandSigma=test_InverseWishartCovariance(k,S);
Sigma1(:,:,1)=RandSigma(:,:,1)/v(1,1);
Sigma1(:,:,2)=RandSigma(:,:,2)/v(1,2);
RandMU=test_GaussianMean(m,Sigma1);

model='gau gau';
type='mixture';

% Saving Hyperparameters used for the simulation
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

% Saving Parameters used to draw samples
FileName2=['Sample_Parameters_' ModelNumber '.def'];
ModelNameParam=fullfile(PathName,FileName2);
data2=struct('ModelInformation',['Model with ' num2str(size(m,2)) '-dimensional space'],...
            'classes',size(m,1),'dim',size(m,2),...
            'folder',ModelNameHypParam(1,1:end-4),...
            'templates',RandMU,...
            'variances',[RandSigma(:,:,1);RandSigma(:,:,2)],...
            'quantities',NumberSamples,...
            'model',model,'type',type,'missingRate',num2str(MissingProb)); 
bb_cl_saveparam(ModelNameParam,data2);

% Generating a random sample (feature vectors and labels)
[Vectors Labels]=test_GaussianSamples(RandMU,RandSigma,NumberSamples);
Labels=Labels';

Missing = test_GenerateMissingPattern(sum(NumberSamples),size(RandMU,2),MissingProb);

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

FileName6='Sample_imputed.def';
ModelNameImputed=fullfile(PathName,FileName6);
save(ModelNameImputed,'Imputed','-ascii');
return