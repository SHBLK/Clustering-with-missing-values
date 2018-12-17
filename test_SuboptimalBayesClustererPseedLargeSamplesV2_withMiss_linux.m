function [EmpiricalErrorBayesSubOptPseedClusterer,...
    TheoreticalErrorBayesSubOptPseedClusterer,...
    TimeOfProcessing_BayesSubOptPseed,...
    TimeComputationPossibleLabelFcnValues,...
    TimeComputationProbabilitiesPossibleLabelFcnValues]=test_SuboptimalBayesClustererPseedLargeSamplesV2_withMiss_linux(GenericExpFolderName,ModelNumber,...
    TypeModel,Vectors,HypParam,n1,n2,MaxHammingDistance,MissMatr)

% test_SuboptimalBayesClustererPseedLargeSamplesV2 implements the suboptimal Pseed
% clustering algorithm for larger point sets (i.e. sets with more than 25 points)
%


%
% To process small point sets use test_SuboptimalBayesClustererPseed 
%
PathName=[GenericExpFolderName '/SAMPLES/Sample_' ModelNumber];

%Computing the partition seed for suboptimal algorithm
n=n1+n2;
if mod(n,2)==0
    M=0.5;
else
    M=(n-1)/(2*n);
end
if abs(n1-n2)==1
    HammingDistance=1;
else
    HammingDistance=2;
end
TStartFindingPseed=tic;
[PartitionSeed,~]=test_GeneratePartitionSeed_withMiss(TypeModel,Vectors,HypParam,n1,n2,HammingDistance,MissMatr,5); %Shahin
tFindingPseed=toc(TStartFindingPseed);

%Searching possible partitions
[ListOfPossiblePartitions,...
    ListOfHammingDistances,...
    ListTimeOfProcessing,...
    ListRAMUsage,...
    tRankingProbabilities]=test_FindListOfPossiblePartitions_withMiss(TypeModel,...
    Vectors,...
    HypParam,...
    PartitionSeed,...
    MaxHammingDistance,MissMatr);

EmpiricalErrorBayesSubOptPseedClusterer=[];
TheoreticalErrorBayesSubOptPseedClusterer=[];
CandidateLabelFcns.phsi_i1=[];
CandidateLabelFcns.phsi_i2=[];
Pr_phi1_NoNorm=[];
Pr_phi2_NoNorm=[];
tComputationProbabilities=0;
tComputationLabelFcnValues=0;
count=0;
% CandidatePartitions=[];
% TStartSpaceOfSearch=tic;
% MaxHammingDistanceCandidatePartitions=2;
% for ii=0:MaxHammingDistanceCandidatePartitions
%     CandidatePartitions=[CandidatePartitions;test_GenerateCandidatePartitionsWithGivenHammingDistance(n1,n2,PartitionSeed,ii)];
% end
% AllLabelFcnValues.phsi=CandidatePartitions;
% TFindingSpaceOfSearch=toc(TStartSpaceOfSearch);
for i=0:MaxHammingDistance
        
%     TStartProbabilities=tic;
%     %Generating possibles partitions with a given Hamming distance   
%     [ListOfPossiblePartitions.phsi_i1,PossibleHammingDistances]=test_GeneratePartitionsWithGivenHammingDistance(n1,n2,PartitionSeed,i);
%     
%     CandidateLabelFcns.phsi_i1=[CandidateLabelFcns.phsi_i1;ListOfPossiblePartitions.phsi_i1];
%     CandidateLabelFcns.phsi_i2=1-CandidateLabelFcns.phsi_i1;
%     
%     tComputationLabelFcnValues=tComputationLabelFcnValues+toc(TStartProbabilities);
%     
%     TStartProbabilities=tic;
%     % Computing probabilities of all possible partitions
%     [Prphi1_S,Prphi2_S]=test_ComputeProbabilityOfPartitions(TypeModel,Vectors,HypParam,ListOfPossiblePartitions,'NoNormalize');
%     
%     Pr_phi1_NoNorm=[Pr_phi1_NoNorm;Prphi1_S];
%     Pr_phi2_NoNorm=[Pr_phi2_NoNorm;Prphi2_S];
%     
%     %Normalizing Probabilities
%     Den=sum(Pr_phi1_NoNorm+Pr_phi2_NoNorm);
%     Pr_phi1=Pr_phi1_NoNorm/Den;
%     Pr_phi2=Pr_phi2_NoNorm/Den;
%     tComputationProbabilities=tComputationProbabilities+toc(TStartProbabilities);
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(find(ListOfHammingDistances==i))
    TStartProbabilities=tic;
    count=count+1;
    if i==0
        CandidateLabelFcns.phsi_i1=ListOfPossiblePartitions{1};
        CandidateLabelFcns.phsi_i2=1-CandidateLabelFcns.phsi_i1;
        TFindingSpaceOfSearch=ListTimeOfProcessing(1);
    else
        CandidateLabelFcns.phsi_i1=[CandidateLabelFcns.phsi_i1;ListOfPossiblePartitions{count}];
        CandidateLabelFcns.phsi_i2=1-CandidateLabelFcns.phsi_i1;
        TFindingSpaceOfSearch=sum(ListTimeOfProcessing(1:count));
    end
    TStartProbabilities=tic;
    % Computing probabilities of all possible partitions
    aux.phsi_i1=ListOfPossiblePartitions{count};
    [Prphi1_S,Prphi2_S]=test_ComputeProbabilityOfPartitions_withMiss(TypeModel,Vectors,HypParam,aux,MissMatr,'NoNormalize');
    Pr_phi1_NoNorm=[Pr_phi1_NoNorm;Prphi1_S];
    Pr_phi2_NoNorm=[Pr_phi2_NoNorm;Prphi2_S];
    
    %Normalizing Probabilities
    Den=sum(Pr_phi1_NoNorm+Pr_phi2_NoNorm);
    Pr_phi1=Pr_phi1_NoNorm/Den;
    Pr_phi2=Pr_phi2_NoNorm/Den;
    tComputationProbabilities=tComputationProbabilities+toc(TStartProbabilities);
    
    %Candidate Label Functions
    AllLabelFcnValues.phsi=CandidateLabelFcns.phsi_i1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tSTART2=tic;
    %[perfData,tPredicteLabels]=test_StartMeasuringRAMUsage;
    [~,ErrorPartitions,PredictedLabels]=test_ComputeErrorClustererV2(AllLabelFcnValues,CandidateLabelFcns,n1,n2,Pr_phi1,Pr_phi2);
    %[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,tPredicteLabels,'off');
    % Filenames for files that contain probabilities for the optimal operator
    FileName1=['Probabilities1_LabelFcns_SubOptimalPseed_HamDist_' num2str(i) '.def'];
    ModelProbabilities1=fullfile(PathName,FileName1);
    FileName2=['Probabilities2_LabelFcns_SubOptimalPseed_HamDist_' num2str(i) '.def'];
    ModelProbabilities2=fullfile(PathName,FileName2);
    % Saving Probabilities
    save(ModelProbabilities1,'Pr_phi1','-ascii');
    save(ModelProbabilities2,'Pr_phi2','-ascii');
    % Saving Errors of Partitions
    FileName4=['Error_Partitions_SubOptimalPseed_HamDist_' num2str(i) '.def'];
    ModelErrorPartitions=fullfile(PathName,FileName4);
    save(ModelErrorPartitions,'ErrorPartitions','-ascii');
    % Saving Optimal vector of labels
    FileName5=['Predicted_Labels_BayesSubOptPseed_HamDist_' num2str(i) '.def'];
    ModelSubOptimalLabelsPseed=fullfile(PathName,FileName5);
    save(ModelSubOptimalLabelsPseed,'PredictedLabels','-ascii');
    % Saving Label function values
    FileName6=['PossibleLabelFcnValues_' num2str(i) '.mat'];
    ListOfPossibleLabelFunctionValues=fullfile(PathName,FileName6);
    save(ListOfPossibleLabelFunctionValues,'CandidateLabelFcns','-mat');
    % Computing Empirical Error
    Mask=['Predicted_Labels_BayesSubOptPseed_HamDist_' num2str(i)];
    EmpiricalErrorBayesSubOptPseedClusterer(i+1)=test_ComputeEmpiricalClusteringError(GenericExpFolderName,ModelNumber,Mask);
    TheoreticalErrorBayesSubOptPseedClusterer(i+1)=test_ComputeTheoreticalClusteringErrorValue(Pr_phi1,Pr_phi2,PredictedLabels'-1,CandidateLabelFcns,M);
%     InfoAllLabelFcnValues=whos('AllLabelFcnValues');
%     RAMMemory_AllLabelFcnValues=InfoAllLabelFcnValues.bytes/(10^6);
%     InfoPossibleLabelFcnValues=whos('CandidateLabelFcns');
%     RAMMemory_PossibleLabelFcnValues=InfoPossibleLabelFcnValues.bytes/(10^6);
%     InfoProbabilitiesPartitions1=whos('Pr_phi1');
%     RAMMemory_ProbabilitiesPartitions1=InfoProbabilitiesPartitions1.bytes/(10^6);
%     InfoProbabilitiesPartitions2=whos('Pr_phi2');
%     RAMMemory_ProbabilitiesPartitions2=InfoProbabilitiesPartitions2.bytes/(10^6);
% %     RAMMemory_BayesSubOptPseed(i+1)=MaxRAMUsage+RAMMemory_AllLabelFcnValues+...
% %         RAMMemory_PossibleLabelFcnValues+RAMMemory_ProbabilitiesPartitions1+...
% %         RAMMemory_ProbabilitiesPartitions2;
%     RAMMemory_BayesSubOptPseed(i+1)=RAMMemory_AllLabelFcnValues+...
%         RAMMemory_PossibleLabelFcnValues+RAMMemory_ProbabilitiesPartitions1+...
%         RAMMemory_ProbabilitiesPartitions2;
    Time2=toc(tSTART2);
    TimeOfProcessing_BayesSubOptPseed(i+1)=tFindingPseed+tComputationLabelFcnValues+tComputationProbabilities+Time2+TFindingSpaceOfSearch;
end
TimeComputationPossibleLabelFcnValues=tComputationLabelFcnValues;
TimeComputationProbabilitiesPossibleLabelFcnValues=tComputationProbabilities;
return