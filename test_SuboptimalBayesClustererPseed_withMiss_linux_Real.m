function [EmpiricalErrorBayesSubOptPseedClusterer,...
    TheoreticalErrorBayesSubOptPseedClusterer,...
    TimeOfProcessing_BayesSubOptPseed,...
    SizeVectorOfProbabilities]=test_SuboptimalBayesClustererPseed_withMiss_linux_Real(GenericExpFolderName,ModelNumber,...
    TypeModel,Vectors,HypParam,n1,n2,ValidLabelFcnValues,AllLabelFcnValues,MaxHammingDistance,MissMatr)
% test_SuboptimalBayesClustererPseed implements the suboptimal Pseed
% clustering algorithm for small point sets (e.g. sets with maximum 25 points)
%

%
% To process large point sets use test_SuboptimalBayesClustererPseedLargeSamplesV2 

PathName=[GenericExpFolderName '/SAMPLES/Sample_' ModelNumber];
try
    % Filenames for files that contain probabilities for the optimal operator
    FileName1='Probabilities1_LabelFcns_Optimal.def';
    ModelProbabilities1=fullfile(PathName,FileName1);
    FileName2='Probabilities2_LabelFcns_Optimal.def';
    ModelProbabilities2=fullfile(PathName,FileName2);
    % Saving Probabilities
    Prphi1_S=load(ModelProbabilities1,'Prphi1_S');
    Prphi2_S=load(ModelProbabilities2,'Prphi2_S');
    clear ModelProbabilities1;
    clear ModelProbabilities2;
catch
    % Computing probabilities of all possible partitions
    [Prphi1_S,Prphi2_S]=test_ComputeProbabilityOfPartitions_withMiss(TypeModel,Vectors,HypParam,ValidLabelFcnValues,MissMatr);
end
%Computing the partition seed for suboptimal algorithm
n=n1+n2;
if mod(n,2)==0
    M=0.5;
    HammingDistance=2;
else
    M=(n-1)/(2*n);
    HammingDistance=3;
end
tSTART1=tic;
[PartitionSeed,~,~]=test_GeneratePartitionSeed_withMiss_Real(TypeModel,Vectors,HypParam,n1,n2,HammingDistance,MissMatr,5);
Time1=toc(tSTART1);
EmpiricalErrorBayesSubOptPseedClusterer=[];
TheoreticalErrorBayesSubOptPseedClusterer=[];
for i=0:MaxHammingDistance
    tSTART2=tic;
    if i~=MaxHammingDistance
        IdxOfCandidateLabelFcns=test_FindLabelFcnsHammingDistance(ValidLabelFcnValues.phsi_i1,PartitionSeed,i,'LesserOrEqual','off');
        SizeVectorOfProbabilities(i+1)=length(IdxOfCandidateLabelFcns);
        CandidateLabelFcns.phsi_i1=ValidLabelFcnValues.phsi_i1(IdxOfCandidateLabelFcns,:);
        CandidateLabelFcns.phsi_i2=ValidLabelFcnValues.phsi_i2(IdxOfCandidateLabelFcns,:);
        Pr_phi1=Prphi1_S(IdxOfCandidateLabelFcns);
        Pr_phi2=Prphi2_S(IdxOfCandidateLabelFcns);
        %Normalizing Probabilities
        Den=sum(Pr_phi1+Pr_phi2);
        Pr_phi1=Pr_phi1/Den;
        Pr_phi2=Pr_phi2/Den;
    else
        IdxOfCandidateLabelFcns=(1:size(ValidLabelFcnValues.phsi_i1,1))';
        SizeVectorOfProbabilities(i+1)=length(IdxOfCandidateLabelFcns);
        CandidateLabelFcns.phsi_i1=ValidLabelFcnValues.phsi_i1;
        CandidateLabelFcns.phsi_i2=ValidLabelFcnValues.phsi_i2;
        Pr_phi1=Prphi1_S(IdxOfCandidateLabelFcns);
        Pr_phi2=Prphi2_S(IdxOfCandidateLabelFcns);
    end
    if isempty(AllLabelFcnValues)==1
        AllLabelFcnValues.phsi=CandidateLabelFcns.phsi_i1;
    end
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
    % Saving Indices of partitions with a given Hamming Distance
    FileName3=['Indices_Partitions_SubOptimalPseed_HamDist_' num2str(i) '.def'];
    ModelIndices=fullfile(PathName,FileName3);
    save(ModelIndices,'IdxOfCandidateLabelFcns','-ascii');
    % Saving Errors of Partitions
    FileName4=['Error_Partitions_SubOptimalPseed_HamDist_' num2str(i) '.def'];
    ModelErrorPartitions=fullfile(PathName,FileName4);
    save(ModelErrorPartitions,'ErrorPartitions','-ascii');
    % Saving Optimal vector of labels
    FileName5=['Predicted_Labels_BayesSubOptPseed_HamDist_' num2str(i) '.def'];
    ModelSubOptimalLabelsPseed=fullfile(PathName,FileName5);
    save(ModelSubOptimalLabelsPseed,'PredictedLabels','-ascii');
    % Computing Empirical Error
    Mask=['Predicted_Labels_BayesSubOptPseed_HamDist_' num2str(i)];
    EmpiricalErrorBayesSubOptPseedClusterer(i+1)=test_ComputeEmpiricalClusteringError(GenericExpFolderName,ModelNumber,Mask);
    TheoreticalErrorBayesSubOptPseedClusterer(i+1)=test_ComputeTheoreticalClusteringErrorValue(Pr_phi1,Pr_phi2,PredictedLabels'-1,CandidateLabelFcns,M);
    InfoAllLabelFcnValues=whos('AllLabelFcnValues');
    RAMMemory_AllLabelFcnValues=InfoAllLabelFcnValues.bytes/(10^6);
    InfoPossibleLabelFcnValues=whos('CandidateLabelFcns');
    RAMMemory_PossibleLabelFcnValues=InfoPossibleLabelFcnValues.bytes/(10^6);
    InfoProbabilitiesPartitions1=whos('Pr_phi1');
    RAMMemory_ProbabilitiesPartitions1=InfoProbabilitiesPartitions1.bytes/(10^6);
    InfoProbabilitiesPartitions2=whos('Pr_phi2');
    RAMMemory_ProbabilitiesPartitions2=InfoProbabilitiesPartitions2.bytes/(10^6);
% %     RAMMemory_BayesSubOptPseed(i+1)=MaxRAMUsage+RAMMemory_AllLabelFcnValues+...
% %         RAMMemory_PossibleLabelFcnValues+RAMMemory_ProbabilitiesPartitions1+...
% %         RAMMemory_ProbabilitiesPartitions2;
%     RAMMemory_BayesSubOptPseed(i+1)=RAMMemory_AllLabelFcnValues+...
%         RAMMemory_PossibleLabelFcnValues+RAMMemory_ProbabilitiesPartitions1+...
%         RAMMemory_ProbabilitiesPartitions2;
    Time2=toc(tSTART2);
    TimeOfProcessing_BayesSubOptPseed(i+1)=Time2;
    % Deleting variables
    clear IdxOfCandidateLabelFcns;
    clear CandidateLabelFcns;
    clear AllCandidateLabelFcns.phsi;
    clear Pr_phi1;
    clear Pr_phi2;
    clear ErrorPartitions;
    clear PredictedLabels;
    clear Time2;
    clear MaxRAMUsage;
end
return