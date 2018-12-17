function [EmpiricalErrorBayesSubOptPmaxClusterer,...
    TheoreticalErrorBayesSubOptPmaxClusterer,...
    TimeOfProcessing_BayesSubOptPmax,...
    SizeVectorOfProbabilities]=test_SuboptimalBayesClustererPmax_withMiss_linux(GenericExpFolderName,ModelNumber,...
    TypeModel,Vectors,HypParam,n1,n2,ValidLabelFcnValues,AllLabelFcnValues,MaxHammingDistance,MissMatrix)
% test_SuboptimalBayesClustererPmax implements the suboptimal Pmax
% clustering algorithm

%
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
    [Prphi1_S,Prphi2_S]=test_ComputeProbabilityOfPartitions_withMiss(TypeModel,Vectors,HypParam,ValidLabelFcnValues,MissMatrix);
end
n=n1+n2;
if mod(n,2)==0
    M=0.5;
else
    M=(n-1)/(2*n);
end
InfoProbabilitiesPartitions1=whos('Prphi1_S');
RAMMemory_ProbabilitiesPartitions1=InfoProbabilitiesPartitions1.bytes/(10^6);
InfoProbabilitiesPartitions2=whos('Prphi2_S');
RAMMemory_ProbabilitiesPartitions2=InfoProbabilitiesPartitions2.bytes/(10^6);
tSTART1=tic;
%[perfData,t]=test_StartMeasuringRAMUsage;
[~,IdxMaxProbability]=max(Prphi1_S+Prphi2_S);
% Finding the highest probability partition
MaxProbabilityPartition=ValidLabelFcnValues.phsi_i1(IdxMaxProbability,:);
% Computing the suboptimal partition using a set of partitions with a given 
% hamming distance respect the highest probability partition
%[~,RAMMemory_FindingPmax]=test_FinishMeasuringRAMUsage(perfData,t,'off');
Time1=toc(tSTART1);
EmpiricalErrorBayesSubOptPseedClusterer=[];
TheoreticalErrorBayesSubOptPseedClusterer=[];
for i=0:MaxHammingDistance
    if i~=MaxHammingDistance
        IdxOfCandidateLabelFcns=test_FindLabelFcnsHammingDistance(ValidLabelFcnValues.phsi_i1,MaxProbabilityPartition,i,'LesserOrEqual','off');
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
    tSTART2=tic;
    %[perfData,tPredictedLabels]=test_StartMeasuringRAMUsage;
    [~,ErrorPartitions,PredictedLabels]=test_ComputeErrorClustererV2(AllLabelFcnValues,CandidateLabelFcns,n1,n2,Pr_phi1,Pr_phi2);
   % [~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,tPredictedLabels,'off');
    % Filenames for files that contain probabilities for the optimal operator
    FileName1=['Probabilities1_LabelFcns_SubOptimalPmax_HamDist_' num2str(i) '.def'];
    ModelProbabilities1=fullfile(PathName,FileName1);
    FileName2=['Probabilities2_LabelFcns_SubOptimalPmax_HamDist_' num2str(i) '.def'];
    ModelProbabilities2=fullfile(PathName,FileName2);
    % Saving Probabilities
    save(ModelProbabilities1,'Pr_phi1','-ascii');
    save(ModelProbabilities2,'Pr_phi2','-ascii');
    % Saving Indices of partitions with a given Hamming Distance
    FileName3=['Indices_Partitions_SubOptimalPmax_HamDist_' num2str(i) '.def'];
    ModelIndices=fullfile(PathName,FileName3);
    save(ModelIndices,'IdxOfCandidateLabelFcns','-ascii');
    % Saving Errors of Partitions
    FileName4=['Error_Partitions_SubOptimalPmax_HamDist_' num2str(i) '.def'];
    ModelErrorPartitions=fullfile(PathName,FileName4);
    save(ModelErrorPartitions,'ErrorPartitions','-ascii');
    % Saving Optimal vector of labels
    FileName5=['Predicted_Labels_BayesSubOptPmax_HamDist_' num2str(i) '.def'];
    ModelSubOptimalLabelsPmax=fullfile(PathName,FileName5);
    save(ModelSubOptimalLabelsPmax,'PredictedLabels','-ascii');
    % Computing Empirical Error
    Mask=['Predicted_Labels_BayesSubOptPmax_HamDist_' num2str(i)];
    EmpiricalErrorBayesSubOptPmaxClusterer(i+1)=test_ComputeEmpiricalClusteringError(GenericExpFolderName,ModelNumber,Mask);
    TheoreticalErrorBayesSubOptPmaxClusterer(i+1)=test_ComputeTheoreticalClusteringErrorValue(Pr_phi1,Pr_phi2,PredictedLabels'-1,CandidateLabelFcns,M);
    InfoAllLabelFcnValues=whos('AllLabelFcnValues');
    RAMMemory_AllLabelFcnValues=InfoAllLabelFcnValues.bytes/(10^6);
    InfoPossibleLabelFcnValues=whos('CandidateLabelFcns');
    RAMMemory_PossibleLabelFcnValues=InfoPossibleLabelFcnValues.bytes/(10^6);
% %     RAMMemory_BayesSubOptPmax(i+1)=MaxRAMUsage+RAMMemory_AllLabelFcnValues+...
% %         RAMMemory_PossibleLabelFcnValues+RAMMemory_ProbabilitiesPartitions1+...
% %         RAMMemory_ProbabilitiesPartitions2+RAMMemory_FindingPmax;
%     RAMMemory_BayesSubOptPmax(i+1)=RAMMemory_AllLabelFcnValues+...
%         RAMMemory_PossibleLabelFcnValues+RAMMemory_ProbabilitiesPartitions1+...
%         RAMMemory_ProbabilitiesPartitions2;
    Time2=toc(tSTART2);
    TimeOfProcessing_BayesSubOptPmax(i+1)=Time1+Time2;
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