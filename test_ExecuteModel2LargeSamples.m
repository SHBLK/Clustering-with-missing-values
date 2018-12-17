function test_ExecuteModel2LargeSamples(data)
% test_ExecuteModel2LargeSamples tests the performance of the suboptimal
% Pseed clustering algorithm and compares its results with the ones from
% k-POD and Fuzzy c-means with optimal completion strategy directly applied
% to the data with missing values, and also classical clustering algorithms such as: 
% fuzzy c-means, k-means, hierarchical with single and complete linkage,and random clusterer (applied to imputed data). 
%
% The sets used for testing were
% generated using a random labeled point process (RLPP) composed of two
% Gaussian distributions with normal means and inverse Wishart covariances.
% This function is used when working with larger point sets.
%
% See also test_ExecuteModel0, test_ExecuteModel1, test_ExecuteModel2,
% test_ExecuteModel0LargeSamples, and test_ExecuteModel1LargeSamples


m=data.m;
v=data.v;
k=data.k;
S=data.S;
ProjectName=data.ProjectName;
MissingName=strcat('Missing_Rate_',num2str(100*data.MissingRate)); 
Missing_Rate =  data.MissingRate; 
TypeModel='EXP_NormalMean_InverseWishartCovariance';
GenericExpFolderName=fullfile(TypeModel,ProjectName,MissingName); 
NumberSamples=data.NumberSamples;
NumberTestingSets=data.NumberTestingSets;
% Maximum Hamming distance between candidate partitions and the seed
% partition
MaxHammingDistance=data.MaxHammingDistance;

% Loading data from previous crashed simulation
try
    load(fullfile(GenericExpFolderName,'Backup.mat'));
    fprintf('\r\rA previous not completed simulation was found... Loading old data...\r\r');
catch
    IndicatorSamplesdrawing='off';
    IndicatiorLabelFunctions='off';
    IndicatorBayesOptClusterer='off';
    IndicatorBayesSubClustererPmax='off';
    IndicatorBayesSubClustererPseed='off';
    IndicatorFuzzyCMeansClusterer='off';
    IndicatorKMeansClusterer='off';
    IdicatorHierSiClusterer='off';
    IndicatorHierCoClusterer='off';
    IndicatorRandomClusterer='off';
    IndicatorKPodClusterer='off';
    IndicatorFuzzyCMeansOPSClusterer='off';
    IndicatorErrorPlots='off';
    data3=[];
end


%Verifying if parfor has been activated by the user
try
    ParForMode=lower(data.ParForMode);
catch
    ParForMode='on';
end

% Verifying if the samples drawing has beeen activated by the user
try
    GenerateSamples=lower(data.GenerateSamples);
catch
    GenerateSamples='on';
end

% Computing initial varibales for the simulation
n1=NumberSamples(1,1);
n2=NumberSamples(1,2);
n=sum(NumberSamples);
if mod(n,2)==0
    M=0.5;
else
    M=(n-1)/(2*n);
end

% Initializing timer to measure the total time of processing
tStart = tic;

if strcmp(IndicatorSamplesdrawing,'off')
    if strcmp(GenerateSamples,'on')
        % Drawing Testing Samples
        for i=1:NumberTestingSets
            ModelNumber=num2str(i);
            if i<=NumberTestingSets/2
                n1=NumberSamples(1,1);
                n2=NumberSamples(1,2);
            else
                n1=NumberSamples(1,2);
                n2=NumberSamples(1,1);
            end
            test_SimNormalMeanInverseWishartSigma_withMissing(m,v,k,S,NumberSamples,...
                GenericExpFolderName,Missing_Rate,ModelNumber);
        end
    end
    

end
IndicatorSamplesdrawing='on';
save(fullfile(GenericExpFolderName,'Backup.mat'));


if strcmp(IndicatorBayesSubClustererPseed,'off')
    data2=[];
    % Processing Testing Samples Using Algorithm II (Suboptimal algorithm: Pseed)
    tBayesSubOptPseed = tic;
    parfor i=1:NumberTestingSets
        fprintf('\rAlg. II (Suboptimal: Pseed): Processing dataset N: %d of %d\r',i,NumberTestingSets);
        TSTART=tic;
        %[perfData,TSuboptPseed]=test_StartMeasuringRAMUsage;
        ModelNumber=num2str(i);
        PathName=[GenericExpFolderName '/SAMPLES/Sample_' ModelNumber];
        % Loading feature vectors
        FileName1='Sample_Data.def';
        ModelNameDataSamples=fullfile(PathName,FileName1);
        Vectors=load(ModelNameDataSamples,'Vectors');
        
        % Loading Observed Pattern 
        FileName_='Sample_missing.def'; 
        ModelNameObservedPattern=fullfile(PathName,FileName_); 
        Missing=load(ModelNameObservedPattern,'Missing'); 
        % Loading hyperparameters
        FileName2='Model_Hyperparameters.def';
        ModelNameHypParam=fullfile(PathName,FileName2);
        HypParam = bb_cl_loadparam(ModelNameHypParam);
        n1=HypParam.quantities(1,1);
        n2=HypParam.quantities(1,2);
        [val1,val2,val4,val5,val6]=test_SuboptimalBayesClustererPseedLargeSamplesV2_withMiss_linux(GenericExpFolderName,...
            ModelNumber,TypeModel,Vectors,HypParam,n1,n2,MaxHammingDistance,Missing);
        
        EmpiricalErrorBayesSubOptPseedClusterer(i,:)=val1;
%         TheoreticalErrorBayesSubOptPseedClusterer(i,:)=val2;
%         %[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,TSuboptPseed,'off');
%         RAMMemory_BayesSubOptPseed(i)=0;%MaxRAMUsage;
        TimeOfProcessing_BayesSubOptPseed(i)=toc(TSTART);

        TimeOfProcessing_BayesSubOptPseedVector(i,:)=val4;
        
        TimeComputationPossibleLabelFcnValues(i)=val5;
        TimeComputationProbabilitiesPossibleLabelFcnValues(i)=val6;
%         RAMMemory_PossibleLabelFcnValues(i)=val7;
%         RAMMemory_ProbabilitiesPartitions1(i)=val8;
%         RAMMemory_ProbabilitiesPartitions2(i)=val9;
        
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_BayesSubOptPseed(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_BayesSubOptPseed);
    data2.TotalElapsedTimeOfProcessing=toc(tBayesSubOptPseed);
%     data2.RAMMemoryUsedForEachSet=mean(RAMMemory_BayesSubOptPseed);
    data2.EmpiricalError=mean(EmpiricalErrorBayesSubOptPseedClusterer);
%     data2.TheoreticalError=mean(TheoreticalErrorBayesSubOptPseedClusterer);
    test_SaveFinalReport(GenericExpFolderName,'Results_SubOptimalBayesClusterer_Pseed',data2);
    
    % Summary of results
    data3.EmpiricalError_BayesSubOptimal_Pseed=data2.EmpiricalError;
%     data3.TheoreticalError_BayesSubOptimal_Pseed=data2.TheoreticalError;
end
clear AllLabelFcnValues;
IndicatorBayesSubClustererPseed='on';
save(fullfile(GenericExpFolderName,'Backup.mat'));

if strcmp(IndicatorFuzzyCMeansClusterer,'off')
    data2=[];
    % Processing Testing Samples Using Fuzzy C-Means Algorithm
    tFCM = tic;
    parfor i=1:NumberTestingSets
        fprintf('Alg. FCM: Processing dataset N: %d of %d\r',i,NumberTestingSets);
        TSTART=tic;
        %[perfData,t]=test_StartMeasuringRAMUsage;
        ModelNumber=num2str(i);
        PredictedLabels=test_FuzzyCMeans_withMissing(GenericExpFolderName,ModelNumber)'-1;
        Mask2='Predicted_Labels_FCM';
        EmpiricalErrorFCM(i)=test_ComputeEmpiricalClusteringError(GenericExpFolderName,ModelNumber,Mask2);
        TheoreticalErrorFCM(i)=test_ComputeTheoreticalClusteringErrorLargeSamples(PredictedLabels,M,MaxHammingDistance,GenericExpFolderName,ModelNumber);
%         %[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,t,'off');
%         %         RAMMemory_FCM(i)=MaxRAMUsage+RAMMemory_PossibleLabelFcnValues(i)+RAMMemory_ProbabilitiesPartitions1(i)+RAMMemory_ProbabilitiesPartitions2(i);
%         %         TimeOfProcessing_FCM(i)=toc(TSTART)+TimeComputationPossibleLabelFcnValues(i)+TimeComputationProbabilitiesPossibleLabelFcnValues(i);
%         RAMMemory_FCM(i)=0;%MaxRAMUsage;
        TimeOfProcessing_FCM(i)=toc(TSTART);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_FCM(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_FCM);
    data2.TotalElapsedTimeOfProcessing=toc(tFCM);
%     data2.RAMMemoryUsedForEachSet=mean(RAMMemory_FCM);
    data2.EmpiricalError=mean(EmpiricalErrorFCM);
%     data2.TheoreticalError=mean(TheoreticalErrorFCM);
    test_SaveFinalReport(GenericExpFolderName,'Results_Fuzzy_C_Means',data2);
    
    % Summary of results
    data3.EmpiricalError_FuzzyCMeans=data2.EmpiricalError;
%     data3.TheoreticalError_FuzzyCMeans=data2.TheoreticalError;
end
IndicatorFuzzyCMeansClusterer='on';
save(fullfile(GenericExpFolderName,'Backup.mat'));

if strcmp(IndicatorKMeansClusterer,'off')
    data2=[];
    % Processing Testing Samples Using K-Means Algorithm
    tKM = tic;
    parfor i=1:NumberTestingSets
        fprintf('Alg. KM: Processing dataset N: %d of %d\r',i,NumberTestingSets);
        TSTART=tic;
        %[perfData,t]=test_StartMeasuringRAMUsage;
        ModelNumber=num2str(i);
        PredictedLabels=test_KMeans_withMissing(GenericExpFolderName,ModelNumber)'-1;
        Mask3='Predicted_Labels_KM';
        EmpiricalErrorKM(i)=test_ComputeEmpiricalClusteringError(GenericExpFolderName,ModelNumber,Mask3);
        TheoreticalErrorKM(i)=test_ComputeTheoreticalClusteringErrorLargeSamples(PredictedLabels,M,MaxHammingDistance,GenericExpFolderName,ModelNumber);
%         %[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,t,'off');
%         %         RAMMemory_KM(i)=MaxRAMUsage+RAMMemory_PossibleLabelFcnValues(i)+RAMMemory_ProbabilitiesPartitions1(i)+RAMMemory_ProbabilitiesPartitions2(i);
%         %         TimeOfProcessing_KM(i)=toc(TSTART)+TimeComputationPossibleLabelFcnValues(i)+TimeComputationProbabilitiesPossibleLabelFcnValues(i);
%         RAMMemory_KM(i)=0;%MaxRAMUsage;
        TimeOfProcessing_KM(i)=toc(TSTART);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_KM(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_KM);
    data2.TotalElapsedTimeOfProcessing=toc(tKM);
%     data2.RAMMemoryUsedForEachSet=mean(RAMMemory_KM);
    data2.EmpiricalError=mean(EmpiricalErrorKM);
%     data2.TheoreticalError=mean(TheoreticalErrorKM);
    test_SaveFinalReport(GenericExpFolderName,'Results_K_Means',data2);
    
    % Summary of results
    data3.EmpiricalError_KMeans=data2.EmpiricalError;
%     data3.TheoreticalError_KMeans=data2.TheoreticalError;
end
IndicatorKMeansClusterer='on';
save(fullfile(GenericExpFolderName,'Backup.mat'));

if strcmp(IdicatorHierSiClusterer,'off')
    data2=[];
    % Processing Testing Samples Using Hierarchical-Si Algorithm
    tHIEsi = tic;
    parfor i=1:NumberTestingSets
        fprintf('Alg. HIE-si: Processing dataset N: %d of %d\r',i,NumberTestingSets);
        TSTART=tic;
        %[perfData,t]=test_StartMeasuringRAMUsage;
        ModelNumber=num2str(i);
        TypeDist='eu';
        Method='si';
        PredictedLabels=test_Hierarchical_withMissing(GenericExpFolderName,ModelNumber,TypeDist,Method)'-1;
        Mask4=['Predicted_Labels_Hierarchical_' Method];
        EmpiricalErrorHierarchicalSi(i)=test_ComputeEmpiricalClusteringError(GenericExpFolderName,ModelNumber,Mask4);
        TheoreticalErrorHierarchicalSi(i)=test_ComputeTheoreticalClusteringErrorLargeSamples(PredictedLabels,M,MaxHammingDistance,GenericExpFolderName,ModelNumber);
%         %[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,t,'off');
%         %         RAMMemory_HierarchicalSi(i)=MaxRAMUsage+RAMMemory_PossibleLabelFcnValues(i)+RAMMemory_ProbabilitiesPartitions1(i)+RAMMemory_ProbabilitiesPartitions2(i);
%         %         TimeOfProcessing_HierarchicalSi(i)=toc(TSTART)+TimeComputationPossibleLabelFcnValues(i)+TimeComputationProbabilitiesPossibleLabelFcnValues(i);
%         RAMMemory_HierarchicalSi(i)=0;%MaxRAMUsage;
        TimeOfProcessing_HierarchicalSi(i)=toc(TSTART);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_HierarchicalSi(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_HierarchicalSi);
    data2.TotalElapsedTimeOfProcessing=toc(tHIEsi);
%     data2.RAMMemoryUsedForEachSet=mean(RAMMemory_HierarchicalSi);
    data2.EmpiricalError=mean(EmpiricalErrorHierarchicalSi);
%     data2.TheoreticalError=mean(TheoreticalErrorHierarchicalSi);
    test_SaveFinalReport(GenericExpFolderName,'Results_Hierarchical_Single',data2);
    
    % Summary of results
    data3.EmpiricalError_Hierarchical_Single=data2.EmpiricalError;
%     data3.TheoreticalError_Hierarchical_Single=data2.TheoreticalError;
end
IdicatorHierSiClusterer='on';
save(fullfile(GenericExpFolderName,'Backup.mat'));

if strcmp(IndicatorHierCoClusterer,'off')
    data2=[];
    % Processing Testing Samples Using Hierarchical-Co Algorithm
    tHIEco = tic;
    parfor i=1:NumberTestingSets
        fprintf('Alg. HIE-co: Processing dataset N: %d of %d\r',i,NumberTestingSets);
        TSTART=tic;
        %[perfData,t]=test_StartMeasuringRAMUsage;
        ModelNumber=num2str(i);
        TypeDist='eu';
        Method='co';
        PredictedLabels=test_Hierarchical_withMissing(GenericExpFolderName,ModelNumber,TypeDist,Method)'-1;
        Mask5=['Predicted_Labels_Hierarchical_' Method];
        EmpiricalErrorHierarchicalCo(i)=test_ComputeEmpiricalClusteringError(GenericExpFolderName,ModelNumber,Mask5);
        TheoreticalErrorHierarchicalCo(i)=test_ComputeTheoreticalClusteringErrorLargeSamples(PredictedLabels,M,MaxHammingDistance,GenericExpFolderName,ModelNumber);
%         %[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,t,'off');
%         %         RAMMemory_HierarchicalCo(i)=MaxRAMUsage+RAMMemory_PossibleLabelFcnValues(i)+RAMMemory_ProbabilitiesPartitions1(i)+RAMMemory_ProbabilitiesPartitions2(i);
%         %         TimeOfProcessing_HierarchicalCo(i)=toc(TSTART)+TimeComputationPossibleLabelFcnValues(i)+TimeComputationProbabilitiesPossibleLabelFcnValues(i);
%         RAMMemory_HierarchicalCo(i)=0;%MaxRAMUsage;
        TimeOfProcessing_HierarchicalCo(i)=toc(TSTART);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_HierarchicalCo(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_HierarchicalCo);
    data2.TotalElapsedTimeOfProcessing=toc(tHIEco);
%     data2.RAMMemoryUsedForEachSet=mean(RAMMemory_HierarchicalCo);
    data2.EmpiricalError=mean(EmpiricalErrorHierarchicalCo);
%     data2.TheoreticalError=mean(TheoreticalErrorHierarchicalCo);
    test_SaveFinalReport(GenericExpFolderName,'Results_Hierarchical_Complete',data2);
    
    % Summary of results
    data3.EmpiricalError_Hierarchical_Complete=data2.EmpiricalError;
%     data3.TheoreticalError_Hierarchical_Complete=data2.TheoreticalError;
end
IndicatorHierCoClusterer='on';
save(fullfile(GenericExpFolderName,'Backup.mat'));

if strcmp(IndicatorRandomClusterer,'off')
    data2=[];
    % Processing Testing Samples Using a Random Algorithm
    tRnd = tic;
    parfor i=1:NumberTestingSets
        fprintf('Alg. Random Op.: Processing dataset N: %d of %d\r',i,NumberTestingSets);
        TSTART=tic;
        %[perfData,t]=test_StartMeasuringRAMUsage;
        ModelNumber=num2str(i);
        PredictedLabels=test_RandomClustererV2LargeSamples(GenericExpFolderName,ModelNumber,n)'-1;
        Mask6='Predicted_Labels_RandomOp';
        EmpiricalErrorRandomClusterer(i)=test_ComputeEmpiricalClusteringError(GenericExpFolderName,ModelNumber,Mask6);
        TheoreticalErrorRandomClusterer(i)=test_ComputeTheoreticalClusteringErrorLargeSamples(PredictedLabels,M,MaxHammingDistance,GenericExpFolderName,ModelNumber);
%         %[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,t,'off');
%         %         RAMMemory_RandomClusterer(i)=MaxRAMUsage+RAMMemory_PossibleLabelFcnValues(i)+RAMMemory_ProbabilitiesPartitions1(i)+RAMMemory_ProbabilitiesPartitions2(i);
%         %         TimeOfProcessing_RandomClusterer(i)=toc(TSTART)+TimeComputationPossibleLabelFcnValues(i)+TimeComputationProbabilitiesPossibleLabelFcnValues(i);
%         RAMMemory_RandomClusterer(i)=0;%MaxRAMUsage;
        TimeOfProcessing_RandomClusterer(i)=toc(TSTART);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_RandomClusterer(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_RandomClusterer);
    data2.TotalElapsedTimeOfProcessing=toc(tRnd);
%     data2.RAMMemoryUsedForEachSet=mean(RAMMemory_RandomClusterer);
    data2.EmpiricalError=mean(EmpiricalErrorRandomClusterer);
%     data2.TheoreticalError=mean(TheoreticalErrorRandomClusterer);
    test_SaveFinalReport(GenericExpFolderName,'Results_Random_Clusterer',data2);
    
    % Summary of results
    data3.EmpiricalError_Random=data2.EmpiricalError;
%     data3.TheoreticalError_Random=data2.TheoreticalError;
end
IndicatorRandomClusterer='on';
save(fullfile(GenericExpFolderName,'Backup.mat'));


if strcmp(IndicatorKPodClusterer,'off')
    data2=[];
    % Processing Testing Samples Using K-Means Algorithm
    tKP = tic;
    parfor i=1:NumberTestingSets
        fprintf('Alg. KP: Processing dataset N: %d of %d\r',i,NumberTestingSets);
        TSTART=tic;
       % [perfData,t]=test_StartMeasuringRAMUsage;
        ModelNumber=num2str(i);
        PredictedLabels=test_KPod_withMissing(GenericExpFolderName,ModelNumber)'-1; 
        Mask3='Predicted_Labels_KP';
        EmpiricalErrorKP(i)=test_ComputeEmpiricalClusteringError(GenericExpFolderName,ModelNumber,Mask3);
        TheoreticalErrorKP(i)=test_ComputeTheoreticalClusteringErrorLargeSamples(PredictedLabels,M,MaxHammingDistance,GenericExpFolderName,ModelNumber);
%         MaxRAMUsage=0;%[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,t,'off');
%         %         RAMMemory_KM(i)=MaxRAMUsage+RAMMemory_ValidLabelFcnValues+RAMMemory_ComputeProbabilities(i);
%         %         TimeOfProcessing_KM(i)=toc(TSTART)+TimeCompProbabilities(i);
%         RAMMemory_KP(i)=MaxRAMUsage;
        TimeOfProcessing_KP(i)=toc(TSTART);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_KP(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_KP);
    data2.TotalElapsedTimeOfProcessing=toc(tKP);
%     data2.RAMMemoryUsedForEachSet=mean(RAMMemory_KP);
    data2.EmpiricalError=mean(EmpiricalErrorKP);
%     data2.TheoreticalError=mean(TheoreticalErrorKP);
    test_SaveFinalReport(GenericExpFolderName,'Results_K_Pod',data2);
    
    % Summary of results
    data3.EmpiricalError_KPod=data2.EmpiricalError;
%     data3.TheoreticalError_KPod=data2.TheoreticalError;
end
IndicatorKPodClusterer='on';
save(fullfile(GenericExpFolderName,'Backup.mat'));

if strcmp(IndicatorFuzzyCMeansOPSClusterer,'off')
    data2=[];
    % Processing Testing Samples Using Fuzzy C-Means Algorithm
    tFCM_OPS = tic;
    parfor i=1:NumberTestingSets
        fprintf('Alg. FCM-OPC: Processing dataset N: %d of %d\r',i,NumberTestingSets);
        TSTART=tic;
        %[perfData,t]=test_StartMeasuringRAMUsage;
        ModelNumber=num2str(i);
        PredictedLabels=test_FuzzyCMeansOPS_withMissing(GenericExpFolderName,ModelNumber)'-1; 
        Mask2='Predicted_Labels_FCM_OPS';
        EmpiricalErrorFCM_OPS(i)=test_ComputeEmpiricalClusteringError(GenericExpFolderName,ModelNumber,Mask2);
        TheoreticalErrorFCM_OPS(i)=test_ComputeTheoreticalClusteringErrorLargeSamples(PredictedLabels,M,MaxHammingDistance,GenericExpFolderName,ModelNumber);
%         MaxRAMUsage=0;%[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,t,'off');
%         %         RAMMemory_FCM(i)=MaxRAMUsage+RAMMemory_ValidLabelFcnValues+RAMMemory_ComputeProbabilities(i);
%         %         TimeOfProcessing_FCM(i)=toc(TSTART)+TimeCompProbabilities(i);
%         RAMMemory_FCM_OPS(i)=MaxRAMUsage;
        TimeOfProcessing_FCM_OPS(i)=toc(TSTART);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_FCM(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_FCM_OPS);
    data2.TotalElapsedTimeOfProcessing=toc(tFCM_OPS);
%     data2.RAMMemoryUsedForEachSet=mean(RAMMemory_FCM_OPS);
    data2.EmpiricalError=mean(EmpiricalErrorFCM_OPS);
%     data2.TheoreticalError=mean(TheoreticalErrorFCM_OPS);
    test_SaveFinalReport(GenericExpFolderName,'Results_Fuzzy_C_Means_OPS',data2);
    
    % Summary of results
    data3.EmpiricalError_FuzzyCMeans_OPS=data2.EmpiricalError;
%     data3.TheoreticalError_FuzzyCMeans_OPS=data2.TheoreticalError;
end

data3.TotalTimeOfProcessing=toc(tStart);

IndicatorFuzzyCMeansOPSClusterer='on';
save(fullfile(GenericExpFolderName,'Backup.mat'));

if strcmp(IndicatorErrorPlots,'off')
    % Saving Results
    test_SaveFinalReport(GenericExpFolderName,'Empirical and Theoretical Errors',data3);
    
    %Saving Simulation settings
    test_SaveFinalReport(GenericExpFolderName,'Results_SimulationSettings',data);
    

    % Saving empirical errors as a function of the Hamming distance
    EE.EmpiricalError_SubOptimalPseed=mean(EmpiricalErrorBayesSubOptPseedClusterer);
    EE.EmpiricalError_FCM=mean(EmpiricalErrorFCM);
    EE.EmpiricalError_KM=mean(EmpiricalErrorKM);
    EE.EmpiricalError_HieSi=mean(EmpiricalErrorHierarchicalSi);
    EE.EmpiricalError_HieCo=mean(EmpiricalErrorHierarchicalCo);
    EE.EmpiricalError_Random=mean(EmpiricalErrorRandomClusterer);
    
    EE.EmpiricalError_FCM_OPS=mean(EmpiricalErrorFCM_OPS);
    EE.EmpiricalError_KP=mean(EmpiricalErrorKP);
    
    EE.EmpiricalError_Std_SubOptimalPseed=std(EmpiricalErrorBayesSubOptPseedClusterer);
    EE.EmpiricalError_Std_FCM=std(EmpiricalErrorFCM);
    EE.EmpiricalError_Std_KM=std(EmpiricalErrorKM);
    EE.EmpiricalError_Std_HieSi=std(EmpiricalErrorHierarchicalSi);
    EE.EmpiricalError_Std_HieCo=std(EmpiricalErrorHierarchicalCo);
    EE.EmpiricalError_Std_Random=std(EmpiricalErrorRandomClusterer);
    
    EE.EmpiricalError_Std_FCM_OPS=std(EmpiricalErrorFCM_OPS);
    EE.EmpiricalError_Std_KP=std(EmpiricalErrorKP);
    
    test_SaveFinalReport(GenericExpFolderName,'Results_EmpiricalError_HammingDistance',EE);
    
    % Plotting empirical errors
    test_PlotResultsLargeSamples_new(GenericExpFolderName,'Results_EmpiricalError_HammingDistance');
    



end
IndicatorErrorPlots='on';

save(fullfile(GenericExpFolderName,'REPORT','CommandWindowVariables.mat'));
delete(fullfile(GenericExpFolderName,'Backup.mat'));
return