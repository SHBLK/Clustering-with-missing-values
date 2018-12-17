function test_ExecuteModel2Real(data,sample_real0,sample_real1)
% test_ExecuteModel2 tests the performance of the optimal (Bayes),
% suboptimal Pmax, and suboptimal Pseed clustering algorithms for data with missing values and compares
% their results with k-POD and Fuzzy c-means with optimal completion strategy directly applied to the data with missing values.  
% In addition, the comparison with classical clustering
% algorithms on imputed data such as: fuzzy c-means, k-means, hierarchical with single and
% complete linkage,and random clusterer is done. The sets used for testing were
% generated using a random sampling scheme from real normalized 2 component rna-seq data.
% This function is used when working with small point sets.
%


calibSize = data.CalibrationSize; 

ProjectName=data.ProjectName;
MissingName=strcat('Missing_Rate_',num2str(100*data.MissingRate)); 
Missing_Rate =  data.MissingRate; 
TypeModel='EXP_NormalMean_InverseWishartCovariance';
GenericExpFolderName=fullfile(TypeModel,ProjectName,MissingName); 

NumberSamples=data.NumberSamples;
NumberTestingSets=data.NumberTestingSets;

reproduce=data.flag_reproduce;

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

if reproduce
   load(fullfile(GenericExpFolderName,'REPORT','CommandWindowVariables.mat'));
   IndicatorErrorPlots='off';
end

%Verifying if parfor has been activated by the user
try
    ParForMode=lower(data.ParForMode);
catch
    ParForMode='on';
end

%Verifying if the samples drawing has beeen activated by the user
try
    GenerateSamples=lower(data.GenerateSamples);
catch
    GenerateSamples='on';
end

% Computing varibales for simulation
n=sum(NumberSamples);

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
            test_RealDataSampling_withMissing(sample_real0,sample_real1,calibSize,NumberSamples,...
                GenericExpFolderName,Missing_Rate,ModelNumber);
        end
    end
    

end
IndicatorSamplesdrawing='on';
save(fullfile(GenericExpFolderName,'Backup.mat'));

if strcmp(IndicatiorLabelFunctions,'off')
    % Computing Label Functions
    n1=NumberSamples(1,1);
    n2=NumberSamples(1,2);
    [AllLabelFcnValues,ValidLabelFcnValues]=test_ComputeLabelFcnsV2(n1,n2,GenericExpFolderName,ProjectName);
%     InfoAllLabelFcnValues=whos('AllLabelFcnValues');
%     RAMMemory_AllLabelFcnValues=InfoAllLabelFcnValues.bytes/(10^6);
%     InfoValidLabelFcnValues=whos('ValidLabelFcnValues');
%     RAMMemory_ValidLabelFcnValues=InfoValidLabelFcnValues.bytes/(10^6);
end
IndicatiorLabelFunctions='on';
save(fullfile(GenericExpFolderName,'Backup.mat'));






if strcmp(IndicatorBayesOptClusterer,'off')
    data2=[];
    % Processing Testing Samples Using Algorithm II (Bayes Algorithm)
    tBayes = tic;
    n=n1+n2;
    if mod(n,2)==0
        M=0.5;
    else
        M=(n-1)/(2*n);
    end
    % Maximum Hamming distance between the candidate partitions
    MaxHammingDistance=M*n;
    Indicator=zeros(1,M*n+1);
    FreqHammingDistancePredictedPmaxPartitions=zeros(1,MaxHammingDistance+1);
    
    %Added codes for parallel
    Prphi1_S_cell = cell(1,NumberTestingSets);
    Prphi2_S_cell = cell(1,NumberTestingSets);
    %MaxRAMUsage_cell = cell(1,NumberTestingSets);
    ErrorPartitions_cell = cell(1,NumberTestingSets);
    PredictedLabels_cell = cell(1,NumberTestingSets);
    val1 = zeros(1,NumberTestingSets);
    val2 = zeros(1,NumberTestingSets);
    parfor i=1:NumberTestingSets
        fprintf('\rAlg. II: Processing dataset N: %d of %d\r',i,NumberTestingSets);
        TSTART=tic;
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
        TSTARTCompProbabilities=tic;
        disp('Working...')
        tic
        [Prphi1_S_cell{i},Prphi2_S_cell{i}]=test_ComputeProbabilityNormalMeanInverseWishartSigma_withMiss(Vectors,HypParam,ValidLabelFcnValues,Missing);
        toc

        TimeCompProbabilities(i)=toc(TSTARTCompProbabilities);
        [val1(i),ErrorPartitions_cell{i},PredictedLabels_cell{i},val2(i)]=test_ComputeErrorClustererV2(AllLabelFcnValues,ValidLabelFcnValues,n1,n2,Prphi1_S_cell{i},Prphi2_S_cell{i});

        TimeOfProcessing_BayesOpt(i)=toc(TSTART);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_BayesOpt(i));
    end
    
    %Post Loop processing for parallel 
    
    for i=1:NumberTestingSets
        ModelNumber=num2str(i);
        PathName=[GenericExpFolderName '/SAMPLES/Sample_' ModelNumber];
        Prphi1_S = Prphi1_S_cell{i};
        Prphi2_S = Prphi2_S_cell{i};
        FileName3='Probabilities1_LabelFcns_Optimal.def';
        ModelProbabilities1=fullfile(PathName,FileName3);
        FileName4='Probabilities2_LabelFcns_Optimal.def';
        ModelProbabilities2=fullfile(PathName,FileName4);
        % Saving Probabilities
        save(ModelProbabilities1,'Prphi1_S','-ascii');
        save(ModelProbabilities2,'Prphi2_S','-ascii');
        
        FileName5='Error_Partitions_Optimal.def';
        ModelErrorPartitions=fullfile(PathName,FileName5);
        ErrorPartitions=ErrorPartitions_cell{i};
        save(ModelErrorPartitions,'ErrorPartitions','-ascii');
 
        FileName6='Predicted_Labels_BayesOpt.def';
        ModelOptimalLabels=fullfile(PathName,FileName6);
        PredictedLabels=PredictedLabels_cell{i};
        save(ModelOptimalLabels,'PredictedLabels','-ascii');
        
        Indicator(1,val1(i)+1)=Indicator(1,val1(i)+1)+1;
        FreqHammingDistancePredictedPmaxPartitions(1,val2(i)+1)=FreqHammingDistancePredictedPmaxPartitions(1,val2(i)+1)+1;
        
%         MaxRAMUsage = 0;%MaxRAMUsage_cell{i};
%         InfoProbabilitiesPartitions1=whos('Prphi1_S');
%         RAMMemory_ProbabilitiesPartitions1=InfoProbabilitiesPartitions1.bytes/(10^6);
%         InfoProbabilitiesPartitions2=whos('Prphi2_S');
%         RAMMemory_ProbabilitiesPartitions2=InfoProbabilitiesPartitions2.bytes/(10^6);
%         RAMMemory_ComputeProbabilities(i)=RAMMemory_ProbabilitiesPartitions1+RAMMemory_ProbabilitiesPartitions2;
%         
%         RAMMemory_BayesOpt(i)=MaxRAMUsage+RAMMemory_AllLabelFcnValues+RAMMemory_ValidLabelFcnValues+RAMMemory_ComputeProbabilities(i);
        
        Mask1='Predicted_Labels_BayesOpt';
        EmpiricalErrorBayesOptClusterer(i)=test_ComputeEmpiricalClusteringError(GenericExpFolderName,ModelNumber,Mask1);
%        TheoreticalErrorBayesOptClusterer(i)=test_ComputeTheoreticalClusteringError(PredictedLabels'-1,ValidLabelFcnValues,M,GenericExpFolderName,ModelNumber);
% 
    end
    clear Prphi1_S_cell;
    clear Prphi1_S;
    clear Prphi2_S_cell;
    clear Prphi2_S;
   % clear MaxRAMUsage_cell;
    clear PredictedLabels_cell;
    clear ErrorPartitions_cell;
    clear val1;
    clear val2;

    
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_BayesOpt);
    data2.TotalElapsedTimeOfProcessing=toc(tBayes);
    %data2.RAMMemoryUsedForEachSet=mean(RAMMemory_BayesOpt);
    data2.EmpiricalError=mean(EmpiricalErrorBayesOptClusterer);
    %data2.TheoreticalError=mean(TheoreticalErrorBayesOptClusterer);
    test_SaveFinalReport(GenericExpFolderName,'Results_OptimalBayesClusterer',data2);
    
    % Summary of results
    data3.EmpiricalError_BayesOptimal=data2.EmpiricalError;
    %data3.TheoreticalError_BayesOptimal=data2.TheoreticalError;
end
try
    SpaceSearchSuboptimal=data.SpaceSearchSuboptimal;
    if strcmp(data.SpaceSearchSuboptimal,'Possible')
        clear AllLabelFcnValues;
        AllLabelFcnValues=[];
        RAMMemory_AllLabelFcnValues=0;
    else
        fprintf('\r\rWARNING: All candidate partitions are going to be tested with the Suboptimal Algorithms\r\r');
    end
catch
    clear AllLabelFcnValues;
    AllLabelFcnValues=[];
    RAMMemory_AllLabelFcnValues=0;
end
IndicatorBayesOptClusterer='on';
save(fullfile(GenericExpFolderName,'Backup.mat'));

if strcmp(IndicatorBayesSubClustererPmax,'off')
    data2=[];
    % Processing Testing Samples Using Algorithm II (Suboptimal algorithm: Pmax)
    tBayesSubOptPmax = tic;
    parfor i=1:NumberTestingSets
        fprintf('\rAlg. II (Suboptimal: Pmax): Processing dataset N: %d of %d\r',i,NumberTestingSets);
        TSTART=tic;
        
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
        [val1,val2,val4,val5]=test_SuboptimalBayesClustererPmax_withMiss_linux(GenericExpFolderName,...
            ModelNumber,TypeModel,Vectors,HypParam,n1,n2,ValidLabelFcnValues,AllLabelFcnValues,MaxHammingDistance,Missing); 
        EmpiricalErrorBayesSubOptPmaxClusterer(i,:)=val1;
        %TheoreticalErrorBayesSubOptPmaxClusterer(i,:)=val2;
        %MaxRAMUsage=0;%[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,TSuboptPmax,'off');
        %RAMMemory_BayesSubOptPmax(i)=MaxRAMUsage+RAMMemory_AllLabelFcnValues+RAMMemory_ValidLabelFcnValues+RAMMemory_ComputeProbabilities(i);
        TimeOfProcessing_BayesSubOptPmax(i)=toc(TSTART);
        %RAMMemory_BayesSubOptPmaxVector(i,:)=val3;
        TimeOfProcessing_BayesSubOptPmaxVector(i,:)=val4+TimeCompProbabilities(i);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_BayesSubOptPmax(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_BayesSubOptPmax);
    data2.TotalElapsedTimeOfProcessing=toc(tBayesSubOptPmax);
    %data2.RAMMemoryUsedForEachSet=mean(RAMMemory_BayesSubOptPmax);
    data2.EmpiricalError=mean(EmpiricalErrorBayesSubOptPmaxClusterer);
    %data2.TheoreticalError=mean(TheoreticalErrorBayesSubOptPmaxClusterer);
    test_SaveFinalReport(GenericExpFolderName,'Results_SubOptimalBayesClusterer_Pmax',data2);
    
    % Summary of results
    data3.EmpiricalError_BayesSubOptimal_Pmax=data2.EmpiricalError;
    %data3.TheoreticalError_BayesSubOptimal_Pmax=data2.TheoreticalError;
end
IndicatorBayesSubClustererPmax='on';
save(fullfile(GenericExpFolderName,'Backup.mat'));

if strcmp(IndicatorBayesSubClustererPseed,'off')
    data2=[];
    % Processing Testing Samples Using Algorithm II (Suboptimal algorithm: Pseed)
    tBayesSubOptPseed = tic;
    % Maximum Hamming distance between candidate partitions and the seed
    % partition
    MaxHammingDistance=M*n;
    parfor i=1:NumberTestingSets
        fprintf('\rAlg. II (Suboptimal: Pseed): Processing dataset N: %d of %d\r',i,NumberTestingSets);
        TSTART=tic;
       % [perfData,TSuboptPseed]=test_StartMeasuringRAMUsage;
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
        [val1,val2,val4,val5]=test_SuboptimalBayesClustererPseed_withMiss_linux_Real(GenericExpFolderName,...
            ModelNumber,TypeModel,Vectors,HypParam,n1,n2,ValidLabelFcnValues,[],MaxHammingDistance,Missing);
        EmpiricalErrorBayesSubOptPseedClusterer(i,:)=val1;
        %TheoreticalErrorBayesSubOptPseedClusterer(i,:)=val2;
        %MaxRAMUsage=0;%[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,TSuboptPseed,'off');
        %RAMMemory_BayesSubOptPseed(i)=MaxRAMUsage+RAMMemory_AllLabelFcnValues+RAMMemory_ValidLabelFcnValues+RAMMemory_ComputeProbabilities(i)*val5(end)/max(val5);
        TimeOfProcessing_BayesSubOptPseed(i)=toc(TSTART);
        %RAMMemory_BayesSubOptPseedVector(i,:)=val3;
        TimeOfProcessing_BayesSubOptPseedVector(i,:)=val4+TimeCompProbabilities(i)*val5/max(val5);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_BayesSubOptPseed(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_BayesSubOptPseed);
    data2.TotalElapsedTimeOfProcessing=toc(tBayesSubOptPseed);
    %data2.RAMMemoryUsedForEachSet=mean(RAMMemory_BayesSubOptPseed);
    data2.EmpiricalError=mean(EmpiricalErrorBayesSubOptPseedClusterer);
    %data2.TheoreticalError=mean(TheoreticalErrorBayesSubOptPseedClusterer);
    test_SaveFinalReport(GenericExpFolderName,'Results_SubOptimalBayesClusterer_Pseed',data2);
    
    % Summary of results
    data3.EmpiricalError_BayesSubOptimal_Pseed=data2.EmpiricalError;
    %data3.TheoreticalError_BayesSubOptimal_Pseed=data2.TheoreticalError;
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
        TheoreticalErrorFCM(i)=test_ComputeTheoreticalClusteringError(PredictedLabels,ValidLabelFcnValues,M,GenericExpFolderName,ModelNumber);
%         MaxRAMUsage=0;%[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,t,'off');
%         %         RAMMemory_FCM(i)=MaxRAMUsage+RAMMemory_ValidLabelFcnValues+RAMMemory_ComputeProbabilities(i);
%         %         TimeOfProcessing_FCM(i)=toc(TSTART)+TimeCompProbabilities(i);
%         RAMMemory_FCM(i)=MaxRAMUsage;
        TimeOfProcessing_FCM(i)=toc(TSTART);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_FCM(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_FCM);
    data2.TotalElapsedTimeOfProcessing=toc(tFCM);
    %data2.RAMMemoryUsedForEachSet=mean(RAMMemory_FCM);
    data2.EmpiricalError=mean(EmpiricalErrorFCM);
    %data2.TheoreticalError=mean(TheoreticalErrorFCM);
    test_SaveFinalReport(GenericExpFolderName,'Results_Fuzzy_C_Means',data2);
    
    % Summary of results
    data3.EmpiricalError_FuzzyCMeans=data2.EmpiricalError;
    %data3.TheoreticalError_FuzzyCMeans=data2.TheoreticalError;
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
       % [perfData,t]=test_StartMeasuringRAMUsage;
        ModelNumber=num2str(i);
        PredictedLabels=test_KMeans_withMissing(GenericExpFolderName,ModelNumber)'-1; 
        Mask3='Predicted_Labels_KM';
        EmpiricalErrorKM(i)=test_ComputeEmpiricalClusteringError(GenericExpFolderName,ModelNumber,Mask3);
        TheoreticalErrorKM(i)=test_ComputeTheoreticalClusteringError(PredictedLabels,ValidLabelFcnValues,M,GenericExpFolderName,ModelNumber);
%         MaxRAMUsage=0;%[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,t,'off');
%         %         RAMMemory_KM(i)=MaxRAMUsage+RAMMemory_ValidLabelFcnValues+RAMMemory_ComputeProbabilities(i);
%         %         TimeOfProcessing_KM(i)=toc(TSTART)+TimeCompProbabilities(i);
%         RAMMemory_KM(i)=MaxRAMUsage;
        TimeOfProcessing_KM(i)=toc(TSTART);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_KM(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_KM);
    data2.TotalElapsedTimeOfProcessing=toc(tKM);
    %data2.RAMMemoryUsedForEachSet=mean(RAMMemory_KM);
    data2.EmpiricalError=mean(EmpiricalErrorKM);
    %data2.TheoreticalError=mean(TheoreticalErrorKM);
    test_SaveFinalReport(GenericExpFolderName,'Results_K_Means',data2);
    
    % Summary of results
    data3.EmpiricalError_KMeans=data2.EmpiricalError;
    %data3.TheoreticalError_KMeans=data2.TheoreticalError;
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
        TheoreticalErrorHierarchicalSi(i)=test_ComputeTheoreticalClusteringError(PredictedLabels,ValidLabelFcnValues,M,GenericExpFolderName,ModelNumber);
%         MaxRAMUsage=0;%[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,t,'off');
%         %         RAMMemory_HierarchicalSi(i)=MaxRAMUsage+RAMMemory_ValidLabelFcnValues+RAMMemory_ComputeProbabilities(i);
%         %         TimeOfProcessing_HierarchicalSi(i)=toc(TSTART)+TimeCompProbabilities(i);
%         RAMMemory_HierarchicalSi(i)=MaxRAMUsage;
        TimeOfProcessing_HierarchicalSi(i)=toc(TSTART);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_HierarchicalSi(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_HierarchicalSi);
    data2.TotalElapsedTimeOfProcessing=toc(tHIEsi);
    %data2.RAMMemoryUsedForEachSet=mean(RAMMemory_HierarchicalSi);
    data2.EmpiricalError=mean(EmpiricalErrorHierarchicalSi);
    %data2.TheoreticalError=mean(TheoreticalErrorHierarchicalSi);
    test_SaveFinalReport(GenericExpFolderName,'Results_Hierarchical_Single',data2);
    
    % Summary of results
    data3.EmpiricalError_Hierarchical_Single=data2.EmpiricalError;
    %data3.TheoreticalError_Hierarchical_Single=data2.TheoreticalError;
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
        TheoreticalErrorHierarchicalCo(i)=test_ComputeTheoreticalClusteringError(PredictedLabels,ValidLabelFcnValues,M,GenericExpFolderName,ModelNumber);
%         MaxRAMUsage=0;%[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,t,'off');
%         %         RAMMemory_HierarchicalCo(i)=MaxRAMUsage+RAMMemory_ValidLabelFcnValues+RAMMemory_ComputeProbabilities(i);
%         %         TimeOfProcessing_HierarchicalCo(i)=toc(TSTART)+TimeCompProbabilities(i);
%         RAMMemory_HierarchicalCo(i)=MaxRAMUsage;
        TimeOfProcessing_HierarchicalCo(i)=toc(TSTART);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_HierarchicalCo(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_HierarchicalCo);
    data2.TotalElapsedTimeOfProcessing=toc(tHIEco);
    %data2.RAMMemoryUsedForEachSet=mean(RAMMemory_HierarchicalCo);
    data2.EmpiricalError=mean(EmpiricalErrorHierarchicalCo);
    %data2.TheoreticalError=mean(TheoreticalErrorHierarchicalCo);
    test_SaveFinalReport(GenericExpFolderName,'Results_Hierarchical_Complete',data2);
    
    % Summary of results
    data3.EmpiricalError_Hierarchical_Complete=data2.EmpiricalError;
    %data3.TheoreticalError_Hierarchical_Complete=data2.TheoreticalError;
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
        PredictedLabels=test_RandomClustererV2(GenericExpFolderName,ModelNumber,n)'-1;
        Mask6='Predicted_Labels_RandomOp';
        EmpiricalErrorRandomClusterer(i)=test_ComputeEmpiricalClusteringError(GenericExpFolderName,ModelNumber,Mask6);
        TheoreticalErrorRandomClusterer(i)=test_ComputeTheoreticalClusteringError(PredictedLabels,ValidLabelFcnValues,M,GenericExpFolderName,ModelNumber);
%         MaxRAMUsage=0;%[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,t,'off');
%         %         RAMMemory_RandomClusterer(i)=MaxRAMUsage+RAMMemory_ValidLabelFcnValues+RAMMemory_ComputeProbabilities(i);
%         %         TimeOfProcessing_RandomClusterer(i)=toc(TSTART)+TimeCompProbabilities(i);
%         RAMMemory_RandomClusterer(i)=MaxRAMUsage;
        TimeOfProcessing_RandomClusterer(i)=toc(TSTART);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_RandomClusterer(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_RandomClusterer);
    data2.TotalElapsedTimeOfProcessing=toc(tRnd);
    %data2.RAMMemoryUsedForEachSet=mean(RAMMemory_RandomClusterer);
    data2.EmpiricalError=mean(EmpiricalErrorRandomClusterer);
    %data2.TheoreticalError=mean(TheoreticalErrorRandomClusterer);
    test_SaveFinalReport(GenericExpFolderName,'Results_Random_Clusterer',data2);
    
    % Summary of results
    data3.EmpiricalError_Random=data2.EmpiricalError;
    %data3.TheoreticalError_Random=data2.TheoreticalError;
end
IndicatorRandomClusterer='on';
save(fullfile(GenericExpFolderName,'Backup.mat'));


if strcmp(IndicatorKPodClusterer,'off')
    data2=[];
    % Processing Testing Samples Using K-POD Algorithm
    tKP = tic;
    parfor i=1:NumberTestingSets
        fprintf('Alg. KP: Processing dataset N: %d of %d\r',i,NumberTestingSets);
        TSTART=tic;
       % [perfData,t]=test_StartMeasuringRAMUsage;
        ModelNumber=num2str(i);
        PredictedLabels=test_KPod_withMissing(GenericExpFolderName,ModelNumber)'-1; 
        Mask3='Predicted_Labels_KP';
        EmpiricalErrorKP(i)=test_ComputeEmpiricalClusteringError(GenericExpFolderName,ModelNumber,Mask3);
        TheoreticalErrorKP(i)=test_ComputeTheoreticalClusteringError(PredictedLabels,ValidLabelFcnValues,M,GenericExpFolderName,ModelNumber);
%         MaxRAMUsage=0;%[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,t,'off');
%         %         RAMMemory_KM(i)=MaxRAMUsage+RAMMemory_ValidLabelFcnValues+RAMMemory_ComputeProbabilities(i);
%         %         TimeOfProcessing_KM(i)=toc(TSTART)+TimeCompProbabilities(i);
%         RAMMemory_KP(i)=MaxRAMUsage;
        TimeOfProcessing_KP(i)=toc(TSTART);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_KP(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_KP);
    data2.TotalElapsedTimeOfProcessing=toc(tKP);
    %data2.RAMMemoryUsedForEachSet=mean(RAMMemory_KP);
    data2.EmpiricalError=mean(EmpiricalErrorKP);
    %data2.TheoreticalError=mean(TheoreticalErrorKP);
    test_SaveFinalReport(GenericExpFolderName,'Results_K_Pod',data2);
    
    % Summary of results
    data3.EmpiricalError_KPod=data2.EmpiricalError;
    %data3.TheoreticalError_KPod=data2.TheoreticalError;
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
        TheoreticalErrorFCM_OPS(i)=test_ComputeTheoreticalClusteringError(PredictedLabels,ValidLabelFcnValues,M,GenericExpFolderName,ModelNumber);
%         MaxRAMUsage=0;%[~,MaxRAMUsage]=test_FinishMeasuringRAMUsage(perfData,t,'off');
%         %         RAMMemory_FCM(i)=MaxRAMUsage+RAMMemory_ValidLabelFcnValues+RAMMemory_ComputeProbabilities(i);
%         %         TimeOfProcessing_FCM(i)=toc(TSTART)+TimeCompProbabilities(i);
%         RAMMemory_FCM_OPS(i)=MaxRAMUsage;
        TimeOfProcessing_FCM_OPS(i)=toc(TSTART);
        fprintf('TOTAL ELAPSED TIME: %d\r\r',TimeOfProcessing_FCM(i));
    end
    data2.ElapsedTimeForEachSet=mean(TimeOfProcessing_FCM_OPS);
    data2.TotalElapsedTimeOfProcessing=toc(tFCM_OPS);
    %data2.RAMMemoryUsedForEachSet=mean(RAMMemory_FCM_OPS);
    data2.EmpiricalError=mean(EmpiricalErrorFCM_OPS);
    %data2.TheoreticalError=mean(TheoreticalErrorFCM_OPS);
    test_SaveFinalReport(GenericExpFolderName,'Results_Fuzzy_C_Means_OPS',data2);
    
    % Summary of results
    data3.EmpiricalError_FuzzyCMeans_OPS=data2.EmpiricalError;
    %data3.TheoreticalError_FuzzyCMeans_OPS=data2.TheoreticalError;
end
IndicatorFuzzyCMeansOPSClusterer='on';
save(fullfile(GenericExpFolderName,'Backup.mat'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(IndicatorErrorPlots,'off')
    % Saving Results
    test_SaveFinalReport(GenericExpFolderName,'Empirical and Theoretical Errors',data3);
    
    % Plotting statistics
    figure
    XValues=0:length(Indicator)-1;
    hb=bar(XValues,Indicator);
    ylabel('Frequency');
    xlabel('Maximum Hamming distance of partitions of the space of search from Pmax');
    set(gca,'XTick',XValues);
    bar_child=get(hb,'Children');
    set(bar_child,'CData',Indicator);
    mycolor=[0 0 0;0 0 0;0 0 0];
    colormap(mycolor);
    % Saving graphic files
    PathName1=[GenericExpFolderName '/GRAPHICS'];
    mkdir(PathName1);
    FileName1=fullfile(PathName1,'STATISTICS_Optimal');
    saveas(gcf,[FileName1 'GREY'],'jpeg');
    saveas(gcf,[FileName1 'GREY'],'eps');
    % Saving numeric files
    PathName2=[GenericExpFolderName '/REPORT'];
    mkdir(PathName2);
    FileName2=fullfile(PathName2,'STATISTICS_Optimal.def');
    save(FileName2,'Indicator','-ascii');
    
    % Plotting statistics about the Hamming distance between the predicted and
    % the maximum probability partitions
    figure
    hb=bar(XValues,FreqHammingDistancePredictedPmaxPartitions);
    ylabel('Frequency');
    xlabel('Hamming Distance between the optimal and Pmax partitions');
    set(gca,'XTick',XValues);
    bar_child=get(hb,'Children');
    set(bar_child,'CData',Indicator);
    mycolor=[0 0 0;0 0 0;0 0 0];
    colormap(mycolor);
    % Saving graphic files
    mkdir(PathName1);
    FileName12=fullfile(PathName1,'STATISTICS_HammingDistance_OptimalvsPmax_');
    saveas(gcf,[FileName12 'GREY'],'jpeg');
    saveas(gcf,[FileName12 'GREY'],'eps');
    % Saving numeric files
    FileName22=fullfile(PathName2,'STATISTICS_HammingDistance_OptimalvsPmax_.def');
    save(FileName22,'Indicator','-ascii');
    

    
    %Saving Simulation settings
    test_SaveFinalReport(GenericExpFolderName,'Results_SimulationSettings',data);
    

    
    % Saving empirical errors as a function of the Hamming distance
    EE.EmpiricalError_Optimal=mean(EmpiricalErrorBayesOptClusterer);
    EE.EmpiricalError_SubOptimalPmax=mean(EmpiricalErrorBayesSubOptPmaxClusterer);
    EE.EmpiricalError_SubOptimalPseed=mean(EmpiricalErrorBayesSubOptPseedClusterer);
    EE.EmpiricalError_FCM=mean(EmpiricalErrorFCM);
    EE.EmpiricalError_KM=mean(EmpiricalErrorKM);
    EE.EmpiricalError_HieSi=mean(EmpiricalErrorHierarchicalSi);
    EE.EmpiricalError_HieCo=mean(EmpiricalErrorHierarchicalCo);
    EE.EmpiricalError_Random=mean(EmpiricalErrorRandomClusterer);
    
    EE.EmpiricalError_FCM_OPS=mean(EmpiricalErrorFCM_OPS);
    EE.EmpiricalError_KP=mean(EmpiricalErrorKP);
    
    EE.EmpiricalError_Std_Optimal=std(EmpiricalErrorBayesOptClusterer);
    EE.EmpiricalError_Std_SubOptimalPmax=std(EmpiricalErrorBayesSubOptPmaxClusterer);
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
    test_PlotResults_new(GenericExpFolderName,'Results_EmpiricalError_HammingDistance');
    


end
IndicatorErrorPlots='on';

save(fullfile(GenericExpFolderName,'REPORT','CommandWindowVariables.mat'));
delete(fullfile(GenericExpFolderName,'Backup.mat'));
return