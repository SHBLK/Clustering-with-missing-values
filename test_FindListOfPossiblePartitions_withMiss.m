function [ListOfPossiblePartitions,ListOfValidHammingDistances,ListTimeOfProcessing,ListRAMUsage,tRankingProbabilities]=test_FindListOfPossiblePartitions_withMiss(TypeModel,Vectors,HypParam,SeedLabelFcn,MaximumHammingDistance,MissMatr)
%Flipping the points of the vector of labels and computing the
%probabilities
n1=HypParam.quantities(1,1);
n2=HypParam.quantities(1,2);
n=n1+n2;

tStart1=tic;
for jj=1:length(SeedLabelFcn)
    TemporalPartition.phsi_i1=SeedLabelFcn;
    TemporalPartition.phsi_i1(jj)=1-TemporalPartition.phsi_i1(jj);
    [val1,val2]=test_ComputeProbabilityOfPartitions_withMiss(TypeModel,Vectors,HypParam,TemporalPartition,MissMatr,'NoNormalize');
    Prphi1_S(jj)=val1;
    Prphi2_S(jj)=val2;
    clear TemporalPartition;
end
ProbabilitiesTemporalPartitions=Prphi1_S+Prphi2_S;
[~,IdxSortedProbabilitiesTemporalPartitions]=sort(ProbabilitiesTemporalPartitions,'descend');
RankedVectorOfPoints=SeedLabelFcn(IdxSortedProbabilitiesTemporalPartitions);

[ListOfPossibleHammingDistance,...
    ListOfHammingDistance,...
    Num0sModify,~,...
    Num1sModify,~]=test_FindListOfPossibleHammingDistances(SeedLabelFcn);
tRankingProbabilities=toc(tStart1);

% Finding the maximum number of points to change in each class for having
% the given hamming distance
Idxs1=ListOfHammingDistance<=MaximumHammingDistance;
ListOfNum0sModify1=Num0sModify(Idxs1);
ListOfNum1sModify1=Num1sModify(Idxs1);
param=15;  %Maximum number of points to modify for each class
Idxs2=(ListOfNum0sModify1<=param)&(ListOfNum1sModify1<=param);
ListOfNum0sModify2=ListOfNum0sModify1(Idxs2);
ListOfNum1sModify2=ListOfNum1sModify1(Idxs2);
MaxNum0sModify=max(ListOfNum0sModify2);
MaxNum1sModify=max(ListOfNum1sModify2);

cont=1;
Num0s=1;
Num1s=0;
IdxOf0s=IdxSortedProbabilitiesTemporalPartitions(1);
IdxOf1s=[];
while ~(Num0s>=MaxNum0sModify && Num1s>=MaxNum1sModify) && cont<length(SeedLabelFcn)
    cont=cont+1;
    if RankedVectorOfPoints(cont)==RankedVectorOfPoints(1)
        Num0s=Num0s+1;
        IdxOf0s=[IdxOf0s IdxSortedProbabilitiesTemporalPartitions(cont)];
    else
        Num1s=Num1s+1;
        IdxOf1s=[IdxOf1s IdxSortedProbabilitiesTemporalPartitions(cont)];
    end
end

ListOfCandidatePartitions=[];
ListOfValidHammingDistances=ListOfPossibleHammingDistance(ListOfPossibleHammingDistance<=MaximumHammingDistance);

for tt=1:length(ListOfValidHammingDistances)
    if ListOfValidHammingDistances(tt)==0
        tStart3=tic;
        tCombinationsOfIdx(tt)=0;
        Out=SeedLabelFcn;
    else
        if SeedLabelFcn(1)==1
            SeedLabelFcn=1-SeedLabelFcn;
        end
        IdxOfGivenHammingDistance=find((ListOfHammingDistance==ListOfValidHammingDistances(tt))&((Num0sModify<=param)&(Num1sModify<=param)));
        if isempty(IdxOfGivenHammingDistance)
            tStart3=tic;
            tCombinationsOfIdx(tt)=0;
            Out=[];
        else
            tStart2=tic;
            if length(IdxOfGivenHammingDistance)==1
                CombinationsOf0s{1}=combnk(IdxOf0s,Num0sModify(IdxOfGivenHammingDistance));
                CombinationsOf1s{1}=combnk(IdxOf1s,Num1sModify(IdxOfGivenHammingDistance));
                aux0s=CombinationsOf0s{1};
                aux1s=CombinationsOf1s{1};
                InfoCombinationsOf0s=whos('aux0s');
                RAMMemory_Comb0s(tt)=InfoCombinationsOf0s.bytes/(10^6);
                InfoCombinationsOf1s=whos('aux1s');
                RAMMemory_Comb1s(tt)=InfoCombinationsOf1s.bytes/(10^6);
            else
                CombinationsOf0s=[];
                CombinationsOf1s=[];
                RAMMemory_Comb0s(tt)=0;
                RAMMemory_Comb1s(tt)=0;
                for i=1:length(IdxOfGivenHammingDistance)
                    CombinationsOf0s{i}=combnk(IdxOf0s,Num0sModify(IdxOfGivenHammingDistance(i)));
                    CombinationsOf1s{i}=combnk(IdxOf1s,Num1sModify(IdxOfGivenHammingDistance(i)));
                    aux0s=CombinationsOf0s{i};
                    aux1s=CombinationsOf1s{i};
                    InfoCombinationsOf0s=whos('aux0s');
                    RAMMemory_Comb0s(tt)=RAMMemory_Comb0s(tt)+InfoCombinationsOf0s.bytes/(10^6);
                    InfoCombinationsOf1s=whos('aux1s');
                    RAMMemory_Comb1s(tt)=RAMMemory_Comb1s(tt)+InfoCombinationsOf1s.bytes/(10^6);
                end
            end
            tCombinationsOfIdx(tt)=toc(tStart2);
            
            tStart3=tic;
            NumOfOutputSeedLabelFcns=0;
            for i=1:size(CombinationsOf0s,2)
                NumOfOutputSeedLabelFcns=NumOfOutputSeedLabelFcns+size(CombinationsOf0s{i},1)*size(CombinationsOf1s{i},1);
            end
            
            Out=repmat(SeedLabelFcn,NumOfOutputSeedLabelFcns,1);
            
            count=0;
            index=true(size(Out,1),1);
            for k=1:length(IdxOfGivenHammingDistance)
                for i=1:size(CombinationsOf0s{k},1)
                    for j=1:size(CombinationsOf1s{k},1)
                        count=count+1;
                        Out(count,[CombinationsOf0s{k}(i,:) CombinationsOf1s{k}(j,:)])=1-Out(count,[CombinationsOf0s{k}(i,:) CombinationsOf1s{k}(j,:)]);
                        if (count>=2) && (n1==n2)
                            Sum1s=sum(bsxfun(@ne,Out(count,:),Out(1:count-1,:)),2);
                            Sum0s=sum(bsxfun(@ne,Out(count,:),1-Out(1:count-1,:)),2);
                            Indicator=sum((Sum1s==n)|(Sum0s==n));
                            if Indicator~=0
                                index(count)=false;
                            end
                        end
                    end
                end
            end
            Out=Out(index,:);
            Out(Out(:,1)==1,:)=1-Out(Out(:,1)==1,:);
        end
    end
    ListOfPossiblePartitions{tt}=Out;
    tListOfPossiblePartitions(tt)=toc(tStart3);
    InfoListOfPossiblePartitions=whos('Out');
    RAMMemory_ListOfPossiblePartitions(tt)=InfoListOfPossiblePartitions.bytes/(10^6);
end
ListTimeOfProcessing=tCombinationsOfIdx+tListOfPossiblePartitions;
ListRAMUsage=RAMMemory_Comb0s+RAMMemory_Comb1s+RAMMemory_ListOfPossiblePartitions;
return