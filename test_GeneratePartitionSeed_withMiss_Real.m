function [BestPartitionSeed,BestProbabilityPartitionSeed,mp_]=test_GeneratePartitionSeed_withMiss_Real(TypeModel,Vectors,HypParam,n1,n2,HammingDistance,MissMatr,varargin)
% test_GeneratePartitionSeed_withMiss_Real finds the maximum probability partition for a
% given point set drawn from a Model whose hyperparameters are 'HypParam'
%

mp_=0;

n=n1+n2;

if nargin==7
    NumIterations=1;
else
    NumIterations=varargin{1};
end

for i=1:NumIterations
InitialPartition.phsi_i1=test_GenerateRandomPossiblePartition(n1,n2);

%Computation of probabilities
if i==1
    [Prphi1_S,Prphi2_S,mp_]=test_ComputeProbabilityOfPartitions_withMiss_Real(TypeModel,Vectors,HypParam,InitialPartition,MissMatr,0,'NoNormalize');
else
    [Prphi1_S,Prphi2_S,~]=test_ComputeProbabilityOfPartitions_withMiss_Real(TypeModel,Vectors,HypParam,InitialPartition,MissMatr,mp_,'NoNormalize');
end

%Initilizing the while loop
MaxProbability(1)=0;
IdxMaxProbability(1)=0;
MaxProbabilityPartition(1,:)=NaN(1,n);
MaxProbability(2)=Prphi1_S+Prphi2_S;
IdxMaxProbability(2)=1;
MaxProbabilityPartition(2,:)=InitialPartition.phsi_i1;
cont=2;

if MaxProbability(cont)~=1
    while MaxProbability(cont)>MaxProbability(cont-1)
        cont=cont+1;
        [ListOfPartitions.phsi_i1,PossibleHammingDistances]=test_GeneratePartitionsWithGivenHammingDistance(n1,n2,MaxProbabilityPartition(cont-1,:),HammingDistance);
        if isempty(ListOfPartitions.phsi_i1)
            fprintf('\r\rWARNING: Hamming distance to find the partition seed changed from %d to %d\r\r',HammingDistance,PossibleHammingDistances(2));
            HammingDistance=PossibleHammingDistances(2);
            [ListOfPartitions.phsi_i1,PossibleHammingDistances]=test_GeneratePartitionsWithGivenHammingDistance(n1,n2,MaxProbabilityPartition(cont-1,:),HammingDistance);
        end
        %[Prphi1_S,Prphi2_S]=test_ComputeProbabilityOfPartitions_withMiss(TypeModel,Vectors,HypParam,ListOfPartitions,MissMatr,'NoNormalize');
        [Prphi1_S,Prphi2_S,~]=test_ComputeProbabilityOfPartitions_withMiss_Real(TypeModel,Vectors,HypParam,ListOfPartitions,MissMatr,mp_,'NoNormalize');
        [val,idx]=max(Prphi1_S+Prphi2_S);
        MaxProbability(cont)=val;
        IdxMaxProbability(cont)=idx;
        MaxProbabilityPartition(cont,:)=ListOfPartitions.phsi_i1(idx,:);
    end
    PartitionSeed(i,:)=MaxProbabilityPartition(cont-1,:);
    ProbabilityPartitionSeed(i)=MaxProbability(cont-1);
else
    PartitionSeed(i,:)=MaxProbabilityPartition(cont,:);
    ProbabilityPartitionSeed(i)=MaxProbability(cont);    
end
end
[BestProbabilityPartitionSeed,IdxBestPartitionSeed]=max(ProbabilityPartitionSeed);
BestPartitionSeed=PartitionSeed(IdxBestPartitionSeed,:);
if any(isnan(BestPartitionSeed))
   disp('Debugging') 
end
return