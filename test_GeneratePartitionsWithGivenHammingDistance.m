function [Out,ListOfPossibleHammingDistance]=test_GeneratePartitionsWithGivenHammingDistance(n1,n2,InitialPartition,NumPointsDifferent)
% test_GeneratePartitionsWithGivenHammingDistance generates all vetors of
% labels with n1 1s and n2 0s having a Hamming distance
% 'NumPointsDifferent' from the 'InitialPartition'
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%
n=n1+n2;
if NumPointsDifferent==0
    Out=InitialPartition;
else
    if InitialPartition(1)==1
        InitialPartition=1-InitialPartition;
    end
    
    if mod(n,2)==0
        M=0.5;
    else
        M=(n-1)/(2*n);
    end
    
    if NumPointsDifferent<=M*n
        [ListOfPossibleHammingDistance,...
            ListOfHammingDistance,...
            Num0sModify,~,...
            Num1sModify,~]=test_FindListOfPossibleHammingDistances(InitialPartition);
        IdxOfGivenHammingDistance=find(ListOfHammingDistance==NumPointsDifferent);
        if isempty(IdxOfGivenHammingDistance)
            fprintf('\rInvalid Hamming distance. Possible values are: [%s]\r',num2str(ListOfPossibleHammingDistance));
            Out=[];
            return
        end
        IdxOf0s=find(InitialPartition==0);
        IdxOf1s=find(InitialPartition==1);
        if length(IdxOfGivenHammingDistance)==1
            CombinationsOf0s{1}=combnk(IdxOf0s,Num0sModify(IdxOfGivenHammingDistance));
            CombinationsOf1s{1}=combnk(IdxOf1s,Num1sModify(IdxOfGivenHammingDistance));
        else
            CombinationsOf0s=[];
            CombinationsOf1s=[];
            for i=1:length(IdxOfGivenHammingDistance)
                CombinationsOf0s{i}=combnk(IdxOf0s,Num0sModify(IdxOfGivenHammingDistance(i)));
                CombinationsOf1s{i}=combnk(IdxOf1s,Num1sModify(IdxOfGivenHammingDistance(i)));
            end
        end
        
        NumOfOutputLabelFcns=0;
        for i=1:size(CombinationsOf0s,2)
            NumOfOutputLabelFcns=NumOfOutputLabelFcns+size(CombinationsOf0s{i},1)*size(CombinationsOf1s{i},1);
        end
        
        Out=repmat(InitialPartition,NumOfOutputLabelFcns,1);
        
        cont=0;
        index=true(size(Out,1),1);
        for k=1:length(IdxOfGivenHammingDistance)
            for i=1:size(CombinationsOf0s{k},1)
                for j=1:size(CombinationsOf1s{k},1)
                    cont=cont+1;
                    Out(cont,CombinationsOf0s{k}(i,:))=1;
                    Out(cont,CombinationsOf1s{k}(j,:))=0;
                    if (cont>=2) && (n1==n2)
                        OutCmp=repmat(Out(cont,:),cont-1,1);
                        Sum1s=sum(OutCmp~=Out(1:cont-1,:),2);
                        Sum0s=sum(OutCmp~=(1-Out(1:cont-1,:)),2);
                        Indicator=sum((Sum1s==n)|(Sum0s==n));
                        if Indicator~=0
                            index(cont)=false;
                        end
                    end
                end
            end
        end
        Out=Out(index,:);
    else
        fprintf('\rInvalid Hamming distance. Its value should be in the range: [%d %d]\r',0,M*n);
        Out=[];
        return;
    end
end
Out(Out(:,1)==1,:)=1-Out(Out(:,1)==1,:);
return