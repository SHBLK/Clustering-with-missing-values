function [ListOfPossibleHammingDistance,...
    ListOfHammingDistance,...
    Num0sModify,Num0sKeep,...
    Num1sModify,Num1sKeep]=test_FindListOfPossibleHammingDistances(InitialPartition)
% test_FindListOfPossibleHammingDistances finds the number of labels to be
% flipped in a vector to generate partitions with a given Hamming distance
% from the partition 'InitialPartition'
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%

Num1sSeed=sum(InitialPartition);
Num0sSeed=sum(1-InitialPartition);

% Computing a list of possible combinations
if Num1sSeed<=Num0sSeed
    Num1sKeep=0:Num1sSeed;
    Num1sModify=Num1sSeed-Num1sKeep;
    Num0sKeep=Num0sSeed-Num1sModify;
    Num0sModify=Num0sSeed-Num0sKeep;
else
    Num0sKeep=0:Num0sSeed;
    Num0sModify=Num0sSeed-Num0sKeep;
    Num1sKeep=Num1sSeed-Num0sModify;
    Num1sModify=Num1sSeed-Num1sKeep;
end
% This part of the code is used just for debugging
% %         disp('#0s-Mod  #1s-Mod #0s-Keep #1s-Kepp')
% %         disp([Num0sModify' Num1sModify' Num0sKeep' Num1sKeep'])

%Computing the appropiate combination for the given Hamming distance
ListOfHammingDistance=min(Num1sKeep+Num0sKeep,Num1sModify+Num0sModify);
ListOfPossibleHammingDistance=unique(ListOfHammingDistance);
return