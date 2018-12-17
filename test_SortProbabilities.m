function [ValidLabelFcnValuesSortedDescendent_phsi_i1,...
          ValidLabelFcnValuesSortedDescendent_phsi_i2,...
          LabelFcnValues,...
          ProbabiliesPartitionsSortedDescendent]=test_SortProbabilities(IDX,...
          ValidLabelFcnValues,AllLabelFcnValues,...
          ProbabiliesPartitionsSortedDescendent,n1,n2)  
% test_SortProbabilities sorts in descendent order the probabilities of the 
% vector of labels that compose the set of reference partitions 
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%
ValidLabelFcnValuesSortedDescendent_phsi_i1=ValidLabelFcnValues.phsi_i1(IDX,:);
ValidLabelFcnValuesSortedDescendent_phsi_i2=ValidLabelFcnValues.phsi_i2(IDX,:);

% Finding Impossible partitions
Numb1sPerRow=sum(AllLabelFcnValues.phsi,2);
IdxOfImpossibleLabelFcnValues=logical((1-double(Numb1sPerRow==n1)).*(1-double(Numb1sPerRow==n2)));
ImpossibleLabelFcnValues=AllLabelFcnValues.phsi(IdxOfImpossibleLabelFcnValues,:);
LabelFcnValues=[ValidLabelFcnValuesSortedDescendent_phsi_i1;ImpossibleLabelFcnValues];
return