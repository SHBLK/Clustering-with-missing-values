function IdxOfLabelFcnsOut=test_FindLabelFcnsHammingDistance(LabelFcnsIn,LabelFcnTarget,distance,varargin)
% test_FindLabelFcnsHammingDistance finds the indices of the partitions
% in the list 'LabelFcnsIn' that have a given Hamming distance 'distance'
% from the partition 'LabelFcnTarget'
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%
if nargin==4
    try
        aux=varargin{1};
    end
    if isempty(aux)==1
        varargin{1}='Equal';
        varargin{2}='on';
    else
        if iscell(aux)==1
            varargin=varargin{1};
        else
            varargin{1}=aux;
            varargin{2}='on';
        end
    end
end
if nargin==3
    varargin{1}=[];
    varargin{2}='on';
end
if strcmp(varargin{2},'on')
    Tstart=tic;
end
[row,~]=size(LabelFcnsIn);
ListOfIndices=(1:row)';
% Matrix where each row is the target label function
LabelFcnTargetMatrix=repmat(LabelFcnTarget,row,1);
if nargin==3
IdxOfLabelFcnsOut=ListOfIndices(min(sum(LabelFcnsIn~=LabelFcnTargetMatrix,2),sum(LabelFcnsIn~=(1-LabelFcnTargetMatrix),2))==distance);
else
    if strcmp(varargin{1},'LesserOrEqual')==1
        IdxOfLabelFcnsOut=ListOfIndices(min(sum(LabelFcnsIn~=LabelFcnTargetMatrix,2),sum(LabelFcnsIn~=(1-LabelFcnTargetMatrix),2))<=distance);
    elseif strcmp(varargin{1},'Lesser')==1
        IdxOfLabelFcnsOut=ListOfIndices(min(sum(LabelFcnsIn~=LabelFcnTargetMatrix,2),sum(LabelFcnsIn~=(1-LabelFcnTargetMatrix),2))<distance);
    elseif strcmp(varargin{1},'GreaterOrEqual')==1
        IdxOfLabelFcnsOut=ListOfIndices(min(sum(LabelFcnsIn~=LabelFcnTargetMatrix,2),sum(LabelFcnsIn~=(1-LabelFcnTargetMatrix),2))>=distance);
    elseif strcmp(varargin{1},'Greater')==1
        IdxOfLabelFcnsOut=ListOfIndices(min(sum(LabelFcnsIn~=LabelFcnTargetMatrix,2),sum(LabelFcnsIn~=(1-LabelFcnTargetMatrix),2))>distance);
    elseif strcmp(varargin{1},'Equal')==1
        IdxOfLabelFcnsOut=ListOfIndices(min(sum(LabelFcnsIn~=LabelFcnTargetMatrix,2),sum(LabelFcnsIn~=(1-LabelFcnTargetMatrix),2))==distance);
    else
        fprintf('ERROR: The option %s is incorrect\r',varargin{1});
    end
end
if strcmp(varargin{2},'on')
    fprintf('Computing the Hamming Distance. Elapsed time: %d s\r',toc(Tstart));
end
return