function W=test_ComputeMatrixWv2(Rows,Columns,LabelFcns,ValidLabelFcnValuesSortedDescendent_phsi_i1,ValidLabelFcnValuesSortedDescendent_phsi_i2)
% test_ComputeMatrixWv2 computes the mismatch error (cost) between each vector
% of labels from the list 'LabelFcns' and each vector of labels from the list
% 'ValidLabelFcnValuesSortedDescendent_phsi'
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%

% Computing Matrix W
n=size(LabelFcns,2);
row=length(Rows);
W=zeros(row,length(Columns));
for j=1:length(Columns)
    W(:,j)=(1/n)*min(sum(bsxfun(@ne,LabelFcns(Columns(j),:),ValidLabelFcnValuesSortedDescendent_phsi_i1(Rows,:)),2),...
                     sum(bsxfun(@ne,LabelFcns(Columns(j),:),ValidLabelFcnValuesSortedDescendent_phsi_i2(Rows,:)),2));    
end
return