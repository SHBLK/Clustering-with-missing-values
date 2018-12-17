function Out=test_GenerateRandomPossiblePartition(n1,n2)
% test_GenerateRandomPossiblePartition generates a random vector of labels
% with n1 1s and n2 0s.
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%

n=n1+n2;
vect=[n1,n2];
RandIdx=randi(2,1,1);
Out=randerr(1,n,vect(RandIdx));
if Out(1)==1
    Out=1-Out;
end
return