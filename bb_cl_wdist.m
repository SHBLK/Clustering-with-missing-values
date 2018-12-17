function D = bb_cl_wdist(Metric,Data,Weights)
% B_CL_WDIST   Pairwise weighted distance matrix
%
%     D = BB_CL_WDIST(METRIC,DATA,WEIGHTS)
%
%    Input:
%        METRIC   - Number or string. Distance function to be used.
%                   The possible values for METRIC are:
%                   1 -  Euclidean Distance (DEFAULT)
%                   2 -  1-abs(Pearson's correlation)
%                   3 -  1-Pearson's correlation
%                   4 -  1-abs(uncentered correlation)
%                   5 -  1-uncentered correlation
%        DATA     - Matrix NxM with the numeric data to be clustered.
%                   Each column is a vector and each row a variable
%        WEIGHTS  - Matrix NxM with the weight for each value of DATA.
%                   The values of this matrix must be betwee 0 and 1.
%
%    Output:
%        D        - MxM Matrix. Pairwise distance between the vectors
%                   in the DATA matrix.
%
%    D = bb_cl_wdist(A,METRIC) returns in D a simetric matrix with
%    the pairwise distance of the vectors listed in A, with the metric
%    defined in METRIC, using the weights defined in WEIGHTS.
%     Each column in A is a vector, and the dimension 
%    is defined by the number of rows in A.
%
%    Observation: In this case the vectors must be in the columns of
%    the data matrix. In case of using rows for the vectors, the 
%    matrix data must be transposed before the call to this routine.
%
%
%    Examples
%    --------
%     
%       A = [1 2 3 ; 1 2 3];
%       Q = [1 1 1 ; 1 1 1];
%       bb_cl_wdist(1,A,Q)
%
%       A = [1 1 ; 2 2 ; 3 3];
%       Q = [1 1 ; 1 1 ; 1 1];
%       bb_cl_wdist(1,A',Q')
%

%   C-Mex function
%   By Marcel Brun, August 1 2001
%   mbrun@vision.ime.usp.br


if Metric==1
   D = pdist(Data','euclidean')./sqrt(size(Data,1));
elseif Metric==3
   D = pdist(Data','correlation');
else
   error('Metric not defined yet');
end
D = D';

