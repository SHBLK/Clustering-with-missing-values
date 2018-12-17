function x = bb_cl_laplacian(Mean,Var,Cases);
% BB_CL_LAPLACIAN   Generates random number with laplacian distribution
%
%    DATA = BB_CL_LAPLACIAN(MEAN,VARS,NUM);
%
%    Input:
%        MEAN  - 1xN Array. Array with mean values
%        VARS  - 1xN Array. Array with variances
%        NUM   - Number of vectors to generate
%    Output:
%        DATA  - NUMxN Matrix. Matrix with NUM rows and
%                N columns.
%
%    DATA = BB_CL_LAPLACIAN(MEAN,VARS,NUM) generates NUM vectors
%    with Laplacian distribution with parameters MEAN(i) and
%    VARS(i) for each variable "i" in the vector
%
%    Examples
%    --------
%
%    Data = gen_laplacian([1 10],[1 1],5)
%
%    See also BB_CL_SIMDATA

%=================================================================
% The generation usses the next algorithm:
%
% y = unifrnd(0,1);
% if y<0.5 
%    x = ((1/b)*log(2*x)+a);
% else
%    x = ((-1/b)*log(2-2*x)+a);
% end
% The parameters (a,b) are computed from MEAN and VARS
%=================================================================
  
if size(Mean,1)>1
   Mean = Mean';
end
if size(Var,1)>1
   Var = Var';
end
 
a = ones(Cases,1)*Mean;
b = ones(Cases,1)*sqrt(2./Var);

y = unifrnd(zeros(Cases,size(a,2)),ones(Cases,size(a,2)));

x = ((1./b).*log(2*y)+a).*(y<0.5) + ((-1./b).*log(2-2*y)+a).*(y>=0.5);
