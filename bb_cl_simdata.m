function [DATA,C] = bb_cl_simdata(Centers,Vars,Cants,Repl,Distr)
% BB_CL_SIMDATA   Generate data for tests
%
%    [DATA,C] = BB_CL_SIMDATA(CENTERS,VARS,CANTS,REPL,DISTR)
%
%    Input:
%        CENTERS - NxM Matrix. Templates for the classes. Each row is a 
%                  different template, and the number of columns
%                  is the dimension of the data.
%        VARS    - Scalar, 1xM Array or NxM Matrix. 
%                  a)  If scalar, all the covariance matrices will be 
%                      diagonal with VARS in all the non null elements. 
%                  b)  If 1xM array, all covariances matrices will be
%                      diagonal, with VARS(i) in the element (i,i)
%                  c)  If MxN Matrix, each row indicates the N elements
%                      in the diagonal for the covariance matrix for 
%                      the associated template 
%                  d)  If cell, each element has to be the covariance matrix
%                      When using laplacian noise, will be used only
%                      the diagonal of the matrices.
%        CANTS   - Scalar or 1xN array
%                  a)  If scalar, will be created CANTS representants 
%                      for each template
%                  b)  If 1xN array, will be created CANTS(i) representants
%                      for the template "i"
%        REPL    - Replications. Each representant will be generated REPL
%                  times and the result will be the mean value. DEFAULT:1
%        DISTR   - Parametrical distribution to use. DEFAULT:'ga' (gaussian)
%                  The possible values for DISTR are
%                      'ga' -  Gaussian distribution (default)
%                      'la' -  Laplacian distribution 
%
%    Output:
%        DATA    - Matrix. Vectors simulated from templates. The number of
%                  columns is M, and the number of rows depends on N and
%                  the values in CANTS. Each row is a vector and each column 
%                  a variable. The vectors are placed randomly in the matrix
%        C       - Array. Class labels for each vector. If the vector in 
%                  position "n" came from template "i", then C(n)=i .
%
%    [DATA,C] = bb_cl_simdata(CENTERS,VARS,CANTS) generates simulation with
%    CANTS element for each class. The templates are listed in the matrix 
%    CENTERS, where the each row is a template and the number of columns is
%    the dimension of the space.The covariance matrix will be diagonal with VARS
%    in the diagonal. If CANTS is a vector with the same size that the number of 
%    templates, then CANTS(i) indicates the number of elements for the template i.
%    If VARS is a vector with the same number than the dimension, then VARS(i)
%    indicates the value for the element (i,i) of the covariance matrix.
%    The data simulated is in the matrix DATA, and the vector C has the class
%    identifier for each row of the data.
%
%    [DATA,C] = bb_cl_simdata(...,REPL,DISTR) generates simulation with REPL
%    replications, the values are averaged over the REPL replicates. DISTR 
%    indicates the distribution to use.
%
%    Examples
%
%    Centers = [ 1 4 7 ; 7 4 1 ];
%    Vars    = [ 1 1 1 ; 1 1 1 ];
%    Cants   = [ 10 10 ];
%    [Data,Cl] = bb_cl_simdata(Centers, Vars, Cants);
%    C = bb_cl_kmeans(Data,2);
%    bb_cl_misclass(C,Cl)
%
%    Centers = [ 1 4 7 ; 7 4 1 ];
%    Vars = {[ 1 0 0 ; 0 1 0 ; 0 0 1 ];[ 1 0 0 ; 0 1 0 ; 0 0 1]}
%    Cants   = [ 10 10 ];
%    [Data,Cl] = bb_cl_simdata(Centers, Vars, Cants);
%    C = bb_cl_kmeans(Data,2);
%    bb_cl_misclass(C,Cl)
%
% See also BB_CL_CMEANS, BB_CL_KMEANS, BB_CL_SOM, BB_CL_HIERARCHICAL, BB_CL_PROFILE
%
  
% Changes may 1 2001
%     a) Now accepts a complete covariance matrix for each template
% 

if nargin<5
   Distr = 'gauss';
end
if nargin<4
   Repl = 1;
end
if nargin<3 
   Cants = ones(1,size(Centers,1));
end
if nargin<2
   Vars = ones(size(Centers));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     If Vars is scalar     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscell(Vars) & size(Vars,1)==1 & size(Vars,2)==1
   Vars = ones(size(Centers)).*Vars;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     If Vars is an array     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscell(Vars) & size(Vars,1)==1 & size(Vars,2)==size(Centers,2)
   Vars = ones(size(Centers,1),1)*Vars;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     If Cants is scalar     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(Cants,1)==1 & size(Cants,2)==1
   Cants = ones(1,size(Centers,1)).*Cants;
end


% Dimension of data vectors
DimeData = size(Centers,2);

% Number of classes
CantClas = size(Centers,1);   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     Covariance Matrix     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~iscell(Vars)
   for i=1:CantClas
      Cov{i}= diag(Vars(i,:));
   end
else
   %if size(Vars,1)>1
   %   Vars = Vars';
   %end
   Cov = Vars;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     Assign classes to vectors     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Inic = 0;
for k=1:CantClas
   C1(1,Inic+1:Inic+Cants(k))=ones(1,Cants(k))*k;
   Inic = Inic+Cants(k);
end

TotalCant = Inic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     mixtures the classes     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k=1:TotalCant
   Pos = unidrnd(TotalCant-k+1,1,1);
   C(k)=C1(Pos);
   C1(Pos)=C1(TotalCant-k+1);
   C1(TotalCant-k+1)=C(k);
end


%------------------------------------------%
%   Creates samples and take means values  %
%------------------------------------------%
DATA = zeros(TotalCant,DimeData);
Distr = [Distr '  '];
Distr = Distr(1,1:2);
Distr = lower(Distr);


if strcmp(Distr,'ga')==1

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%     Gaussian Distribution     %%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if Repl==1
      for k=1:TotalCant
         %DATA(k,1:DimeData) = mvnrnd(Centers(C(k),:),diag(Vars(C(k),:)),1);
         DATA(k,1:DimeData) = mvnrnd(Centers(C(k),:),Cov{C(k)},1);
      end
   else
      for k=1:TotalCant
         %DATA(k,1:DimeData) = sum(mvnrnd(Centers(C(k),:),diag(Vars(C(k),:)),Repl))/Repl;
         DATA(k,1:DimeData) = sum(mvnrnd(Centers(C(k),:),Cov{C(k)},Repl))/Repl;
      end
   end

elseif strcmp(Distr,'la')==1

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%     Laplacian Distribution     %%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if Repl==1
      for k=1:TotalCant
         %DATA(k,1:DimeData) = bb_cl_laplacian(Centers(C(k),:),Vars(C(k),:),1);
         DATA(k,1:DimeData) = bb_cl_laplacian(Centers(C(k),:),diag(Cov{C(k)}),1);
      end
   else
      for k=1:TotalCant
         %DATA(k,1:DimeData) = sum(bb_cl_laplacian(Centers(C(k),:),Vars(C(k),:),Repl))/Repl;
         DATA(k,1:DimeData) = sum(bb_cl_laplacian(Centers(C(k),:),diag(Cov{C(k)}),Repl))/Repl;
      end
   end
end