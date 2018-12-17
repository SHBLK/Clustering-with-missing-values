function [Cl,CEN,LINK,Dist] = bb_cl_hierarchical(DATA,NUMCLUS,DISTANCE,LINKAGE,WEIGHTS)
% BB_CL_HIERARCHICAL   Clustering using Hierarchical Clustering algorithm 
%
%    [CL,CENTERS,LINKAGE] = BB_CL_HIERARCHICAL(DATA,NUMCLUS,DISTANCE,METHOD,WEIGHTS)
%
%    Input:
%        DATA     - Matrix NxM with the numeric data to be clustered.
%                   Each row is a vector and each column a variable
%        NUMCLUS  - Positive integer. Number of clusters to be created
%        DISTANCE - Number or string. Distance function to be used.
%                   The possible values for DISTANCE are:
%                   'eu' or 1 -  Euclidean Distance (DEFAULT)
%                   'ce' or 2 -  1-abs(centered correlation)
%                   'c2' or 3 -  1-centered correlation
%                   'un' or 4 -  1-abs(uncentered correlation)
%                   'u2' or 5 -  1-uncentered correlation
%        METHOD  - Number or string. Algorithm used to compute the linkage
%                   The possible values for METHOD are:
%                   'co' or 1 -  Complete (furthest distance, DEFAULT)
%                   'si' or 2 -  Single (nearest distance)
%                   'av' or 3 -  Average (average distance)
%        WEIGHTS  - Matrix NxM with the weight for each value of DATA.
%                   The values of this matrix must be betwee 0 and 1.
%                  (DEFAULT = empty array - no weights)
%
%    Output:
%        CL      - Array 1xN with the cluster label for each vector
%                  Label values go from 1 to NUMCLUS
%        CENTERS - Array NUMCLUSxN with the centers of the clusters
%                  For empty clusters, the value is the zero vector
%        LINKAGE - Nx3 Matrix. Linkage values obtained from the pairwise 
%                  distances matrix, used to obtain the clusters.
%                  
%
%    BB_CL_HIERARCHICAL finds NUMCLUS clusters for the data in DATA 
%    using Hierarchical clustering algorithms, from the pairwise distance
%    matrix betwee the vectors in DATA, computed using the distance defined
%    by DISTANCE and the linkage method defined in METHOD. The matrix 
%    returned in LINK is used to compute the dendrogram.
%
%    Examples
%    --------
%
%    Data = [-1 2 3 ; 2 -3 1 ; -2 3 2 ; -3 3 -3 ; 2 -2 1; -3 2 1];
%    [C,CEN,LINK] = bb_cl_hierarchical(Data,3,'eu','co')
%    [Order,ImDend] = bb_cl_dend(LINK,5,30,C,3);
%    Data2 = Data(Order,:) % Reorder data to fit dendrogram
%    Im = bb_cl_profile(Data2,5,15);                
%    subplot(1,2,1);image(Im);subplot(1,2,2);image(ImDend);
%
% See also BB_CL_CMEANS, BB_CL_SOM, BB_CL_KMEANS, BB_CL_DEND
%

%    Changes December 11 200 - By Marcel Brun
%       a) Now the output parameters are [Cl,CEN,LINK] to have 
%          coherence with the other algorithms.

%    Changes August 1 2001 - By Marcel Brun
%       a) Now it accepts Quality matrix
%       b) Now it calls BB_CL_WDIST in place of BB_CL_DISTANCE
%       c) Now it computes the centers ussing the weight

%
%    Changes October 10 2002 - By Marcel Brun
%       a) Default value empty for quality matrix

if nargin<5
   %WEIGHTS = ones(size(DATA));
   WEIGHTS  = '';
end
if nargin<4
   LINKAGE = 'co';
end
if nargin<3
   DISTANCE = 'eu';
end

if isempty(WEIGHTS)
   WEIGHTS = ones(size(DATA));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%     Convert string to numeric distance     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isnumeric(DISTANCE)
   DISTANCE = bb_cl_findstring(DISTANCE,{'eu','ce','c2','un','u2'});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%     Convert string to numeric linkage     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isnumeric(LINKAGE)
   LINKAGE = bb_cl_findstring(LINKAGE,{'co','si','av'});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     Classification     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Dist    = bb_cl_wdist(DISTANCE,DATA',double(WEIGHTS'));
LINK    = bb_cl_linkage(Dist',LINKAGE)';
Cl      = bb_cl_hierclus(LINK',NUMCLUS);

if min(Cl)==0
   Cl = Cl+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     Centers     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
CEN = [];
Dim = size(DATA,2);
for k=1:NUMCLUS
   Qual = sum(WEIGHTS(find(Cl==k),:),1);
   Sum  = sum(DATA(find(Cl==k),:).*WEIGHTS(find(Cl==k),:),1);
   for i=1:Dim
       if Qual(i)==0
          CEN(k,i)=10^20;
       else
          CEN(k,i)=Sum(i)./Qual(i);
       end
   end
   %if sum(Cl'==k)>0
   %   CEN(k,1:Dim)= sum(DATA.*( (Cl'==k)*ones(1,Dim) ) )/sum(Cl'==k);
   %else
   %   CEN(k,1:Dim)= zeros(1,Dim);
   %end
end

return