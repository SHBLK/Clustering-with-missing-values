function Z = bb_cl_linkage(D,Method)
% BB_CL_LINKAGE   Linkage from pairwise distance matrix
%
%    L = BB_CL_LINKAGE(D,METHOD) 
%
%    Input:
%        D       - MxM Matrix. Pairwise distance between the vectors
%                  in the DATA matrix.
%        METHOD  - Number. Algorithm used to compute the linkage
%                  The possible values for METHOD are:
%                  1 -  Complete (furthest distance, DEFAULT)
%                  2 -  Single (nearest distance)
%                  2 -  Average (average distance)%
%
%    Output:
%        LINKAGE - 3xN Matrix. Linkage values obtained from the pairwise 
%                  distances matrix, used to obtain the clusters.
%                  
%
%    L = BB_CL_LINKAGE(D,METHOD) Returns in L the linkage
%    information from the pairwise distance matrix D, using
%    the aglomerative method defined by METHOD. The linkage
%    information in L is used in bb_cl_hierclus to define 
%    the clusters and in bb_cl_dend to design the dendrogram
%
%    Examples
%    --------
%     
%       A = [1 1 ; 2 2 ; 3 3;];
%       D = bb_cl_distance(A',1);
%       L = bb_cl_linkage(D,1);
%       [O,Im] = bb_cl_dend(L,5,30);
%       image(Im);
%
% See alseo BB_CL_DISTANCE, BB_CL_DEND, BB_CL_HIERARCHICAL


if Method==1
   Z = linkage(D,'complete');
elseif Method==2
   Z = linkage(D,'single');
elseif Method==3
   Z = linkage(D,'average');
else
   error('Linkage method not available');
end
