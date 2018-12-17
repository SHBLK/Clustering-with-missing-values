function CL = bb_cl_hierclus(Link,NumClus)
% BB_CL_HIERCLUS   Compute clusters from linkage information
%
%    CL = BB_CL_HIERCLUS(LINK,NUMCLUS);
%
%    Input:
%        LINK    - 3xN Matrix. Linkage values obtained from the pairwise 
%                  distances matrix
%        NUMCLUS - Positive integer. Number of clusters to be created
%
%    Output:
%        CL      - Array 1xN with the cluster label for each vector
%                  Label values go from 1 to NUMCLUS
%                  
%    BB_CL_HIERCLUS finds NUMCLUS clusters from the linkage information
%    in LINK, obtained from the pairwise distance matrix.
%
%    Examples
%    --------
%     
%       A = [1 1 ; 2 1 ; 4 3; 4 4 ;5 4];
%       D = bb_cl_distance(A',1);
%       L = bb_cl_linkage(D,1);
%       CL = bb_cl_hierclus(L,2)
%
% See alseo BB_CL_DISTANCE, BB_CL_DEND, BB_CL_HIERARCHICAL, BB_CL_LINKAGE



CL = cluster(Link,'MAXCLUST',NumClus)';