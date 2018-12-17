function Pos = bb_cl_findstring(String,List)
% BB_CL_FINDSTRING   return poistion of the string in the list 
%
%    POS = BB_CL_FINDSTRING(STRING,LIST) 
%
%    Input:
%        STRING   - String. String to search
%        LIST     - Cell Array. List with strings
%
%    Output:
%        POS      - Position of the string in the ist
%
%    POS = BB_CL_FINDSTRING(STRING,LIST) returns in POS the 
%    position of the string STRING in the cell list of strings
%    LIST. It only compares the two first characters of the 
%    string in lower case.
%    This routine is useful when we want to have parameters
%    passed like string to routines that uses numeric parameters
%
%    Examples
%    --------
%
%    bb_cl_findstring('single',{'complete','single','average'})
%
%    bb_cl_findstring('single',{'co','si','av'})
%
%    bb_cl_findstring('si',{'complete','single','average'})
%
% See also BB_CL_HIERARCHICAL

String = [String '  '];
String = lower(String(1:2));
Pos = 0;
for i=1:size(List,2)
   String2 = [List{i} '  '];
   String2 = lower(String2(1:2));
   if strcmp(String,String2)
      Pos = i;
   end
end
