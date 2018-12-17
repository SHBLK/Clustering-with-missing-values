function [Info,arch] = bb_cl_loadparam(File)
% BB_CL_LOADPARAM   Load information from specification file
%
%
%    [Info] = bb_cl_loadparam(File)      - 
%    [Info] = bb_cl_loadparam()          - Ask File
%    [Info,File] = bb_cl_loadparam(..)   - returns name of the file
%
%    Input:
%        File    - Text file while are located the lines whit information
%
%   Output:
%        Info    - Structure with the information
%        File    - Loaded text file
%
% 
%   INFO = bb_cl_loadparam(FILE) returns the structure 'Info' with the 
%  values defined in FILE. Each line of the file has to have the following
%  format
%
%      .<type Id> <tag name> <value1>
%
%  The diferent types are:
%
%      Type    Id    Variable created
%      ------- ---   ----------------------------------------------
%      numeric (n) - numerical matrix, vector or scalar
%      string  (s) - string
%      text    (t) - 1xN Cell vector with strings
%      cell    (c) - Cell matrix with strings
%
%   The difference between 'text' and 'cell' is than in 'cell'
%  lines, the spaces and tabs work as delimiters between the
%  strings in differen positions of the cell
%  
%  Example :
% 
%      Let the file example.def with the next 4 lines
%
%         .n temp 120
%         .s word hello
%         .n numbers 5 4 3 2 1
%         .c list january   february march    april
%         .c list may       june     july     agost
%         .c list september october  november december
%         .t html This a HTML code <b>with bold text</b>
%         .t html with more than one line
%         .n mu 1   2  3  4  5
%         .n mu 6   7  8  9 10
%         .n mu 11 12 13 14 15
%
%      Data = bb_cl_loadparam(Spec,'example.def')
%
%      In this example, Data will be an structure with the next fields:
%
%          Data.temp      Scalar
%          Data.word      String
%          Data.numbers   1x5 Array
%          Data.list      3x4 Cell Array with strings
%          Data.html      2x1 cell array with strings
%          Data.mu        3x5 Matrix
%
%
% By Marcel Brun, e-mail:mbrun@vision.ime.usp.br, June 28 2001
% Further Information: http://www.bioinfo.usp.br/
% 
% SEE ALSO BB_CL_SAVEPARAM
%

% Changes October 9 2002, by Marcel Brun
%    a) Changed line
%             if size(Posi,1)==0
%       by
%             if isempty(Posi)
%       because in new version of Matlab it work different.
%
% Changes October 28 2002, by Marcel Brun
%    a) Now returns the name of the loaded file

% Observation Nov 05 2002, by Marcel Brun
%    a) Incompatibility between bb_cl_loadparam and bb_cl_readparam
%       when there repeated the same field as scalar, having each time
%       different length. This version pad with zeros the matrix, but
%       the previous one created a cell list with the vectors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     File Name if not defined     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin<1
   [f,p] = uigetfile('*.def');
   arch = [p f];
   if f==0
      Info = [];
      return
   end
else
   arch = File;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     Reads file     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Fid,Ferr] = fopen(arch,'rt');
if Fid==-1
   Info = [];
   disp(Ferr);
   return;
end
Text = {};
Line = 0;
while feof(Fid)==0
   Line = Line + 1;
   Text{Line} = sub_filtra(fgetl(Fid));
end
fclose(Fid);
CantLines = Line;

Info    = [];
TagList = {};
CanList = [];
Info    = {};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Reads Parameters   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:CantLines
   String = [Text{i} ' '];
   [Posi,Type,Tag,String]=sub_parse(String,TagList);
   if ~strcmp(Type,'')

      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%     new Tag     %%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%
      if Posi==0
         Posi = length(CanList)+1;
         TagList{Posi} = Tag;
         CanList(Posi)=0;
         Info(Posi).Cant=0;
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%     Update Tag info     %%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      CanList(Posi)=CanList(Posi)+1;;
      Info(Posi).Cant = CanList(Posi);
      Info(Posi).Tag  = Tag;
      Value = sub_assign(String,Type);

      if length(Value)>0

         mLine = CanList(Posi);
         if strcmp(Type,'n')
                for j=1:length(Value)
                    Info(Posi).Value(mLine,j) = Value(1,j);
                end
         elseif strcmp(Type,'x')
                for j=1:length(Value)
                    Info(Posi).Value(mLine,j) = Value(1,j);
                end
         elseif strcmp(Type,'s') 
                Info(Posi).Value = Value;
         elseif strcmp(Type,'c') 
                for j=1:length(Value)
                    Info(Posi).Value{mLine,j} = Value{1,j};
                end
         elseif strcmp(Type,'t')
                Info(Posi).Value{mLine,1} = Value;
         end
      end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     Store fields in strutcure     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Info2 = [];
for i=1:length(Info)
   Info2 = setfield(Info2,Info(i).Tag,Info(i).Value);
end
Info = Info2;

return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Sub_assign    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Value = sub_assign(String,Type)
   String2 = fliplr(deblank(fliplr(deblank(String))));
   c = [findstr(String2,' ') size(String2,2)+1 ];
   b = [ 0 c(1:size(c,2)-1) ] + 1;
   c = c-1;
   if strcmp(Type,'n') 
      Value = [];
      for k=1:size(c,2);
          Value(1,k) = str2double(deblank(String2(b(k):c(k))));
      end
   elseif strcmp(Type,'x') 
      Value = [];
      for k=1:size(c,2);
          Value(1,k) = str2double(deblank(String2(b(k):c(k))));
      end
   elseif strcmp(Type,'s')
      Value = deblank(String2(b(1):size(String2,2)));
   elseif strcmp(Type,'c')
      Value = {};
      if size(c,2)>=1
         for k=1:size(c,2);
             Value{1,k} = deblank(String2(b(k):c(k)));
         end
      end
   elseif strcmp(Type,'t')
      Value = deblank(String2(b(1):size(String2,2)));
   else
      Value = [];
   end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Determines the tag position in the tag list     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Posi,Type,Tag,String2] = sub_parse(String,TagList)

   String2 = fliplr(deblank(fliplr(deblank(String))));
   StringTag = [String2 '     '];
   if strcmp(StringTag(1:3),'.s ');
      Type = 's';
   elseif strcmp(StringTag(1:3),'.t ');
      Type = 't';
   elseif strcmp(StringTag(1:3),'.n ');
      Type = 'n';
   elseif strcmp(StringTag(1:3),'.x ');
      Type = 'x';
   elseif strcmp(StringTag(1:3),'.c ');
      Type = 'c';
   else
      Type = '';
   end
   if ~strcmp(Type,'')
      String2 = fliplr(deblank(fliplr(deblank(String2(length(Type)+2:length(String2))))));
      c = [findstr(String2,' ') size(String2,2)+1 ];
      c = c-1;
      Tag = String2(1:c(1));
      String2 = fliplr(deblank(fliplr(deblank(String2(c(1)+1:length(String2))))));
      if length(TagList)==0
         Posi = 0;
      else
         Posi = find(strcmp(Tag,TagList)==1);
         if isempty(Posi)
            Posi = 0;
         else
            Posi = Posi(1);
         end
      end
   else
      Posi   = 0;
      Tag    = '';
      String2 = '';
   end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    Deletes extra spaces or tabs     %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = sub_filtra(X)
i = 1;
j = 1;
Can = size(X,2);
Y = '';
while i<=Can
   Y(j) = X(i);
   if Y(j)==9
      Y(j) = 32;
   end
   j = j + 1;
   if X(i)==32 | X(i)==9
      while i<=Can & (X(i)==32  | X(i)==9 )
         i = i + 1;
      end
   else
      i = i + 1;
   end
end
return

