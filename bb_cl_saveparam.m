function bb_cl_saveparam(mFile, Data)
% BB_CL_SAVEPARAM   Writes information to specification file
%
%    bb_cl_saveparam(File,Data) 
%
%    Input:
%        File    - Text file where will be saved the information
%        Data    - Structure with the information
%
% 
%   bb_cl_saveparam(FILE,DATA) saves the data in the structure DATA
%   to the file FILE. The field names will be used as tag names.
%   The type will be determined by the data type for each field,
%   with the next rules:
%
%      numeric (.n) - numerical matrices (can be scalar or vector)
%      numeric (.x) - Cell array with identical size matrices
%      string  (.s) - Strings
%      text    (.t) - 1xN Cell vector with strings containing spaces
%      cell    (.c) - NXM Cell array with strings (if it is not "text")
%
% If the field is a cell array with matrices, they will be saved
% as a unique matrix, concatenated by the 1st dimention. When loaded
% with bb_cl_loadparam, it will be recognized simply as a matrix.
% Further update in bb_cl_loadparam will enables to recognize this
% kind of field.
%
% By Marcel Brun (mbrun@vision.tamu.edu), June 28 2001
% 
% SEE ALSO BB_CL_LOADPARAM

% Changes August 28 2001, By Marcel Brun
%    a) Corrected bug when saving cell arrays with matrices
%       (was using the variable "i" for the loop


Names = fieldnames(Data);
if strcmp(computer,'PCWIN')
   mFin  = [char(13) char(10)];
else
   mFin  = [char(10)];
end
mText = [];

for i=1:size(Names,1)
   Name  = Names{i};
   Values = getfield(Data,Name);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%     Cell Array      %%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if iscell(Values)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%%%%     Detects type of cell array   %%%%%
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Type = 'c';
      if length(Values)>1 & isnumeric(Values{1})
         Type = 'x';
      else
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%     Detects if it s text list    %%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         IsText = 0;
         if size(Values,2)==1
            for j=1:length(Values)
               if length(findstr(' ',Values{j}))>0 | length(findstr(char(9),Values{j}))>0
                  IsText = 1;
                  j=length(Values);
               end
            end
         end
         if IsText==1
            Type = 't';
         end
      end
      
      if strcmp(Type,'t')
         for j=1:length(Values);
             mText = [mText '.t ' Name ' '  num2str(Values{j}) mFin];
         end
         mText = [mText ' '  mFin];
      elseif strcmp(Type,'x')
         for j=1:length(Values)
             Matr = Values{j};
             for j = 1:size(Matr,1)
                 mText = [mText '.x ' Name ' '  num2str(Matr(j,:)) mFin];
                 pause(0)
             end
             mText = [mText mFin];
         end
      else
         for j = 1:size(Values,1)
             mText = [mText '.c ' Name ' ' ];
             for k = 1:size(Values,2)
                 mText = [mText Values{j,k} ' '];
                 pause(0)
             end
             mText = [mText mFin];
         end
         mText = [mText mFin];
      end
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%     Matrix     %%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%
   elseif isnumeric(Values)
      for j=1:size(Values,1)
          mText = [mText '.n ' Name ' '  num2str(Values(j,:)) mFin];
          pause(0)
      end
      mText = [mText ' '  mFin];
      
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%     Non empty String     %%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   elseif ischar(Values) & length(Values)>0
      mText = [mText '.s ' Name ' '  Values mFin];
      mText = [mText ' '  mFin];
   end
   pause(0)
end

fid = fopen(mFile,'wb');
fwrite(fid,mText,'uchar');
fclose(fid);
