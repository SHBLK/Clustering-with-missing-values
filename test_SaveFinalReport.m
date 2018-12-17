function test_SaveFinalReport(GenericExpFolderName,ProjectName,data)
% test_SaveFinalReport saves a final report with the results of a simulation
% in a text file
%
% Author: Marco Benalcázar
%
% Email: marco_benalcazar@hotmail.com
%
% Date: May 17, 2013.
%
PathName=[GenericExpFolderName '/REPORT'];
mkdir(PathName);
FileName=fullfile(PathName,[ProjectName '.def']);
bb_cl_saveparam(FileName,data);
return