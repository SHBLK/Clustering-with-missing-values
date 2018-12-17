%Run file for the real data example in the paper. For processing larger point sets (e.g. Fig. 3(c)) see 
%the other run file.


clc;
close all;
clear all;
warning off;


data.MissingRate = 0.15;         %Missing Rate

data.flag_reproduce=1;                %When 1, reproduces tha paper results (Fig. 4a). If 0, data_0 and data_1 should be set below.

if ~data.flag_reproduce

 load('brca_data_100_var_mean.mat');  %pre-processed feature data
 dat_all=data_all;
 data_0 = dat_all(1:(size(dat_all,1)/2),:);
 data_1 = dat_all(((size(dat_all,1)/2)+1):end,:);
 
else 
   data_0=[];  
   data_1=[];  
end


data.CalibrationSize = 90;



data.NumberSamples=[10 10];     %Number of points for each class
data.ProjectName='Project_Real_1';   %Name of the folder to save the results
data.NumberTestingSets=40;       %Number of sets used for testing
data.ParForMode='on';          %This variable allows us to active the parfor
                                %function of Matlab to process faster the
                                %script ('on' is used to activate)
data.GenerateSamples='on';      %Allows us to generate new sets for testing
                                %(if it is on) or use previous generated 
                                %sets (if it is off) 
data.SpaceSearchSuboptimal='AllCandidates';%Allows us to control the space
                                %of search of the Pmax suboptimal algorithm
                                %'AllCandidates' tells the Pmax suboptimal
                                %algorithm to use a set of candidate
                                %partitions composed of all the partitions
                                %for a given set of points. 'Possible' uses
                                %a set of candidate partitions composed of
                                %only those partitions that assign the
                                %correct number of points for each cluster
test_ExecuteModel2Real(data,data_0,data_1)