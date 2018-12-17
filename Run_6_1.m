%Run file for the Gaussian means and inverse-Wishart covariances model.
%Simulatin setup corresponding to Fig. 3(c)
clc;
close all;
clear all;
warning off;

data.MissingRate = 0.15;         %Missing Rate
data.m=[0 0 0 0 0;0.45 0.45 0.45 0.45 0.45];                %Values for hyperparameter m   
data.v=[30 5];                   %Values for hyperparameter v
data.k=[75 75];                   %Values for hyperparameter k
data.S(:,:,1)=0.3*(75-5-1)*diag(ones(1,5));    %Scatter Matrices
data.S(:,:,2)=0.3*(75-5-1)*diag(ones(1,5));
data.NumberSamples=[35 35];     %Number of points for each class
data.ProjectName='Project_6_1';   %Name of the folder to save the results
data.NumberTestingSets=4;%40;       %Number of sets used for testing
data.ParForMode='off';          %This variable allows us to active the parfor
                                %function of Matlab to process faster the
                                %script ('on' is used to activate)
data.GenerateSamples='on';      %Allows us to generate new sets for testing
                                %(if it is on) or use previous generated 
                                %sets (if it is off)
data.MaxHammingDistance=2;%10 %in the paper
test_ExecuteModel2LargeSamples(data);
