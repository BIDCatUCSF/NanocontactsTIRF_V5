%% Kyle Marchuk, PhD
% Biological Imaging Development Center at UCSF
% May 2017

% This script is the master control for studying nanocontact dynamics using
% TIRF microscopy and the Imaris output from Peter Beemiller's code.

% User input includes 3 file names, their location, and the method to
% determine +/- TCR contacts.

%% 
tic
clear all 
close all
clc
%% Load the appropriate paths and files
% Directory with the files
pathDir = 'C:\Users\Kyle\Desktop\Research\MATLAB\NanocontactsTIRF\TestData\';
% TCR intensity file
tcrName = 'OTI_ICAM_1_X_pmhc_2-2.tif';
% MapStack file
maskName = 'MapStack_OTI_ICAM_1_X_pmhc_2-2.tif';
% ContourImage file
contourName = 'contourImageStack_OTI_ICAM_1_X_pmhc_2-1.tif';
% Open the files into stacks
TCRStack = openTIFF(pathDir,tcrName);
maskStack = openTIFF(pathDir,maskName);
contourStack = openTIFF(pathDir,contourName);
% Ensure the stacks are ones and zeros
maskStack = makeBinary(maskStack,1);
contourStack = makeBinary(contourStack,1);
toc