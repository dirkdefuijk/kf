function [FileName, PathName] = make_list_of_data_paths()
%% make_list_of_data.m
% !!!! NOTE: SHOULD BE EDITED MANUALLY WHEN NEW DATA IS ADDED !!!!!!

% List of all available stall data sets, where a selection can be made of
% which to include in the parameter estimation routine.
% Joost van Ingen 27-10-2016

% Legend:
% MAN  = stall manoeuvre type (QS = quasi-static, DYN = dynamic)
% CONF = configuration (examples: clean, f15g = flaps 15 deg, gear down)
% M[-] = average Mach number during manoeuvre
% h[ft]= average altitude during manoeuvre

%% Create the actual list
% path where flight path reconstructed data is stored
data_folder = '/Users/Joost/Dropbox/thesis/07_code/matlab/data/processed';

% create cell array with paths to data sets
FullFilePathList = { ...
%
% ... %%% FL0-FL50 %%%                                              %MAN    %CONF   %M[-]   %h[ft]  %notes
% [data_folder,'/FL0-FL50','/M0.1-M0.4','/stall_FL0-FL50_M0.1-M0.4_set001.mat']; ...  QS      clean   0.19    4903    
% [data_folder,'/FL0-FL50','/M0.1-M0.4','/stall_FL0-FL50_M0.1-M0.4_set002.mat']; ...  QS      clean   0.19    4927
% [data_folder,'/FL0-FL50','/M0.1-M0.4','/stall_FL0-FL50_M0.1-M0.4_set003.mat']; ...  QS      clean   0.20    4983
% [data_folder,'/FL0-FL50','/M0.1-M0.4','/stall_FL0-FL50_M0.1-M0.4_set004.mat']; ...  QS      clean   0.18    4853
%
% ... %%% FL50-FL80 %%%                                             %MAN    %CONF   %M[-]   %h[ft]  %notes
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set001.mat']; ... QS      clean   0.20    6838
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set002.mat']; ... QS      clean   0.20    6853
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set003.mat']; ... QS      clean   0.20    6903
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set004.mat']; ... QS      clean   0.19    7011
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set005.mat']; ... QS      clean   0.21    6914
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set006.mat']; ... QS      clean   0.20    6906
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set007.mat']; ... QS      clean   0.21    6931
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set008.mat']; ... QS      clean   0.20    6914
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set009.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set010.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set011.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set012.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set013.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set014.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set015.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set016.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set017.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set018.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set019.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set020.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set021.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set022.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set023.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set024.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set025.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set026.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set027.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set028.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set029.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set030.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set031.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set032.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set033.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set034.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set035.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set036.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set037.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set038.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set039.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set040.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set041.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set042.mat']; ... QS      clean
% [data_folder,'/FL50-FL80','/M0.1-M0.4','/stall_FL50-FL80_M0.1-M0.4_set043.mat']; ... QS      clean
%
% ... %%% FL80-FL110 %%%                                                               %MAN     %CONF   %M[-]   %h[ft]  %notes
% [data_folder,'/FL80-FL110','/M0.1-M0.4','/stall_FL80-FL110_M0.1-M0.4_set001.mat']; ...QS      clean 
% [data_folder,'/FL80-FL110','/M0.1-M0.4','/stall_FL80-FL110_M0.1-M0.4_set002.mat']; ...QS      clean 
% [data_folder,'/FL80-FL110','/M0.1-M0.4','/stall_FL80-FL110_M0.1-M0.4_set003.mat']; ...QS      clean 
% [data_folder,'/FL80-FL110','/M0.1-M0.4','/stall_FL80-FL110_M0.1-M0.4_set004.mat']; ...QS      clean 
% [data_folder,'/FL80-FL110','/M0.1-M0.4','/stall_FL80-FL110_M0.1-M0.4_set005.mat']; ...QS      clean 
% [data_folder,'/FL80-FL110','/M0.1-M0.4','/stall_FL80-FL110_M0.1-M0.4_set006.mat']; ...QS      clean 
% [data_folder,'/FL80-FL110','/M0.1-M0.4','/stall_FL80-FL110_M0.1-M0.4_set007.mat']; ...QS      clean 
% [data_folder,'/FL80-FL110','/M0.1-M0.4','/stall_FL80-FL110_M0.1-M0.4_set008.mat']; ...QS      clean 
% [data_folder,'/FL80-FL110','/M0.1-M0.4','/stall_FL80-FL110_M0.1-M0.4_set009.mat']; ...QS      clean 
% [data_folder,'/FL80-FL110','/M0.1-M0.4','/stall_FL80-FL110_M0.1-M0.4_set010.mat']; ...QS      clean 
% [data_folder,'/FL80-FL110','/M0.1-M0.4','/stall_FL80-FL110_M0.1-M0.4_set011.mat']; ...QS      clean 
[data_folder,'/FL80-FL110','/M0.1-M0.4','/stall_FL80-FL110_M0.1-M0.4_set012.mat']; ...DYN1.0g clean 
[data_folder,'/FL80-FL110','/M0.1-M0.4','/stall_FL80-FL110_M0.1-M0.4_set013.mat']; ...DYN1.0g clean 
%
% ... %%% FL110-FL150 %%%                                           %MAN     %CONF   %M[-]   %h[ft]  %notes
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set001.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set002.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set003.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set004.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set001.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set002.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set003.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set004.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set005.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set006.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set007.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set008.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set009.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set010.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set011.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set012.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set013.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set014.mat'];...QS      clean 
% [data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set015.mat'];...QS      clean 
[data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set016.mat'];...DYN1.0g clean 
[data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set017.mat'];...DYN1.0g clean 
[data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set018.mat'];...DYN1.0g clean 
[data_folder,'/FL110-FL150','/M0.1-M0.4','/stall_FL110-FL150_M0.1-M0.4_set019.mat'];...DYN1.0g clean 
%
% ... %%% FL150-FL200 %%%                                           %MAN     %CONF   %M[-]   %h[ft]  %notes
% [data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set001.mat'];...QS      clean 
% [data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set002.mat'];...QS      clean 
% [data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set003.mat'];...QS      clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set004.mat'];...DYN1.0g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set005.mat'];...DYN1.0g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set006.mat'];...DYN1.0g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set007.mat'];...DYN1.0g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set008.mat'];...DYN1.0g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set009.mat'];...DYN1.0g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set010.mat'];...DYN1.0g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set011.mat'];...DYN1.0g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set012.mat'];...DYN1.0g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set013.mat'];...DYN1.1g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set014.mat'];...DYN1.1g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set015.mat'];...DYN1.1g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set016.mat'];...DYN1.1g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set017.mat'];...DYN1.3g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set018.mat'];...DYN1.3g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set019.mat'];...DYN1.3g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set020.mat'];...DYN1.1g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set021.mat'];...DYN1.1g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set022.mat'];...DYN1.1g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set023.mat'];...DYN1.1g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set024.mat'];...DYN1.1g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set025.mat'];...DYN1.1g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set026.mat'];...DYN1.3g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set027.mat'];...DYN1.3g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set028.mat'];...DYN1.3g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set029.mat'];...DYN1.3g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set030.mat'];...DYN1.0g  clean 
[data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL150-FL200_M0.1-M0.4_set031.mat'];...DYN1.0g  clean 
%
% ... %%% FL200-FL400 %%%                                           %MAN     %CONF   %M[-]   %h[ft]  %notes
% [data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL200-FL400_M0.1-M0.4_set001.mat'];...QS       clean 
% [data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL200-FL400_M0.1-M0.4_set002.mat'];...QS       clean 
% [data_folder,'/FL150-FL200','/M0.1-M0.4','/stall_FL200-FL400_M0.1-M0.4_set003.mat'];...QS       clean 
};

% create empty cells for paths and filenames
FileName = cell(length(FullFilePathList),1);
PathName = cell(length(FullFilePathList),1);

% split file paths and file names
    for i = 1:length(FullFilePathList)
        [PathName{i},FileName{i}] = fileparts(FullFilePathList{i});
        FileName{i} = [FileName{i},'.mat'];
    end
        
end
