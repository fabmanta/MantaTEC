% load the files
% close all;clear all;clc;
% addpath '/Users/fabio/Fabio_data/PhD_Singa_2014/2TEC_Project/3_code/matlab' %from cluster
% addpath '/Users/fabio/Fabio_data/PhD_Singa_2014/2TEC_Project/3_code/usefull_scripts'

% if ~exist('dateString')
%     
%     dateString = '20200124';
%     evntName = 'Mw6.7'; %Merapi or Mw7.8
%     altitutePPI = 300; % in Km
%     delta = 30; %aquisition period in seconds
%     outputFolder = sprintf('output_%s_%s_%dkm_%dsec',dateString,evntName,altitutePPI,delta);
%     outputDirectory = sprintf('/Users/fabio/Fabio_data/PostdocCNES_2019/2TECII_Project/4_output/EQ_output/%s',outputFolder);
%     
% end


files =dir(fullfile(outputDirectory,'*.mat'));
staz = [];

cd (outputDirectory)

for ii =1:length(files)  
    if size(files(ii).name,2)<= 12
    staz = [staz; files(ii).name(1:12)];
    load(files(ii).name);
    end
end
