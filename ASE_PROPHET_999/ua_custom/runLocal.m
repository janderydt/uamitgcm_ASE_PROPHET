function runLocal

addpath(genpath('/media/janryd69/mainJDeRydt/UaMITgcm_v2/UaMITgcm_git/UaSource_beta_19Jan2024'));
addpath(genpath('/media/janryd69/mainJDeRydt/UaMITgcm_v2/UaMITgcm_git/coupling'));

fileID = fopen('options_for_ua','w');
fprintf(fileID,'ASE_PROPHET_999\n');
fprintf(fileID,'/mnt/SSD1/Documents/Projects/PROPHET/uamitgcm_ASE_PROPHET/ASE_PROPHET_999/output/\n');
fprintf(fileID,'calendar\nNewMeltrate.mat\nDataForMIT.mat\nmatlab\nxy');
fclose(fileID);

restartfile = string(input("Path to restart file to read: "));
copyfile(restartfile,'./ASE_PROPHET_999-RestartFile.mat');

callUa;