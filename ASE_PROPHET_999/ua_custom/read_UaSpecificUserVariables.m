function UserVar=read_UaSpecificUserVariables(UserVar)

%InputDataDirectory='/work/n02/n02/janryd69/UaMITgcm/Ua_InputData';
InputDataDirectory='/media/wchm8/mainJDRydt2/UaMITgcm_v2/Ua_InputData';

RunTable=readtable([InputDataDirectory,'/RunTable.csv']);

I = find(strcmp(RunTable{:,'ID'},UserVar.UaMITgcm.Experiment)) ;
if isempty(I)
    error('Experiments unknown');
end

switch RunTable{I(1),'GeometryInterpolants'}{:}
    case 'Bedmachine_Bamber2009'
        UserVar.Geometry = 'Bedmachine_Bamber2009';
        UserVar.GeometryInterpolants = [InputDataDirectory,'/GriddedInterpolants_sBh_Bedmachine_Bamber2009.mat'];
        UserVar.FirnInterpolants = [InputDataDirectory,'/GriddedInterpolants_Firn_RACMO.mat'];
    case 'Bedmachine_Cryosat2'
        UserVar.Geometry = 'Bedmachine_Cryosat2';
        UserVar.GeometryInterpolants = [InputDataDirectory,'/GriddedInterpolants_sBh_Bedmachine_Cryosat2.mat'];
        UserVar.FirnInterpolants = [InputDataDirectory,'/GriddedInterpolants_Firn_RACMO.mat'];
    case 'Bedmap2'
        UserVar.Geometry = 'Bedmap2';
        UserVar.GeometryInterpolants = [InputDataDirectory,'/GriddedInterpolants_sBh_Bedmap2.mat'];
        UserVar.FirnInterpolants = [InputDataDirectory,'/GriddedInterpolants_Firn_RACMO.mat'];
    case 'Cryosat2'
        UserVar.Geometry = 'Cryosat2';
        UserVar.GeometryInterpolants = [InputDataDirectory,'/GriddedInterpolants_sBh_Millan_Cryosat2.mat'];
        UserVar.FirnInterpolants = [InputDataDirectory,'/GriddedInterpolants_Firn_RACMO.mat'];
    case 'REMA'
        UserVar.Geometry = 'REMA';
        UserVar.GeometryInterpolants = [InputDataDirectory,'/GriddedInterpolants_sBh_Millan_REMA.mat'];
        UserVar.FirnInterpolants = [InputDataDirectory,'/GriddedInterpolants_Firn_RACMO.mat'];
    case 'Bedmachine'
        UserVar.Geometry = 'Bedmachine';
        UserVar.GeometryInterpolants = [InputDataDirectory,'/GriddedInterpolants_sBh_Bedmachine.mat'];
        UserVar.FirnInterpolants = 'Bedmachine';
    case 'Bedmachine_Bamber2009_modifiedThwaites'
        UserVar.Geometry = 'Bedmachine_Bamber2009_modifiedThwaites';
        UserVar.GeometryInterpolants = [InputDataDirectory,'/GriddedInterpolants_sBh_Bedmachine_Bamber2009_modifiedThwaites'];
        UserVar.FirnInterpolants = [InputDataDirectory,'/GriddedInterpolants_Firn_RACMO.mat'];
    case 'Bedmachine20190905_Bamber2009'
	UserVar.Geometry = 'Bedmachine_Bamber2009';
        UserVar.GeometryInterpolants = [InputDataDirectory,'/GriddedInterpolants_sBh_Bedmachine2019-09-05_Bamber2009.mat'];
        UserVar.FirnInterpolants = [InputDataDirectory,'/GriddedInterpolants_Firn_RACMO.mat'];
    case 'Bedmachine20200715_Bamber2009'
        UserVar.Geometry = 'Bedmachine_Bamber2009';
        UserVar.GeometryInterpolants = [InputDataDirectory,'/GriddedInterpolants_sBh_Bedmachine2020-07-15_Bamber2009.mat'];
        UserVar.FirnInterpolants = [InputDataDirectory,'/GriddedInterpolants_Firn_RACMO.mat'];
    case 'Bedmachinev3_2000'
        UserVar.Geometry = 'Bedmachinev3_2000';
        UserVar.GeometryInterpolants = [InputDataDirectory,'/GriddedInterpolants_Geometry_01-Jun-2000_EXTRUDED.mat'];
        UserVar.FirnInterpolants = UserVar.GeometryInterpolants;        
    otherwise
        error('Geometry unknown');    
end

Experiment_fullname = RunTable{I(1),'Row'}{:}; 

UserVar.NameOfFileForReadingSlipperinessEstimate=[InputDataDirectory,'/',Experiment_fullname,'_C-Estimate.mat'];
UserVar.NameOfFileForReadingAGlenEstimate=[InputDataDirectory,'/',Experiment_fullname,'_AGlen-Estimate.mat'];

UserVar.RACMO_SMB=[InputDataDirectory,'/ScatteredInterpolants_SMB.mat'];

UserVar.SlidingLaw = RunTable{I(1),'SlidingLaw'}{:};
