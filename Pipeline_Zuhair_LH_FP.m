%% Pipeline_Zuhair_LH_FP
% IMPORTANT NOTES:

%9/29/21: Modified for Zuhair, lots of legacy stuff in all parts of this
%pipline and scripts it calls.

%10/6/21: edited to work more easily between my and Zuhair's computer.
%Changing lots of stuff to do cooler things.

%1) Before running this code please make sure that all the
% pipeline scripts are in the same directory, otherwise this isn't
% guaranteed to run! This should be the case if you downloaded everything
% and didn't rearrange anything.

%2) Make sure MATLAB's current folder is this pipeline's parent folder.
%Don't add all the folders to the path because there are different versions of the same
%script.

%3) Change the FP_PARENT_DIRECTORY and FP_RAW_DIRECTORY below


% PLEASE NOTE THAT MOUSE 850 HAD A BROKEN REFERENCE CHANNEL SO IF YOU ARE
% COMPILING THE DATA USING REFERENCE INSTEAD OF SIGNAL THIS MOUSE WILL BE
% EXCLUDED

%clear MATLAB's memory to prevent issues
clear;



%% Things to change as needed

%Brain area
brain_area = 'BLA';

%Treatment Groups
CONTROL_MICE = ["129";"130";"133";"134";"139";"140";"143";"144"];
SHOCK_MICE = ["131";"132";"135";"136";"137";"141";"142";"145";"146";"147"];

%Don't use these mice in the averaged plots
EXCLUDE_FROM_COLLAPSE = ["137";"144"];

%Either designate here or we can call a script in the pipeline to auto
%determine and designate.
RESILIENT_MICE = [];
HELPLESS_MICE = [];

%zuhair or rick
%if mac = true, use mac slashses, / instead of \ and use zuhair's directories
mac = true;

DOWNSAMPLE_FP_DATA = true; %set this to true to downsample your data from native sampling rate to downsample rate

DESIRED_YLIMS = [-4 10]; %set the lower and upper ylims to be used in all graphs, e.g.: [-4 10]
%set to false to let code match within subplots
%or when you want to compare

DESIRED_CLIMS = [-4 10]; %set the lower and upper clims to be used in heatmaps, e.g.: [-4 10]
%set to false to let code match within subplots
%or when you want to compare

%NOTE: if REF_VS_SIG > 0, lims designated here are overwritten (failsafe)

REF_VS_SIG = 0; 
% 0 = Zscore true, compile_with_ref false (default)
% 1 = Zscore false, compile_with_ref false
% 2 = Zscore false, compile_with_ref true
% Note: If going to ref_vs_sig, must run 1 first, to save
% the ylims for running 2.
% Doing so will overwrite Desired Y and Clims designated
% here to allow for it to run as designed (failsafe)

%Scenarios for related switches: DESIRED_YLIMS, DESIRED_CLIMS, REF_VS_SIG

%A) Most common: You know limits you want to use for all mice
% DESIRED_YLIMS = [-4 10]; % DESIRED_CLIMS = [-4 10]; % REF_VS_SIG = 0;

%B) You want to let it decide the lims individually for each mouse.
% DESIRED_YLIMS = false; % DESIRED_CLIMS = false; % REF_VS_SIG = 0;
% check the saved var IndivLims.mat to see which lims you wanna use

%C) You want to compare sig and ref, non-zscored, but using the same lims.
%First, run: DESIRED_YLIMS = false; % DESIRED_CLIMS = false; % REF_VS_SIG = 1;
%Second, run: DESIRED_YLIMS = false; % DESIRED_CLIMS = false; % REF_VS_SIG = 2;



VISIBLE_GRAPHS = false; %Put as false if you don't want the graphs popping up while code is running

DUMMY_TTL = CONTROL_MICE; %designate which mice need dummy TTLS this to true if control mice don't
%e.g. DUMMY_TTL = ["420"; "421"; "321"];
%put DUMMY_TTL = CONTROL_MICE if all need it
%[], not false, if not needed (i.e. all mice have their own ttl)
%NOTE: currently using a manually saved
%dummy_shock_col.mat

PRELIM_GRAPHS = true; %turn off if you don't want the prelim graphs (for speed)



%% switch flippin

% Make sure to change this directory to the parent folder of this pipeline.
% ALSO MAKE SURE TO CHANGE RAW DIRECTORY TO YOUR RAW FP DATA
% e.g.: 'C:\Users\rbc52\Documents\MATLAB\Crouse et al\2019-06'
%'C:\Users\egira\Documents\GitHub\CueReward\2019-14 ACh3.0';
%'C:\Users\egira\Desktop\raw'

if mac %mac = true, usehair's directory
    slash = '/';
    FP_PARENT_DIRECTORY = '/Users/zia2/OneDrive - Yale University/Desktop/Picciotto Lab/C21 LH + FP/FP/BLA/10-14 code/Ind/';
    FP_RAW_DIRECTORY = '/Users/zia2/OneDrive - Yale University/Desktop/Picciotto Lab/C21 LH + FP/FP/BLA/10-14 code/Ind/raw';
    FP_RAW_PRE_SPLIT_DIRECTORY = '/Users/zia2/OneDrive - Yale University/Desktop/Picciotto Lab/C21 LH + FP/FP/BLA/10-14 code/Ind/raw_pre_split';
    %if FP_PARENT_DIRECTORY wasn't designated, stop script and alert the user
    if isempty(FP_PARENT_DIRECTORY)
        fprintf('Please designate FP_PARENT_DIRECTORY in this pipeline and make sure the current folder for MATLAB is pointed in it');
        return
    end
else %mac = false, use pc slashes and rick's directory
    slash = '\';
    FP_PARENT_DIRECTORY = 'C:\Users\rbc52\Desktop\Zuhair\Ind';
    FP_RAW_DIRECTORY = 'C:\Users\rbc52\Desktop\Zuhair\Ind\raw';
    FP_RAW_PRE_SPLIT_DIRECTORY = 'C:\Users\rbc52\Desktop\Zuhair\Ind\raw_pre_split';
    %if FP_PARENT_DIRECTORY wasn't designated, stop script and alert the user
    if isempty(FP_PARENT_DIRECTORY)
        fprintf('Please designate FP_PARENT_DIRECTORY in this pipeline and make sure the current folder for MATLAB is pointed in it');
        return
    end
end

%Downsampling
NATIVE_SAMPLING_RATE = 120.5; %set to native sampling rate

if DOWNSAMPLE_FP_DATA
    SAMPLING_RATE = 30; %Set this to your desired sampling rate of your FP data, in data points per second.
else %If you are not downsampling, it will be native sampling rate
    SAMPLING_RATE = NATIVE_SAMPLING_RATE;
end


%REF_VS_SIG switching
if REF_VS_SIG == 0
    ZSCORE_FP = true;
    COMPILE_WITH_REF = false;
elseif REF_VS_SIG == 1
    ZSCORE_FP = false;
    COMPILE_WITH_REF = false;
    DESIRED_YLIMS = false; %overrite y and clims to let the code determine for each mouse and save to be loaded in 2
    DESIRED_CLIMS = false;
elseif REF_VS_SIG == 2
    ZSCORE_FP = false;
    COMPILE_WITH_REF = true;
    DESIRED_YLIMS = false; %overrite y and clims to let the code determine for each mouse and save to be loaded in 2
    DESIRED_CLIMS = false;
end

%compile with reference or signal channel output folder
if COMPILE_WITH_REF
    FP_OUTPUT_DIRECTORY = [FP_PARENT_DIRECTORY slash 'generated output (reference)'];
else
    FP_OUTPUT_DIRECTORY = [ FP_PARENT_DIRECTORY slash 'generated output' ];
end

%individual limits
if REF_VS_SIG == 0 && any(DESIRED_YLIMS) == 0 && any(DESIRED_CLIMS) == 0
    FP_OUTPUT_DIRECTORY = [FP_PARENT_DIRECTORY slash 'generated output indiv lims'];
else
    FP_OUTPUT_DIRECTORY = [ FP_PARENT_DIRECTORY slash 'generated output' ];
end

%zscoring
if ZSCORE_FP
    %keep output directory the same
else
    FP_OUTPUT_DIRECTORY = [ FP_OUTPUT_DIRECTORY ' nonzscored' ];
end




%DUMMY_TTL
if ~isempty(DUMMY_TTL)
    load('dummy_shock_col.mat')
end

%% directory making

FP_PROC_DIRECTORY = [ FP_PARENT_DIRECTORY slash 'processed' ];

FP_PRELIM_DIRECTORY = [ FP_OUTPUT_DIRECTORY slash 'generated prelim' ];

FP_TIMESTAMP_FILE = [ FP_OUTPUT_DIRECTORY slash 'pipeline timestamps.xlsx' ];


save(getPipelineVarsFilename);

make_directory(FP_OUTPUT_DIRECTORY)


%%
split_files(FP_RAW_PRE_SPLIT_DIRECTORY,FP_RAW_DIRECTORY)
% 
Basic_FP_processing(FP_RAW_DIRECTORY,FP_PROC_DIRECTORY)

FP_Compile_Zuhair
% FP_Compile_Zuhair_120hz
