%% FP_Shock_Test


clear;
load(getPipelineVarsFilename);

folder = [FP_PARENT_DIRECTORY '\generated output\generated processed'];
outputfolder = [FP_PROC_DIRECTORY '\Shock testing output'];
outputfile = 'FP Shock Test MATLAB Output';

make_directory(outputfolder)

codename = 'FP_Shock_Test';


%% Import Doric

%Auto Import Data

C = dir([folder, '\*.mat']);
filenames = {C(:).name}.';

%exclude any temp files
filenames = filenames(~startsWith(filenames,'~'));
%remove red for now

%
data = filenames; %data columns: filenames, nonzscored, zscored, [medsum medata], ???
raw = cell(length(filenames),1);

for i = 1:length(filenames)%% TODO: output PROCESSED but save (do this in subtractreferenceandsave)
    % Create the full file name and partial filename
    fullname = [folder '\' C(i).name];
    load(fullname);
    % Read in the data
    %variable loaded from .mat file containing non-zscored data
    if size(myData, 2) == 6 %6 is the size with an extra column for red
        %columns are time, ref, sig (green), red, corrected (sig-ref), DIO
        %Note: changed it to data{i,2} from {i,3} to match how i call it
        %later, other code puts it in 3.
        data{i,2} = myData;
    else
        data{i,2} = myData(:,[1 3 5]);
    end
    columnLabels = cHeader; %variable loaded from .mat file
    clear myData;
end




%Cutting headers off when pulling now so I can make then matrices
% %cut off column headings
% for row = 1:length(data)
%     data{row,2} = data{row,2}(2:end,:);
% end



%% Loop through all data files
for file = 1:length(data)
    
    
    
    
    %% Scrub weird spikes
    %turn weird spikes to NaN

    %TODO check to see if any of the cols have something greater than |100|
    %before scrubbing
    
    %what is datarow even doing here? looping for no reason? 
%     for datarow=1:size(data{file,2},1)
%         for scrubrow=[2:5]
%             data{file,2}(data{file,2}(:,scrubrow)<-100,2) = NaN;
%             data{file,2}(data{file,2}(:,scrubrow)>100,2) = NaN;
%         end
%     end
    
    %% find pulses
    %find when pulses start and end
    
    %basic_FP_processing already trims to the first pulse
    %start 3 sec after beginning to find the next pulses
    
    %change to use sampling rate instead, but after downsample bc of 120.5
    %sampling rate
    %for idx = 3*SAMPLING_RATE:size(data{file,2}(:,3))
    
    pulse_start_idx = 0;
    pulse_end_idx = 0; 
    for idx = 4*SAMPLING_RATE:size(data{file,2}(:,6),1)
        %find start of pulses
        if data{file,2}(idx,6) < 1 && data{file,2}(idx-1,6) == 1
            pulse_start_idx = pulse_start_idx + 1;
            pulse_start(pulse_start_idx,1) = idx;
        end
        
        %find end of pulses
        if data{file,2}(idx,6) == 1 && data{file,2}(idx-1,6) < 1
            pulse_end_idx = pulse_end_idx + 1;
            pulse_end(pulse_end_idx,1) = idx;
        end
    end
    
    %     num_of_pulses = round(size(pulses,1)/61);
    num_of_pulses = 2;
    
    mouse = data{file,1}(11:14);
    
    
    
    %% plot sig only
    
    figure
    hold on
    
    plot(data{file,2}(:,1), data{file,2}(:,3), 'Color', [0.2549 0.8314  0.1647 0.4])
    
    %this must have been based on the old/manual shock testing
    %     for pulse_idx = 0:num_of_pulses-1
    %         xline(data{file,2}(pulses(pulse_idx*61+1),1), '--c');
    %
    %     end
    %
    
    %only draw line at shock for now (10 sec after pulse start, which is
    %when
    xline(data{file,2}(pulse_end(1),1), '--c', 'LineWidth', 1.25, 'Alpha', 0.3);
    
    
    print([outputfolder '\' mouse], '-dpng');
    
    
    xlim([data{file,2}(pulse_start(1)-(SAMPLING_RATE*10),1) data{file,2}(pulse_start(1)+(SAMPLING_RATE*30),1)]);
    
    print([outputfolder '\' mouse ' shock zoom' ], '-dpng');
    
    
    %% red only
    figure
    hold on
    
    plot(data{file,2}(:,1), data{file,2}(:,4), 'Color',[1     0     1  0.4])
    
    %only draw line at shock for now
    xline(data{file,2}(pulse_end(1),1), '--c', 'LineWidth', 1.25, 'Alpha', 0.3);
    
    
    print([outputfolder '\' mouse ' red only'], '-dpng');
    
    %had to change this up to match the indexing method here
    xlim([data{file,2}(pulse_start(1)-(SAMPLING_RATE*10),1) data{file,2}(pulse_start(1)+(SAMPLING_RATE*30),1)]);
    
    print([outputfolder '\' mouse ' red only shock zoom' ], '-dpng');
    
    %% sig + ref + red
    %ref = (blue) [ 0     0     1 0.3]
    %sig = green) [0.2549    0.8314    0.1647 0.3]
    %red channel = (magenta) [1     0     1  0.4]
    
    figure
    hold on
    %ref
    plot(data{file,2}(:,1), data{file,2}(:,2), 'Color',[ 0     0     1 0.3]);
    
    %sig
    plot(data{file,2}(:,1), data{file,2}(:,3), 'Color',[0.2549    0.8314    0.1647 0.3]);
    
    %red ch
    plot(data{file,2}(:,1), data{file,2}(:,4), 'Color', [1     0     1  0.3]);
    
    %only draw line at shock for now
    xline(data{file,2}(pulse_end(1),1), '--c', 'LineWidth', 1.25, 'Alpha', 0.3);
    
    print([outputfolder '\' mouse ' sig+ref+red'], '-dpng');
    
    xlim([data{file,2}(pulse_start(1)-(SAMPLING_RATE*10),1) data{file,2}(pulse_start(1)+(SAMPLING_RATE*30),1)]);
    
    print([outputfolder '\' mouse ' sig+ref+red shock zoom' ], '-dpng');
    
    close all
    
    %save excel files for Vernon 
    save_datacell = num2cell(data{file,2});
    writecell([columnLabels; save_datacell], [outputfolder '\' mouse '.xlsx']);
    clear save_datacell
end




%% Save data in file

%save all variables together
save([outputfolder '\' outputfile '.mat']);



%% Print code version text file

%print the version of the code used
fileID = fopen([outputfolder '\codeused.txt'],'w');
fprintf(fileID, codename);




