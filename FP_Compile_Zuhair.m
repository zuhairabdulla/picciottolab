%% FP_Compile_Zuhair
%modified my code to use structs to match Ian's plotting

clear;
load(getPipelineVarsFilename);

%set the seconds before and after shock start you want.
%Note: time_after works by subtracting from shockendinds, so you'll need to
%change that to adding if you want to look past the shock, or just make it
%a negative number here.

time_before = 1;
time_after = 3;

folder = FP_PROC_DIRECTORY;
outputfolder = [FP_OUTPUT_DIRECTORY slash 'compiled'];
outputfile = 'FP Compile Zuhair MATLAB Output';

make_directory(outputfolder)

codename = 'FP_Compile_Zuhair';


%% Import Doric

%Auto Import Data
if mac
    C = dir([folder, '/*.mat']);
else
    C = dir([folder, '\*.mat']);
end

filenames = {C(:).name}.';

%exclude any temp files
filenames = filenames(~startsWith(filenames,'~'));
%remove red for now

%preallocate vars
data(length(filenames)).ID = [];


for i = 1:length(filenames)%% TODO: output PROCESSED but save (do this in subtractreferenceandsave)
    %add filename
    data(i).filenames = filenames{i};
    
    % Create the full file name and partial filename
    fullname = [folder slash C(i).name];
    load(fullname);
    % Read in the data
    %variable loaded from .mat file containing non-zscored data
    
    if size(myData, 2) == 6 %6 is the size with an extra column for red
        %columns are 1: time, 2: ref, 3: sig (green), 4: red, 5: corrected (sig-ref), DIO
        %Note: changed it to data{i,2} from {i,3} to match how i call it
        %later, other code puts it in 3.
        if COMPILE_WITH_REF %Note: This wasn't done in Eric's swag version of the code, may need to udpate it
            data(i).raw = myData(:,[1 2 6 4]);
        else %signal channel, green only
            data(i).raw  = myData(:,[1 3 6 4]);
        end
    end
    
    if COMPILE_WITH_REF %Note: This wasn't done in Eric's swag version of the code, may need to udpate it
        data(i).raw  = myData(:,[1 2 5]);
    else %signal channel, green only
        data(i).raw  = myData(:,[1 3 5]);
    end
    columnLabels = cHeader; %variable loaded from .mat file
    clear myData;
end




%% Loop through all data files

%preallocate vars for speed
all_mice_trials = cell(length(CONTROL_MICE)+length(SHOCK_MICE),2);
control_mice_trials = cell(length(CONTROL_MICE),2);
shock_mice_trials = cell(length(SHOCK_MICE),2);

all_mice_ID = strings(length(CONTROL_MICE)+length(SHOCK_MICE),1);

all_mice_mean = all_mice_trials;
all_mice_std = all_mice_trials;
all_mice_n = all_mice_trials;
all_mice_sqrt_n = all_mice_trials;
all_mice_sem = all_mice_trials;

control_mice_mean = control_mice_trials;
control_mice_std = control_mice_trials;
control_mice_n = control_mice_trials;
control_mice_sqrt_n = control_mice_trials;
control_mice_sem = control_mice_trials;

shock_mice_mean = shock_mice_trials;
shock_mice_std = shock_mice_trials;
shock_mice_n = shock_mice_trials;
shock_mice_sqrt_n = shock_mice_trials;
shock_mice_sem = shock_mice_trials;

%Note: doing cells for n just in case trials ever vary later instead of
%hardcoding 120

for file = 1:length(data)
    %Does this mouse have a red channel?
    if size(data(file).raw, 2) == 4 %4 is the size with an extra column for red
        currMouseRed = true;
    else
        currMouseRed = false;
    end
    
    %Mouse ID and day
    data(file).ID = filenames{file,1}(11:13);
    data(file).day = str2num(filenames{file,1}(17));
    
    %determine group based on group arrays and
    if any(strcmpi(data(file).ID, CONTROL_MICE)) % conrol = 0
        data(file).group = 0;
    elseif any(strcmpi(data(file).ID, SHOCK_MICE)) % shock = 1
        data(file).group = 1;
    end
    
    %exclude or not?
    if any(strcmpi(data(file).ID, EXCLUDE_FROM_COLLAPSE)) % yes = 1
        data(file).exclude = 1;
    else % no = 0
        data(file).exclude = 0;
    end
    
    
    %% prelim graphing
    if PRELIM_GRAPHS
        if VISIBLE_GRAPHS
            figure
        else
            figure('Visible', 'off')
        end
        plot(data(file).raw(:,1), data(file).raw(:,2))
        title([filenames{file,1}(11:17) ' dff0 signal'])
        print([FP_PRELIM_DIRECTORY slash filenames{file,1}(11:17) ' dff0 signal'],'-dpng')
        
        close all
    end
    
    %% Make zscore cell entry (file,3)
    dffcolumn = data(file).raw(:,2);
    
    %z-score switch
    if ZSCORE_FP
        if currMouseRed
            data(file).raw(:,4) = nanzscore(data(file).raw(:,4));
        end
        zdff = nanzscore(dffcolumn);
        
    else %no z-score
        zdff = dffcolumn;
    end
    
    data(file).raw(:,2) = zdff;
    
    if PRELIM_GRAPHS
        %prelim graph
        if VISIBLE_GRAPHS
            figure
        else
            figure('Visible', 'off')
        end
        %same col called here as prev prelim bc of zdff overrite
        plot(data(file).raw(:,1), data(file).raw(:,2))
        title([filenames{file,1}(11:17) ' zscore dff0 signal'])
        print([FP_PRELIM_DIRECTORY slash filenames{file,1}(11:17) ' zscore dff0 signal'],'-dpng')
        
        close all
    end
    
    %% Scrub weird spikes
    %ZA off for now, can put back in
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
    %designate temp variables for easy calling, make sure to clear at end
    %of loop
    
    %trim off NaNs, and overrite raw with trimmed
    NoMoreNaN = 1+max(find(isnan(data(file).raw(:,1))));
    data(file).raw = data(file).raw(NoMoreNaN:end,:);
    
    
    time = data(file).raw(:,1); %time col
    signal = data(file).raw(:,2); %signal column
    
    %set shock based on DUMMY_TTL
    if any(strcmpi(data(file).ID, DUMMY_TTL)) %if a DUMMY_TTL mouse, pull amount of TTLs needed
        shock = ones(length(time),1); %prefill with 1's in case  you're pulling fewer than you'll need
        if length(time) > length(dummy_shock_col) %if current file is bigger than dummy_shock_col, grab it all
            shock(1:length(dummy_shock_col),1) = dummy_shock_col;
        else %if dummy_shock_col is bigger, pull only as what is needed. NOTE: May be a prob if grab a partial shock
            shock(1:1:length(time),1) = dummy_shock_col(1:length(time),1);
        end
    else
        shock = data(file).raw(:,3); %ttl/dio column from own file
    end
    
    
    
    shockindices = find(~shock);
    
    shockstartinds = find(diff(shock)<0)+1;
    consecs = (find(diff(shockstartinds)==1)+1);
    shockstartinds(consecs) = [];
    %     shockendinds = find(diff(shock)>0)+1;
    %     ZA: above line is original way Ian
    %     found shockends, but it has some jitter due to the way TTLs are
    %     sometimes caught at up/downswings, which can further get messed up
    %     with downsampling. Just going to find by adding 4*SAMPLING_RATE. Left
    %     other steps below that are needed, in case you ever turn it to the
    %     old way.
    shockendinds = shockstartinds + 4*SAMPLING_RATE; %just grab 4 sec after
    consecs = (find(diff(shockendinds)==1)+1); %consecs cleans up any sloppy ttl's that are caught between 0 and 1
    shockendinds(consecs) = [];
    shockonoff = [shockstartinds shockendinds];
    
    %testing
    %         shocktime = shockstartinds - shockendinds;
    %         any(shocktime~=-120);
    %         shock_num = 18;
    %         test = data(file).raw(shockstartinds(shock_num)-1:shockendinds(shock_num)+1,:);
    %         data(file).shocktime = shocktime;
    
    %when you want to change how far AUCindices grabs, use time_after and time_before at the top.
    %This changes automatically for if you did or didn't downsample, but this will cause problems with the Native sampling rate
    %of 120.5
    AUCindices = [shockstartinds-SAMPLING_RATE*time_before shockendinds-SAMPLING_RATE*time_after];
    
    
    
    % calculate AUC for each shock window and store in output vector
    %         if contains(filename,'exp')
    AUC = zeros(length(shockonoff),1); %preallocating for speeeeeeeeeed
    for jj = 1:length(shockonoff)
        %         AUC(jj) = trapz(ZscoreStandSig(shockonoff(jj,1):shockonoff(jj,2)));%does this need to be scaled somehow to account for spacing/sampling??
        AUC(jj) = trapz(signal(AUCindices(jj,1):AUCindices(jj,2)));
        %ZA changed to signal bc don't do ZscoreStandSig
        
    end
    
    %Put temp vars into data struct in order to work with Ian's plotting
    %10/13/21: Not used for plotting anymore but used to help with a record
    %of different columns, reminding of data structure/organization, and calling based on file to
    %assing to groups, etc.
    data(file).signal = signal;
    data(file).shockindices = shockindices;
    data(file).shockstartinds = shockstartinds;
    data(file).shockendinds = shockendinds;
    data(file).shockonoff = shockonoff;
    data(file).AUCindices = AUCindices;
    data(file).AUC = AUC;
    data(file).shockonoff = shockonoff;
    
    
    
    %find signal around shock, designated by shockonoff and time_before/time_after
    for shock_idx = 1:length(shockonoff)
        sig_around_shock(:,shock_idx) = signal(shockonoff(shock_idx,1)-SAMPLING_RATE*time_before:shockonoff(shock_idx,2)-SAMPLING_RATE*time_after,1);
    end
    
    all_mice_trials{ceil(file/2),data(file).day} = sig_around_shock; %because 2 days, round up after file/2
    all_mice_ID{ceil(file/2),1} = data(file).ID; %mouse ID for all_mice
    
    %also place the data in the group's matrix, if not excluded. Excluded
    %mice will leave the row in their cell blank
    if data(file).group == 0 && data(file).exclude == 0 %0 = control
        control_mice_idx_test = strfind(CONTROL_MICE, data(file).ID); %find which mouse
        control_mice_idx = find(~cellfun(@isempty,control_mice_idx_test),1) ;
        control_mice_trials{control_mice_idx,data(file).day} = sig_around_shock; %put data into idx row
    elseif data(file).group == 1 && data(file).exclude == 0 %1 = shock. Made this elseif in case you add more groups.
        shock_mice_idx_test = strfind(SHOCK_MICE, data(file).ID); %find which mouse
        shock_mice_idx = find(~cellfun(@isempty,shock_mice_idx_test),1) ;
        shock_mice_trials{shock_mice_idx,data(file).day} = sig_around_shock; %put data into idx row
    end
    
    
    %will need to clear all this after putting the temp files somewhere,
    %putting a pin in this now to try using structs to take what Ian did
    %already for plotting
    clear shock time signal shockindices consecs shockstartinds shockendinds consecs shockonoff AUCindices AUC sig_around_shock
    
end

%% mean and sems
%using cellfun for speed and to utilize cells best.
%@(x) nanmean(x,2) is a way to designate the input args for functions used in
%cellfuns (for instance needing to take the mean over dim 2, not default of 1.
all_mice_mean = cellfun(@(x) nanmean(x, 2),all_mice_trials,'UniformOutput',false);
all_mice_std = cellfun(@(x) nanstd(x,0, 2),all_mice_trials,'UniformOutput',false);
all_mice_n = cellfun(@(x) size(x,2),all_mice_trials, 'UniformOutput',false);
all_mice_sqrt_n = cellfun(@sqrt, all_mice_n, 'UniformOutput',false);
all_mice_sem = cellfun(@rdivide, all_mice_std, all_mice_sqrt_n, 'UniformOutput',false);

control_mice_mean = cellfun(@(x) nanmean(x,2),control_mice_trials,'UniformOutput',false);
control_mice_std = cellfun(@(x) nanstd(x,0,2),control_mice_trials,'UniformOutput',false);
control_mice_n = cellfun(@(x) size(x,2),control_mice_trials, 'UniformOutput',false);
control_mice_sqrt_n = cellfun(@sqrt, control_mice_n, 'UniformOutput',false);
control_mice_sem = cellfun(@rdivide, control_mice_std, control_mice_sqrt_n, 'UniformOutput',false);

shock_mice_mean = cellfun(@(x) nanmean(x,2),shock_mice_trials,'UniformOutput',false);
shock_mice_std = cellfun(@(x) nanstd(x,0,2),shock_mice_trials,'UniformOutput',false);
shock_mice_n = cellfun(@(x) size(x,2),shock_mice_trials, 'UniformOutput',false);
shock_mice_sqrt_n = cellfun(@sqrt, shock_mice_n, 'UniformOutput',false);
shock_mice_sem = cellfun(@rdivide, shock_mice_std, shock_mice_sqrt_n, 'UniformOutput',false);

% %testing
% %with loop (slower but know it works)
% for ms = 1:size(shock_mice_trials,1)
%    for  day = 1:size(shock_mice_trials,2)
%        shock_mice_sem_loop{ms,day} = shock_mice_std{ms,day}/sqrt(size(shock_mice_trials{ms,day},2));
%    end
% end
%
% %spotcheck
% test_mean = nanmean(shock_mice_trials{1,1},2);
% test_sem = nanstd(shock_mice_trials{1,1},0,2)/sqrt(size(shock_mice_trials{1,1},2));


%% plot AUCs against shock number for each animal
%not fixed yet

% shockvec = 1:120;
% for ii = 1:length(datasplit)
%     if ~isempty(datasplit(ii,1).shockAUCs) && ~isempty(datasplit(ii,2).shockAUCs)
%         figure;
%         hold on;
%         plot(shockvec,datasplit(ii,1).shockAUCs,'k.');
%         plot(shockvec,datasplit(ii,2).shockAUCs,'r.');
%         lsline
%         xlabel('Shock number');
%         legend(append(string(datasplit(ii,1).ID),' day ',string(datasplit(ii,1).day)),...
%             append(string(datasplit(ii,2).ID),' day ',string(datasplit(ii,2).day)));
% %         if string(datasplit(ii,1).region)=='v'
% %             region = 'vHPC';
% %         elseif string(datasplit(ii,1).region)=='d'
% %             region = 'dHPC';
% %         else
% %         end
%         if usesmooth=='y'
%             titletag = append('(LOWESS; ',string(spanarg),')');
%         else
%             titletag = '';
%         end
%         title({'AUC during LH induction shocks';string(datasplit(ii,1).ID);titletag});
%
%         filetitle = append('shockAUCs LHbothdays ',string(datasplit(ii,1).ID));
%         figfilename = strrep(strrep(strrep(strrep(strrep(strrep(filetitle,'/','over'),'_',''),'(',''),')',''),' ','_'),'.','');
%         figfile = append(outputdir,'/',figfilename);
%         saveas(gcf, figfile,'fig');
%         print(figfile, '-dpng');
%     else
%     end
% end


%% For each animal, plot zscore curves for all shocks on one plot with color gradient


% Make a vector with the color gradient and designate colors
colorvec = jet(120);
set(0,'DefaultFigureColormap',feval('jet'))
green = [0.4660, 0.6740, 0.1880];
cyan = [0.3010, 0.7450, 0.9330];
gray1 = [.7 .7 .7 .25];

if REF_VS_SIG == 1 %initialize if on signal, nonz run. don't want to overrite otherwise
    SIGNAL_YLIM = NaN(length(all_mice_trials),2);
    SIGNAL_CLIM = NaN(length(all_mice_trials),2);
elseif REF_VS_SIG == 2 %load in and set desireds to signals
    load('SigVSRefLims.mat'); % if you want the same C and YLIM, uncomment below
    %SIGNAL_CLIM = SIGNAL_YLIM;
end

if REF_VS_SIG == 0 && any(DESIRED_YLIMS) == 0 && any(DESIRED_CLIMS) == 0
    ALL_YLIM = NaN(length(all_mice_trials),2);
    ALL_CLIM = NaN(length(all_mice_trials),2);
end

generic_timestamps = [0-time_before:1/SAMPLING_RATE:4-time_after]'; %e.g. -1 to 1, in 1/30 increments


%loop through each mouse in all_mice_trials
for mouse_idx = 1:length(all_mice_trials)
    
    %____original heated line plots_____
    if VISIBLE_GRAPHS
        figure
    else
        figure('Visible', 'off')
    end
    
    hold on
    
    %day 1 subplot
    p1 = subplot(1,2,1);
    hold on;
    plot(generic_timestamps,all_mice_trials{mouse_idx,1});
    p1.ColorOrder = colorvec;
    colorbar('southoutside');
    caxis([0,120]);
    title('LH induction day 1');
    yline(0);
    if ZSCORE_FP
        ylabel('Z-scored dF/F');
    else
        ylabel('dF/F');
    end
    
    yl1 = get(gca,'YLim');
    
    if any(DESIRED_YLIMS ~= 0) && REF_VS_SIG == 0 %you set lims and default run
        ylim(DESIRED_YLIMS);
    else %nothing in this else bc don't change in old way until after both
    end
    
    
    %day 2 subplot
    p2 = subplot(1,2,2);
    hold on;
    plot(generic_timestamps,all_mice_trials{mouse_idx,2}); %day 2 = col 2
    p2.ColorOrder = colorvec;
    colorbar('southoutside');
    caxis([0,120]);
    title('LH induction day 2');
    yline(0);
    if ZSCORE_FP
        ylabel('Z-scored dF/F');
    else
        ylabel('dF/F');
    end
    yl2 = get(gca,'YLim');
    
    %ylim setting
    if any(DESIRED_YLIMS ~= 0) && REF_VS_SIG == 0 %you set lims and default run
        ylim(DESIRED_YLIMS); %same for all
    elseif ~DESIRED_YLIMS && REF_VS_SIG == 0  %default run but letting set ylims for mice individually
        minylim = min(yl1(1),yl2(1));
        maxylim = max(yl1(2),yl2(2));
        ylnew = [minylim,maxylim];
        set(p1,'YLim',ylnew);
        set(p2,'YLim',ylnew);
        ALL_YLIM(mouse_idx,:) = ylnew;
    elseif ~DESIRED_YLIMS && REF_VS_SIG == 1  %default setting of ylim to get each sub
        minylim = min(yl1(1),yl2(1));
        maxylim = max(yl1(2),yl2(2));
        ylnew = [minylim,maxylim];
        set(p1,'YLim',ylnew);
        set(p2,'YLim',ylnew);
        %save ylnew to be used for ref
        SIGNAL_YLIM(mouse_idx,:) = ylnew;
    elseif any(any(SIGNAL_YLIM ~=0)) && REF_VS_SIG == 2 %read limits in after signal (not false) and ref run
        ylnew = SIGNAL_YLIM(mouse_idx,:); %read that mouse's ylims
        set(p1,'YLim',ylnew);
        set(p2,'YLim',ylnew);
    end
    
    
    %titles and save
    sgtitle({[brain_area ' ACh during LH induction shocks'];all_mice_ID(mouse_idx)});
    figfilename = [all_mice_ID{mouse_idx} ' ShockTraces LHbothdays'];
    figfile = append(outputfolder,slash,figfilename);
    saveas(gcf, figfile,'fig');
    print(figfile, '-dpng');
    
    close
    
    %_____heatmap plotting_____
    
    if VISIBLE_GRAPHS
        figure
    else
        figure('Visible', 'off')
    end
    
    
    p1 = subplot(1,2,1);
    hold on;
    
    cf = imagesc(generic_timestamps,1,all_mice_trials{mouse_idx,1}'); % not using the actual time since the number
    ylim([0 size(all_mice_trials{mouse_idx,1},2)+0.5]) %not sure why it wouldn't let me designate this whi
    colormap jet
    nanmap = [0 0 0; colormap]; %add black for NaNs
    colormap(nanmap);
    cb = colorbar;
    xlim([0-time_before 4-time_after]);
    ax = gca;
    ax.XTick = [0-time_before (0-time_before)/2 0 (4-time_after)/2 4-time_after];
    ax.TickDir = 'out';
    ax.XAxis.TickLength = [0.02 0.01];
    ax.XAxis.LineWidth = 1.75;
    ax.YAxis.TickLength = [0.02 0.01];
    ax.YAxis.LineWidth = 1.75;
    
    
    cl1 = get(gca,'CLim');
    
    if any(DESIRED_CLIMS ~= 0) && REF_VS_SIG == 0 %you set lims and default run
        ax.CLim = DESIRED_CLIMS;
    else %nothing in this else bc don't change in old way until after both
    end
    
    
    if ZSCORE_FP
        ylabel(cb, 'Z-scored dF/F');
    else
        ylabel(cb, 'dF/F');
    end
    
    %day 2 subplot
    
    p2 = subplot(1,2,2);
    hold on;
    
    cf = imagesc(generic_timestamps,1,all_mice_trials{mouse_idx,2}'); % not using the actual time since the number
    ylim([0 size(all_mice_trials{mouse_idx,2},2)+0.5]) %not sure why it wouldn't let me designate this whi
    colormap jet
    nanmap = [0 0 0; colormap]; %add black for NaNs
    colormap(nanmap);
    cb = colorbar;
    xlim([0-time_before 4-time_after]);
    ax = gca;
    ax.XTick = [0-time_before (0-time_before)/2 0 (4-time_after)/2 4-time_after];
    ax.TickDir = 'out';
    ax.XAxis.TickLength = [0.02 0.01];
    ax.XAxis.LineWidth = 1.75;
    ax.YAxis.TickLength = [0.02 0.01];
    ax.YAxis.LineWidth = 1.75;
    
    cl2 = get(gca,'CLim');
    
    %clim setting
     if any(DESIRED_CLIMS ~= 0) && REF_VS_SIG == 0 %you set lims and default run
        ax.CLim = DESIRED_CLIMS; %same for all
    elseif ~DESIRED_CLIMS && REF_VS_SIG == 0  %default run but letting set ylims for mice individually
        minclim = min(cl1(1),cl2(1));
        maxclim = max(cl1(2),cl2(2));
        clnew = [minclim,maxclim];
        set(p1,'CLim',clnew);
        set(p2,'CLim',clnew);
        ALL_CLIM(mouse_idx,:) = clnew;
    elseif ~DESIRED_CLIMS && REF_VS_SIG == 1  %default setting of ylim to get each sub
        minclim = min(cl1(1),cl2(1));
        maxclim = max(cl1(2),cl2(2));
        clnew = [minclim,maxclim];
        set(p1,'CLim',clnew);
        set(p2,'CLim',clnew);
        %save ylnew to be used for ref
        SIGNAL_CLIM(mouse_idx,:) = clnew;
    elseif any(any(SIGNAL_CLIM ~= 0)) && REF_VS_SIG == 2 %read limits in after signal (not false) and ref run
        clnew = SIGNAL_CLIM(mouse_idx,:); %read that mouse's ylims
        set(p1,'CLim',clnew);
        set(p2,'CLim',clnew);
    end
    
    if ZSCORE_FP
        ylabel(cb, 'Z-scored dF/F');
    else
        ylabel(cb, 'dF/F');
    end
    
    %titles and save
    sgtitle({[brain_area ' ACh during LH induction shocks'];all_mice_ID(mouse_idx)});
    figfilename = [all_mice_ID{mouse_idx} ' ShockHeatmaps LHbothdays'];
    figfile = append(outputfolder,slash,figfilename);
    saveas(gcf, figfile,'fig');
    print(figfile, '-dpng');
    
    close
    
    %_____mean+sem plotting_______
    
    if VISIBLE_GRAPHS
        figure
    else
        figure('Visible', 'off')
    end
    
    hold on
    
    %Transpose and grab for plotting, day 1
    mean_allSignals = all_mice_mean{mouse_idx,1}';
    sem_allSignals = all_mice_sem{mouse_idx,1}';
    %all_Signals not used here for consistency
    
    % Make a standard deviation fill for mean signal
    xx = [generic_timestamps', fliplr(generic_timestamps')];
    yy = [mean_allSignals + sem_allSignals,...
        fliplr(mean_allSignals - sem_allSignals)];
    
    %day 1 subplot
    p1 = subplot(1,2,1);
    hold on;
    %plot gray lines, comment out if you don't want
    plot(generic_timestamps,all_mice_trials{mouse_idx,1},'color',gray1);
    colorbar('southoutside', 'Visible', 'off'); %use colorbar for constant fig spacing
    
    %plot sem area
    h = fill(xx, yy, 'g'); % plot this now for overlay purposes
    set(h, 'facealpha', 0.25, 'edgecolor', 'none');
    
    % Plot the mean signals
plot(generic_timestamps, mean_allSignals, 'color', green, 'LineWidth', 1);
    
    
    title('LH induction day 1');
    yline(0);
    if ZSCORE_FP
        ylabel('Z-scored dF/F');
    else
        ylabel('dF/F');
    end
    
    yl1 = get(gca,'YLim');
    if any(DESIRED_YLIMS ~= 0) && REF_VS_SIG == 0 %you set lims and default run
        ylim(DESIRED_YLIMS);
    else %nothing in this else bc don't change in old way until after both
    end
    
    %day 2 subplot
    %Transpose and grab for plotting, day 2
    mean_allSignals = all_mice_mean{mouse_idx,2}';
    sem_allSignals = all_mice_sem{mouse_idx,2}';
    %all_Signals not used here for consistency
    
    % Make a standard deviation fill for mean signal
    xx = [generic_timestamps', fliplr(generic_timestamps')];
    yy = [mean_allSignals + sem_allSignals,...
        fliplr(mean_allSignals - sem_allSignals)];
    
    p2 = subplot(1,2,2);
    hold on;
    %plot gray lines, comment out if you don't want
    plot(generic_timestamps,all_mice_trials{mouse_idx,2},'color',gray1);
    colorbar('southoutside', 'Visible', 'off'); %use colorbar for constant fig spacing
    
    %plot sem area
    h = fill(xx, yy, 'g'); % plot this now for overlay purposes
    set(h, 'facealpha', 0.25, 'edgecolor', 'none');
    
    % Plot the mean signals
plot(generic_timestamps, mean_allSignals, 'color', green, 'LineWidth', 1);
    
    
    title('LH induction day 2');
    yline(0);
    if ZSCORE_FP
        ylabel('Z-scored dF/F');
    else
        ylabel('dF/F');
    end
    
    yl2 = get(gca,'YLim');
    
        %ylim setting
    if any(DESIRED_YLIMS ~= 0) && REF_VS_SIG == 0 %you set lims and default run
        ylim(DESIRED_YLIMS); %same for all
    elseif ~DESIRED_YLIMS && REF_VS_SIG == 0  %default run but letting set ylims for mice individually
        minylim = min(yl1(1),yl2(1)); %this probably doesn't need to be redone since it was done for first graph
        maxylim = max(yl1(2),yl2(2));
        ylnew = [minylim,maxylim];
        set(p1,'YLim',ylnew);
        set(p2,'YLim',ylnew);
    elseif ~DESIRED_YLIMS && REF_VS_SIG == 1  %default setting of ylim to get each sub
        minylim = min(yl1(1),yl2(1));
        maxylim = max(yl1(2),yl2(2));
        ylnew = [minylim,maxylim];
        set(p1,'YLim',ylnew);
        set(p2,'YLim',ylnew);
        %don't resave though
    elseif any(any(SIGNAL_YLIM ~=0)) && REF_VS_SIG == 2 %read limits in after signal (not false) and ref run
        ylnew = SIGNAL_YLIM(mouse_idx,:); %read that mouse's ylims
        set(p1,'YLim',ylnew);               %^also prob don't need to reread this
        set(p2,'YLim',ylnew);
    end
    
    %titles and save
    sgtitle({[brain_area ' ACh during LH induction shocks'];all_mice_ID(mouse_idx)});
    figfilename = [all_mice_ID{mouse_idx} ' MeanShockTraces LHbothdays'];
    figfile = append(outputfolder,slash,figfilename);
    saveas(gcf, figfile,'fig');
    print(figfile, '-dpng');

    
    close
    
end

%save lims 
if REF_VS_SIG == 1
    save('SigVSRefLims.mat', 'SIGNAL_YLIM', 'SIGNAL_CLIM')
end

if REF_VS_SIG == 0 && any(DESIRED_YLIMS) == 0 && any(DESIRED_CLIMS) == 0
    save('IndivLims.mat', 'ALL_YLIM', 'ALL_CLIM')
end

%% save excel files for Vernon
%     save_datacell = num2cell(data{file,2});
%     writecell([columnLabels; save_datacell], [outputfolder '\' mouse '.xlsx']);
%     clear save_datacell





%% Save data in file

%save all variables together
save([outputfolder slash outputfile '.mat']);



%% Print code version text file

%print the version of the code used
fileID = fopen([outputfolder slash 'codeused.txt'],'w');
fprintf(fileID, codename);




