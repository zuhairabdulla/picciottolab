function Basic_FP_processing = Basic_FP_processing(directory,output_directory,slash)

load(getPipelineVarsFilename); %load pipeline variables

%ZA not trimming here due to shock happening as soon as program starts,
%can change later
autoTrim = false;
printExcel = false; %to print the processed as excel
%list the numbers of mice with both red and green channels
redAndGreenMice = {};
if not(autoTrim)
    trimTime = 0; %Time to trim off in seconds, ZA put to 0 for no trim
end
start = 0;

% Change this directory to the folder containing your raw doric files!
savetype = 'matlab';

files = dir(directory);
files = files(contains({files.name},{'.csv'}));

processed_files = dir(output_directory);

%remove non csv files
processed_files = [ processed_files(:).name ];

make_directory(output_directory); %RIC i fixed makedirectory to be cooler
make_directory(FP_PRELIM_DIRECTORY);

for file = files'
    
    filename = strcat(file.name);
    %only process .csv files, don't process "PROCESSED" files, and don't
    %process any that already have a 'PROCESSED' version in the folder
    if isempty(strfind(filename, '.csv')) || isempty(strfind(filename, 'PROCESSED_'))==false || sum(strcmp(strcat('PROCESSED_',filename),{files.name}))>0 || ~isempty(strfind(processed_files, filename)) || ( COMPILE_WITH_REF  && contains(filename, '0850') )
        fprintf('Skipping %s\n', filename);
        continue
    end
    
    allData = readmatrix([directory, slash filename]); % 1: skip first two lines line (header); might need to skip more depeding how the file but basically the goal is to scrap the headers.
    firstLine = find(allData(:,1) > 0.1, 1); % Everything before ~100 ms is noise from the lock-in filter calculation; it sounds like this is default in the correction we get wqhen we extract DF/F0
    data = allData(firstLine:end, :);
    
    % Strip the last row of data if it has NaN, as it otherwise causes
    % problems
    lastrow = data(size(data, 1), :);
    if isnan(lastrow(2))
        data = data(1:size(data, 1) - 1, :);
    end
    
    %prelim graphing
    if PRELIM_GRAPHS
        if VISIBLE_GRAPHS
            figure
        else
            figure('Visible', 'off')
        end
        
        plot(data(:,1),data(:,2))
        title([filename(1:7) ' isosbestic'])
        print([FP_PRELIM_DIRECTORY slash filename(1:7) ' isosbestic'],'-dpng')
        
        if VISIBLE_GRAPHS
            figure
        else
            figure('Visible', 'off')
        end
        plot(data(:,1),data(:,3))
        title([filename(1:7) ' signal'])
        print([FP_PRELIM_DIRECTORY slash filename(1:7) ' signal'],'-dpng')
        
        close all
    end
    
    DF_F0 = calculateDF_F0_2nd_order(data);
    
    %Comment out for ZA
    %find index of spaces in filename
    
    %         spaces = strfind(filename,' ');
    %         if any(contains(filename(1:spaces(1)), redAndGreenMice)) %if contains red channel
    %             DIO = data(:,7);
    %             df_f0_red = calculateDF_F0_Red_2nd_order(data);
    %             currMouseRed = true;
    %         else
    %            currMouseRed = false;
    %            DIO = data(:, 5);
    %         end
    
    currMouseRed = false;
    DIO = data(:, 5);
    
    %% Trim here- Trimming is happening earlier in pipeline than in 2020 eLife Paper which will result in
    %%slightly different df/f0 calculation
    
    %ZA don't trim at all
    
    %         %find the start pulse (first 0 in DIO), align t = 0 to that, and trim
    %         %tempdata to remove pre-start pulse rows
    %         if autoTrim
    %             firstIndex = find(DIO<1, 1, 'first'); %trim to medpc start pulse
    %         else
    %             firstIndex = floor(trimTime*(NATIVE_SAMPLING_RATE));
    %         end
    %         starttime = DF_F0(firstIndex,1);
    %         DF_F0(:,1) = DF_F0(:,1)-starttime;
    %         DF_F0 = DF_F0(firstIndex:end,:);
    %         if currMouseRed
    %                 df_f0_red = df_f0_red(firstIndex:end,:);
    %                 df_f0_red(:,1) = df_f0_red(:,1)-starttime;
    %         end
    %         DIO = DIO(firstIndex:end,:);
    
    %%
    if DOWNSAMPLE_FP_DATA %thanks Alexa
        tt = DF_F0(end,1); %total time from begining to end (after trimming)
        t = tt * SAMPLING_RATE; %multiply by points per second to be a multiple of the frame rate
        newTime = linspace(0, tt, t); %make the new time vector
        newDIO = zeros(length(newTime),size(DIO,2)); %probably can remove
        newDF_F0 = zeros(length(newTime),size(DF_F0,2)); %probably can remove
        if currMouseRed
            newdf_f0_red = zeros(length(newTime),size(df_f0_red,2));
            for i = 1:size(df_f0_red,2)
                newdf_f0_red(:,i) = interp1(df_f0_red(:,1),df_f0_red(:,i), newTime); %downsample the data to match the frame rate
            end
            df_f0_red = newdf_f0_red;
        end
        for i = 1:size(DF_F0,2)
            newDF_F0(:,i) = interp1(DF_F0(:,1),DF_F0(:,i), newTime); %downsample the data to match the frame rate
        end
        newDIO = interp1(DF_F0(:,1),DIO, newTime);
        DIO = newDIO.';
        DF_F0 = newDF_F0;
        
    end
    
    
    if printExcel
        make_directory([output_directory ' excel']);
    end
    if currMouseRed
        correctedSignal = subtractReferenceAndSave(horzcat(DF_F0, df_f0_red(:,2)), output_directory, filename, DIO, printExcel, true,slash);
    else
        correctedSignal = subtractReferenceAndSave(DF_F0, output_directory, filename, DIO, printExcel, false, slash);
    end
    
    fprintf('Processed %s\n', filename);
end
Basic_FP_processing = true; %return true at end because why not
end