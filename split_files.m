function Basic_FP_processing = split_files(directory,output_directory,slash)
   
load(getPipelineVarsFilename); %load pipeline variables


    files = dir(directory);
    files = files(contains({files.name},{'.csv'}));
   
    processed_files = dir(output_directory);
    
    %remove non csv files
    processed_files = [ processed_files(:).name ]; 

    make_directory(output_directory); %RIC i fixed makedirectory to be cooler


    for file = files'

        filename = strcat(file.name);
        %only process .csv files, don't process "PROCESSED" files, and don't
        %process any that already have a 'PROCESSED' version in the folder
        if isempty(strfind(filename, '.csv')) || isempty(strfind(filename, 'PROCESSED_'))==false || sum(strcmp(strcat('PROCESSED_',filename),{files.name}))>0 || ~isempty(strfind(processed_files, filename)) || ( COMPILE_WITH_REF  && contains(filename, '0850') ) 
            fprintf('Skipping %s\n', filename);
            continue
        end

        allData = readmatrix([directory,slash filename]); % Just read the numbers here, not headers, for speed. Will convert to cell right before writing
        allDataHeader = readcell([directory,slash filename],'Range','A1:Z2'); %Just read headers, with extra columns, just in case. 
        
        %grabbing based on a 3 digit number, may need to change if mouse
        %numbers vary
        mouse_one = filename(1:3);
        mouse_two = filename(5:7);
        base_filename = filename(8:end);
        
        %slow but works
%         %Grab the columns for data, then header, that correspond to each mouse. Could make this read
%         %through headers to find the ones we want, but hard coded here for
%         %simplicity
%         %Then, join and save
%         mouse_one_data = num2cell(allData(:,[1 2 3 4 8 10 11]));
%         mouse_one_header = allDataHeader(:,[1 2 3 4 8 10 11]);
%         mouse_one_join = [mouse_one_header ; mouse_one_data];
%         writecell(mouse_one_join, [output_directory '\' mouse_one base_filename]);
%         
%         mouse_two_data = num2cell(allData(:,[1 5 6 7 9 12 13]));
%         mouse_two_header = allDataHeader(:,[1 5 6 7 9 12 13]);
%         mouse_two_join = [mouse_two_header; mouse_two_data];
%         writecell(mouse_two_join, [output_directory '\' mouse_two base_filename]);
        
        
                %Grab the columns for data, then header, that correspond to each mouse. Could make this read
        %through headers to find the ones we want, but hard coded here for
        %simplicity
        %Then, join and save
        mouse_one_header = allDataHeader(:,[1 2 3 4 8 10 11]);
        mouse_one_data = allData(:,[1 2 3 4 8 10 11]);
        writecell(mouse_one_header, [output_directory slash mouse_one base_filename]);
        writematrix(mouse_one_data, [output_directory slash mouse_one base_filename],'WriteMode','append');
        
        mouse_two_header = allDataHeader(:,[1 5 6 7 9 12 13]);
        mouse_two_data = allData(:,[1 5 6 7 9 12 13]);
        writecell(mouse_two_header, [output_directory slash mouse_two base_filename]);
        writematrix(mouse_two_data, [output_directory slash mouse_two base_filename],'WriteMode','append');
        
        
        fprintf('Split %s\n', filename);
    end
     Split_Files = true; %return true at end because why not
end