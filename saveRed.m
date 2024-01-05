function myData = saveRed(df_f0, directory, filename, DIO, saveExcel)
        myData = horzcat(df_f0,DIO);
        filename_w = strcat(directory,'\','PROCESSED_', filename); %so basically all it does is to take the raw data and saves the corrected signal with a filename preceded by "processed_" (to avoid overwriting the original file... and the other function calls for that
        cHeader = {'Time' 'Reference (DF/F0)' 'Ca2+ Signal (DF/F0)' 'DIO'};

      
          filename_w = filename_w(1:end-4);
          save(filename_w, 'myData', 'cHeader');
      if saveExcel
          directory = [directory ' excel'];
          filename_w = strcat(directory,'\','PROCESSED_', filename); %so basically all it does is to take the raw data and saves the corrected signal with a filename preceded by "processed_" (to avoid overwriting the original file... and the other function calls for that

          commaHeader = [cHeader;repmat({','},1,numel(cHeader))];
          commaHeader = commaHeader(:)';
          textHeader = cell2mat(commaHeader);
          textHeader = textHeader(1:end-1);
          fid = fopen(filename_w,'w'); 
          fprintf(fid,'%s\n',textHeader);
          fclose(fid);
          dlmwrite(filename_w, myData, 'delimiter', ',', '-append', 'precision','%.5f');
      end

end