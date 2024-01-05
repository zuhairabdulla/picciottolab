function df_f0 = calculateDF_F0(rawData, order)
    time = rawData(:, 1);
    df_f0 = time;
    dataColIndexes = [2, 3]; % Depending how it has been saved from DNS, the first number should be the column index of the reference, and the second one the index of the calcium signal. Change as needed.
      for colIdx = dataColIndexes(1):dataColIndexes(2) %actual calculation of DF/F0
        data = rawData(:, colIdx);       
        if order == 1
            reg = polyfit(time, data, 1);
            f0 = reg(1)*time + reg(2);
        elseif order == 2
            reg = polyfit(time, data, 2);
            f0 = reg(1).*time.*time + reg(2)*time + reg(3);
        else
            error('if you need to do more than 2nd order your data is trash')
        end
        df_f0_tmp = (data - f0)./f0 * 100;
        df_f0 = horzcat(df_f0, df_f0_tmp);
      end
end