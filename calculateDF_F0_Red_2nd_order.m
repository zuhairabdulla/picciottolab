function df_f0_rcamp = calculateDF_F0_Red_2nd_order(rawData)
    time = rawData(:, 1);
    df_f0_rcamp = time;
    dataColIndexes = [5]; % Depending how it has been saved from DNS, the first number should be the column index of the reference, and the second one the index of the calcium signal. Change as needed.
     %modified loop to just be one iteration and only calculate df/f for
     %rcamp signal
    for colIdx = dataColIndexes(1) %actual calculation of DF/F0
        data = rawData(:, colIdx);
        reg = polyfit(time, data, 2);
        f0 = reg(1).*time.*time + reg(2)*time + reg(3);
        df_f0_tmp = (data - f0)./f0 * 100;
        df_f0_rcamp = horzcat(df_f0_rcamp, df_f0_tmp);
      end
         
end