function myFindmouseNumbers = findMouseNumbers(filenames)
    filename = filenames{1};
    spaces = strfind(filename,' ');
    spaces = spaces(1);
    num = filename - '0' ;
    firstNum = find( num >= 0 & num < 10, 1) ;
    mouseNumbers = str2double(filename(firstNum:spaces(1)-1));
    
    for i = 2:size(filenames, 1)
        if isscalar(mouseNumbers)
            s = num2str(mouseNumbers);
        else
            s = num2str(mouseNumbers(length(mouseNumbers)));
        end
            filename = filenames{i};
            if not(contains(filename,s))
                spaces = strfind(filename,' ');
                spaces = spaces(1);
                %find first num
                num = filename - '0' ;
                firstNum = find( num >= 0 & num < 10, 1) ;

                mouseNumbers = [mouseNumbers, str2double(filename(firstNum:spaces(1)-1))];
            end
    end
    myFindmouseNumbers = mouseNumbers;
end