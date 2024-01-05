
% Makes a directory specified by MDIR_DIRECTORY_NAME with some basic error
% handling
function myMakeDir = make_directory(myDir)
    [~, ~, msg_id] =  mkdir(myDir);
    dir_msg_id = msg_id;
    if ~contains(dir_msg_id, 'DirectoryExists')
        fprintf('Created directory %s\n', myDir);
    end
    myMakeDir = myDir;
end