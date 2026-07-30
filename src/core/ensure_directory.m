function ensure_directory(path_name)
%ENSURE_DIRECTORY Create directory if needed.
    if ~exist(path_name, 'dir')
        mkdir(path_name);
    end
end
