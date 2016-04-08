function path = correct_path_ending(path)
    if ispc && path(end) ~='\' % for windows
         path = [path,'\'];
    elseif ~ispc && path(end) ~='/' % for Linux & Mac
        path = [path,'/'];
    end
end