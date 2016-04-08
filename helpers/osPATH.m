function path = osPATH( path )
%OSPATH returns the correct path of the operation system

    path =  strrep(path,'\',filesep);
    path =  strrep(path,'/',filesep);
end

