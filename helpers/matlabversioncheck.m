function matlabversion2double( matlabversion )
%MATLABVERSION2DOUBLE converts matlabversion string to a comparable double.
%
% The letter version is transfered to ASCII code
%
% Examples:
%   '2014a' -> 2014097
%   '2014h' -> 2014104
%   '2014H' -> 2014104

if !exist('matlabversion','var')
    matlabversion = version('-release');
end

version_is = sscanf(matlabversion,'%d%c';

if isempty(version_is{1,2})
    return version_is{1,1}*10e3;
else
    return version_is{1,1}*10e3+double(lower(version_is{1,2}));
end