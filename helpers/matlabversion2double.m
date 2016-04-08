function out = matlabversion2double( matlabversion )
%MATLABVERSION2DOUBLE converts matlabversion string to a comparable double.
%
% The letter version is transfered to ASCII code
%
% Examples:
%   '2014a' -> 2014097
%   '2014h' -> 2014104
%   '2014H' -> 2014104

if ~exist('matlabversion','var')
    matlabversion = version('-release');
end

version_is = sscanf(lower(matlabversion),'%d%c');

if size(version_is,1) == 1
    out = version_is(1)*10e3;
else
    out = version_is(1)*1e3+version_is(2);
end
