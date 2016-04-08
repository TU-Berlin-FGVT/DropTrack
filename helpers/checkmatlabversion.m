function out = checkmatlabversion(matlabversion_required, operator)
%CHECKMALTLABVERSION checks the required matlabversion for downgrade compatibility functions.
%
% matlabversion is a formated string e.g. '2014a' or '2014', given by version('-release')
%
% By default this functions returns true if the required and system matlabversion is the same.
% operator is optional and contains one of  following stings:
%
%    '<','>','==','<=','>='
%

if ~exist('operator','var')
    operator = '==';
end

version_required = matlabversion2double(matlabversion_required);
version_current = matlabversion2double();

switch operator
    case '<'
        out = version_current <  version_required;
    case '>'
        out = version_current >  version_required;
    case '=='
        out = version_current ==  version_required;
    case '<='
        out = version_current <=  version_required;
    case '>='
        out =  version_current >=  version_is;
    otherwise
        error('checkmatlabversion:InvalidInputArgument:operator','operator expects one of `<`,`>`,`==`.');
end
