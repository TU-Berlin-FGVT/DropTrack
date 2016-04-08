function out = get_interval(lin)
%GET_INTERVAL to get the interval limits of 1-dim binary cluster
%
% example:
%
% lin = [0;1;1;1;1;0;0;0;0;1;1;0;0;1;0];
% out = get_interval(lin);
% disp(out)
% out =
%   [2,5 ; 10,11; 14;14 ];
%
%
    if ~isempty(lin)
        lin = logical(lin);
        a = diff(lin);
        index = 1:length(lin);
        b1_i = index(a ==1)+1;
        b2_i = index(a ==-1);
        if lin(1)
            b1_i =[1,b1_i];
        end
        if lin(end)
            b2_i =[b2_i, index(end)];
        end
        out = [b1_i', b2_i'];
    else
        out = [];
    end
end