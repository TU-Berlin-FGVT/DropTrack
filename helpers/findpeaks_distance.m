function [ pks, loc ] = findpeaks_distance(x,y,dx)
%FINDPEAKS_DISTANCE to find local maxima of y-values within defined distance along the x-direction
%
% input:
%   x  - double(n,1), values ordered by ASC
%   y  - double(n,1),
%   dx - double(1,1),
%
%  P(x,y) is element of Real^2
%
% output:
%  pks - local maxima of y-values
%  loc - binary list in which local maxima marked as true
%
    X =x(:);
    Y =y(:);
    LOC = false(size(X)); % marks local maxima as true

    % split values into chunks, if the distance along x is greater than 2*dx
    int = get_interval(diff(X) < 2*dx);

    if ~isempty(int)
        if int(1) ~= 1
            int = [1,1;int];
        end
        if int(end) ~= size(X,1)
            int = [int;size(X,1),size(X,1)];
        end

        % search for local max in every chunk
        for j = 1:size(int,1)

            % get the x,y value of the current chunk
            x = X(int(j,1):int(j,2));
            y = Y(int(j,1):int(j,2));

            loc =  y == max(y); % marks local maxima of the current chunk as true

            % set init marker
            % .*_l -> marker for the search at descending x-value direction
            % .*_h -> marker for the search at ascending x-value direction
            i_l = find(loc,1,'first');
            i_h = i_l;

            % set the limit of the current chunk;
            x_max = max(x);
            x_min = min(x);

            % search for local maxima in current chunk
            while x_min+dx < x(i_l) ||  x(i_h) < x_max-dx

                % search local maximum in descending direction within distance of 2*dx
                tmp_l = x(i_l) - 2*dx <= x  & x < x(i_l);
                if any(tmp_l)
                    % find and add local maximum to the list of the current chunk
                    loc = loc |  ( tmp_l  & y == max( y(tmp_l) ) );
                end

                % search local maxima in ascending direction within distance of 2*dx
                tmp_h = x(i_h)        <  x  & x <= x(i_h) + 2*dx ;
                if any(tmp_h)
                    % find and add local maximum to the list of the current chunk
                    loc = loc |  ( tmp_h  & y == max( y(tmp_h) ) );
                end

                % break the loop if no maximum can be found.
                if ~any(tmp_l) && ~ any(tmp_h)
                    break;
                end

                % preparation for in next loop
                i_l = find(loc,1,'first');
                i_h = i_h+find(loc((i_h+1):end),1,'first');
                if ~any(i_l)
                    i_l = 1;
                end
                if ~any(i_h)
                    i_h = size(loc,1);
                end
            end

            % transfer marked local maxima of the current chunk to the global list
            LOC(int(j,1):int(j,2))= loc ;
        end
        pks = Y(LOC);
        lin = 1:size(LOC,1);
        loc = lin(LOC)';
    else
        pks = Y;
        loc = 1:size(Y,1);
    end
end
