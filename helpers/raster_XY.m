function [X_single,Y_single , X_global ] = raster_XY(X,Y,dx,method)
    %RASTER_XY function to raster Y to a equidistance dx grid of X

    % get the lokal limits
    X_min_l = cellfun(@min,X);
    X_max_l = cellfun(@max,X);

    % get the lokal limits
    X_min_g = min(X_min_l);
    X_max_g = max(X_max_l);

    % get the number of grid points limited to less than 10^4
    n   = min(ceil((X_max_g-X_min_g)/dx),10^4); % limit

    % define the global grid
    X_global = linspace(X_min_g,X_max_g,n)';

    % get the single grid of each data set
    X_single = arrayfun(@(Xmax,Xmin) X_global( Xmin <= X_global & X_global <= Xmax),X_max_l,X_min_l,'UniformOutput', false);

    % eleminate data sets with less than two points
    tmp = cellfun(@(hlin) size(hlin,1) > 2 ,X_single);
    X = X(tmp);
    Y = Y(tmp);
    X_single = X_single(tmp);

    Y_single = cellfun(@(h, v, hlin)                 ...
                        interp1(h,v,hlin,method),    ...
                        X,                           ...
                        Y,                           ...
                        X_single,                    ...
                        'UniformOutput', false);
end