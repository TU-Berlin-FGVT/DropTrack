% DropTrack provides a complete pipeline from data management, image preprocessing,
% drop detection and tracking, event detection, to parameter estimation
% (e.g. velocity)
%
% For details please visit http://www.rhaensch.de/droptrack.html
%
% For questions or remarks, don't hesitate to contact:
% johannes.kamp@tu-berlin.de
%
% If you use this code in your work, please refer to the following paper:
% Kamp, J.; Hänsch, R.; Kendzierski, G.; Kraume, M. & Hellwich, O. 
% Automated image analysis for trajectory determination of single drop collisions 
% Computers & Chemical Engineering, 2016
% http://dx.doi.org/10.1016/j.compchemeng.2016.03.033
%
% DropTrack (c) 2015, Johannes Kamp, Ronny Hänsch, Gregor Kendzierski and contributors.
%
% DropTrack is free software: you can redistribute it and/or modify it under the terms
% of the BSD 2-clause License.
%
% DropTrack is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
% PURPOSE. See the license file provided with the code for the license terms and more
% details.
%
% Created by Gregor Kendzierski, December 2015
%

classdef view_SEQ < handle
    %VIEW_SEQ class is used to generate averaged plots of records with same parameters of rising drops

    properties ( Access = public )
        path   % path in which the records with same parameter can be found

        DATA      % cell array of single records' data
        DATA_path % cell array of paths for each record
        dataFileName = 'exp_data.mat';   % name of the records' datafile

        SIunits = false; % Output Units: false -> image; true -> SI-Units
    end

    methods ( Access = public )
        function obj = view_SEQ(path)
            %VIEW_SEQ constructor of the view_SEQ class

            if ~exist(path,'dir'), error('record_EXP:FalsePath','path does not exist!'); end
            obj.path = correct_path_ending(osPATH(path));
        end

        function out = load_data(obj)
            %LOAD_DATA method to load all data and reduce to required data

            obj.DATA_path = getAllFiles(obj.path,[obj.dataFileName,'$'],true);

            obj.DATA = cell(size(obj.DATA_path,1),1);

            % load data
            for i = 1:size(obj.DATA_path)

                load(obj.DATA_path{i},'data'); % load single data file

                % reduce to required only data
                f = fieldnames(data);
                f = f(~ismember(f,{'drops','parameters','frames','Events','dateTime','commitID','commitCLEAN','error','log'}));
                for i_ = 1:size(f,1)
                    data = rmfield(data,f{i_});
                end

                f = fieldnames(data.drops);
                f = f(~ismember(f,{'centers','velocities','validates'}));
                for i_ = 1:size(f,1)
                    data.drops = rmfield(data.drops,f{i_});
                end
                data.drops(1).centers = [];
                data.drops(1).velocities = [];
                obj.DATA{i} = data;
            end

            out = obj.DATA;
        end

        function out = rm_errorDATA(obj,arraylabel)
            %RM_ERRORDATA method to return labeling array to unmask records with thrown errors
            % out =-1 -> noTemperature is given
            % out = 0 -> not selected parameterSet
            % out > 0 -> same temperatures of selected records have the same numbered label

            if ~exist('arraylabel','var'), arraylabel = true(size(obj.DATA,1),1); end

            out(arraylabel) = cellfun(@(d) isempty(d.error), obj.DATA(arraylabel));
        end

        function [out, parameterSetValues] = label_same_parameterSet(obj,arraylabel)
            %LABEL_SAME_PARAMETERSET method to return labeling array to mask records with same parameterSets
            % out = 0 -> not selected records
            % out > 0 -> same paramterSets have the same numbered label

            if ~exist('arraylabel','var'), arraylabel = true(size(obj.DATA,1),1); end
            para = cell2mat(cellfun(@(d) [d.parameters.d_disp_bot_target,  ...
                                          d.parameters.d_disp_top_target,  ...
                                          d.parameters.drop_edge_distance  ...
                                         ], obj.DATA(arraylabel),'UniformOutput',false));
            parameterSet = zeros(size(obj.DATA,1),1);
            [p ,label,parameterSet(arraylabel)] = unique(para,'rows');

            parameterSetValues = cell2struct(                               ...
                [{label} ; mat2cell(p,size(p,1),ones(size(p,2),1))'],       ...
                {'label'; 'd_oil_bottom'; 'd_oil_top';'drop_edge_distance'} ...
                );

            out = parameterSet;
        end
        function out = label_same_temperature(obj,arraylabel)
            %LABEL_SAME_TEMPERATURE method to return labeling array to mask records with temperature
            % out =-1 -> noTemperature is given
            % out = 0 -> not selected parameterSet
            % out > 0 -> same temperatures of selected records have the same numbered label

            out = double(any(arraylabel,2));

            % select paramterSet with given gitVersions
            T_tmp = cellfun(@(d)isfield(d.parameters,'temperature'),obj.DATA(arraylabel));
            out(out==1) = -1*(~T_tmp);

            if any(T_tmp) && any(~T_tmp), warning('view_SEQ:check_same_temperature:noTemperature','There are following record with no given temperature:\n %s', strcat(obj.DATA_path(~ID_tmp),'\n')), end

            if any(T_tmp)
                % group parameterSet with same temperature by labeling
                [~,~,ci] = unique(cellfun(@(d)d.parameters.temperature,obj.DATA(0<out)));
                out(0<out) = ci;
            end
        end
        function out = label_same_gitVersion(obj,arraylabel)
            %LABEL_SAME_GITVERSION method to return labeling array to mask records with same gitVersions
            % out =-1 -> no gitVersion is given
            % out = 0 -> not selected parameterSet
            % out > 0 -> same gitVersions of selected records have the same numbered label

            out = double(any(arraylabel,2));

            % select paramterSet with given gitVersions
            ID_tmp = ~cellfun(@(d)isempty(d.commitID),obj.DATA(out));
            out(out==1) = -1*(~ID_tmp);

            if any(ID_tmp) && any(~ID_tmp), warning('view_SEQ:check_same_gitVersion:noGitVersions','There are following record with no given gitVersion:\n %s', strcat(obj.DATA_path(~ID_tmp),'\n')), end

            if any(ID_tmp)
                % group parameterSet with same gitVersion by labeling
                [~,~,ci] = unique(cellfun(@(d)d.commitID,obj.DATA(0<out)));
                out(0<out) = ci;
            end
        end

        function set_SIunits(obj, value)
            %SET_SIUNITS method to set output units to SI-Units or Image-Units
            % value = true -> SI-Units
            % value = false-> Image-Units

            obj.SIunits = any(value);
        end
        function [px2mm, fps] = get_mean_SIunitsFactors(obj,arraylabel)
            %GET_MEAN_SIUNITSFACTORS method to estimate the converting factors to SI-Units

            px2mm = 1;
            fps   = 1;
            if obj.SIunits
                fps = cellfun(@(d) d.parameters.Framerate, obj.DATA(arraylabel),'UniformOutput',false);
                fps = unique(cell2mat(fps(~cellfun(@isempty,fps))));
                if size(fps,1)~=1, error('view_SEQ:get_mean_SIunitsFactors:diffFPS','Different framerates in given arraylabel!'); end

                px2mm = cellfun(@(d) d.frames.px2mm, obj.DATA(arraylabel),'UniformOutput',false);
                px2mm = mean(cell2mat(px2mm(~cellfun(@isempty,px2mm))));
            end
        end

%% PLOTTING methods
        function [fh, data] = subplot_averaged(obj, arraylabel,varargin)
            %SUBPLOT_AVERAGED method to plot the averaged data of selected records
            % ah - axis_handle of the plot
            % See:  plot_averaged_data(), plot_averaged_height_over_time(),
            %       plot_averaged_risingVelocity_over_time, plot_averaged_risingVelocity_over_height
            fh = figure();
            n = 3;
            data = struct();
            try
                p3 = subplot(3,1,3);
                [~, data_avVH] = obj.plot_averaged_risingVelocity_over_height(arraylabel,p3,varargin);
                data.averaged_risingVelocity_over_height = data_avVH;
            catch e
                if strcmp('view_SEQ:get_strict_monotone_height:noData',e.identifier)
                    delete(p3);
                    n = 2;
                end
                rethrow(e);
            end

            p1 = subplot(n,1,1);
            varargin = ['LegendOFF','warningOFF',varargin];
            [~, data_avHT] = obj.plot_averaged_height_over_time(arraylabel,p1,varargin);
            data.averaged_height_over_time = data_avHT;

            p2 = subplot(n,1,2);
            [~, data_avVT] = obj.plot_averaged_risingVelocity_over_time(arraylabel,p2,varargin);
            data.averaged_risingVelocity_over_time = data_avVT;

            if exist('e','var'), rethrow(e); end;
        end

        function ah = plot_averaged_data( ~ , ah, data, varargin)
            %PLOT_AVERAGED_DATA method to plot averaged data generally
            % ah - axis_handle of the plot
            % 'OverwriteOFF' - switch off to overwrite the plot, just adding to the plot
            % 'LegendOFF' - switch off the the legend
            % 'NumberOfRecordOFF' - switch off the number of record which is used for averaging of each timestep
            % See:  plot_averaged_height_over_time(), plot_averaged_risingVelocity_over_time,
            %       plot_averaged_risingVelocity_over_height

            if ~exist('ah','var')
                ah = gca(figure());
            end
            if iscell(varargin{:})
                varargin = varargin{:}{:};
            end
            if ismember('OverwriteOFF', varargin), hold on; else hold off; end
            Ylimits = [min([data.Y_min;data.Y_mean-data.Y_std]); max([data.Y_max;data.Y_mean+data.Y_std])];
            if ismember('NumberOfRecordOFF', varargin)
                ph_mean = plot(ah, data.X,data.Y_mean,'-b');
                ylim(ah, Ylimits);
                ph_number = [];
            else
                [ax, ph_mean, ph_number] =  plotyy(ah, data.X,data.Y_mean,data.X ,data.numberOfRecords );
                ylim(ax(1), Ylimits);
                set(ax(1),'ycolor','k');

                % second right axes
                ylabel(ax(2),'number of records');
                yRightLimits = [min(data.numberOfRecords),max(data.numberOfRecords)];
                yRightTicks = linspace(yRightLimits(1),yRightLimits(2),3);
                yRightTicksLabel = arrayfun(@num2str, yRightTicks,'UniformOutput',false);
                set(ax(2),'YTick',yRightTicks,'YTickLabel', yRightTicksLabel);
                ylim(ax(2),yRightLimits+[0,1]);
                set(ph_mean, 'LineStyle','-','Color','b');
                ph_number = {ph_number, 'averaged number of record'};
            end
            ph_mean =   {ph_mean, 'average'};

            title(data.title)

            hold on;
                ph_maxmin = {plot(ah,data.X,data.Y_max,'-k'),'max/min'};
                            plot(ah,data.X,data.Y_min,'-k');
                ph_std    = {plot(ah,data.X,data.Y_mean+data.Y_std,':r'),'+/- std'};
                            plot(ah,data.X,data.Y_mean-data.Y_std,':r');
                extrema =  [max([data.Y_mean; data.Y_mean+data.Y_std]); min([data.Y_min; data.Y_mean-data.Y_std])];
                ph_bott = [];
                if isfield(data, 'bottomDetachmentEvent')
                    ph_bott   = {plot(ah,data.bottomDetachmentEvent*ones(2,1), extrema ,'--r'),'averaged bottomDetachment event'};
                end
                ph_coll = [];
                if isfield(data, 'collisionEvent')
                    ph_coll   = {plot(ah,data.collisionEvent*ones(2,1), extrema ,'-.k'),'averaged collision event'};
                end

            hold off;

            if  ~ismember('LegendOFF', varargin)
                ph_legend = [ph_mean; ph_maxmin; ph_std; ph_bott; ph_coll; ph_number];

                ph_legend_linehandle = ph_legend(:,1);
                legend([ph_legend_linehandle{:}]', ph_legend(:,2));
            end
        end
        function [ah, data] = plot_averaged_height_over_time(obj,arraylabel, ah, varargin)
            %PLOT_AVERAGED_HEIGH_OVER_TIME method to plot the averaged height over time
            % ah - axis_handle of the plot
            % See:  plot_averaged_data()

            obj.plot_warnings(arraylabel, varargin);

            data = obj.get_averaged_height_over_time(arraylabel);

            data.title = 'Averaged height / time of the rising droplet';

            ah = obj.plot_averaged_data(ah, data, varargin);

            if obj.SIunits
                data.xlabel = 't [s]';
                data.ylabel = 'x [mm]';
            else
                data.xlabel = 't [frame]';
                data.ylabel = 'x [px]';
            end
            xlabel(data.xlabel);
            ylabel(data.ylabel);
        end
        function [ah, data] = plot_averaged_risingVelocity_over_time(obj,arraylabel, ah, varargin)
            %PLOT_AVERAGED_RISINGVELOCITY_OVER_TIME method to plot the averaged risingVelocity over time
            % ah - axis_handle of the plot
            % See:  plot_averaged_data()

            obj.plot_warnings(arraylabel, varargin);

            data = obj.get_averaged_risingVelocity_over_time(arraylabel);

            data.title = 'Averaged risingVelocity / time of the rising droplet';

            ah = obj.plot_averaged_data(ah, data, varargin);

            if obj.SIunits
                data.xlabel = 't [s]';
                data.ylabel = 'dx/dt [mm/s]';
            else
                data.xlabel = 't [frame]';
                data.ylabel = 'dx/dt [px/frame]';
            end
            xlabel(data.xlabel);
            ylabel(data.ylabel);
        end
        function [ah, data] = plot_averaged_risingVelocity_over_height(obj,arraylabel, ah, varargin)
            %PLOT_AVERAGED_RISINGVELOCITY_OVER_HEIGHT method to plot the averaged risingVelocity over height
            % ah - axis_handle of the plot
            % See:  plot_averaged_data()

            obj.plot_warnings(arraylabel, varargin);

            data = obj.get_averaged_risingVelocity_over_height(arraylabel);

            data.title = 'Averaged risingVelocity / height of the rising droplet';

            ah = obj.plot_averaged_data(ah, data, varargin);

            if obj.SIunits
                data.xlabel = 'x [mm]';
                data.ylabel = 'dx/dt [mm/s]';
            else
                data.xlabel = 'x [px]';
                data.ylabel = 'dx/dt [px/frame]';
            end
            xlabel(data.xlabel);
            ylabel(data.ylabel);
        end

%% EVENT methods
        function out = get_bottomDetachmentTime(obj, arraylabel)
            %GET_BOTTOMDETACHMENTTIME method to get the mean of the bottomDetachmentTime intervall for each selected record

            out = cellfun(@(d) round(mean(d.Events.bottomDetachmentTime)),obj.DATA(arraylabel),'UniformOutput',false);
        end
        function [timeshift, arraylabel] = get_timeshift_bottomDetachmentEvent(obj,arraylabel )
            %GET_TIMESHIFT_BOTTOMDETACHMENTEVENT method to get the alginement in time with the bottomDetachmentEvent for each selected record, records with no detected event will be ignored

            % get the bottomDetachmentTime of each selected record
            timeshift = obj.get_bottomDetachmentTime(arraylabel);

            % eleminate records with no detected bottomDetachmentEvent
            tmp = cellfun(@(t) 1 < t, timeshift);
            timeshift = timeshift(tmp);
            arraylabel(arraylabel) = tmp;
        end
        function out = get_averaged_bottomDetachmentHeight(obj,arraylabel, heightshift)
            %GET_AVERAGED_BOTTOMDETACHMENTHEIGHT method to get the averaged bottomDetachmentHeight of selected records

            % get time of bottomDetachmentTime for each record
            bottomDetachmentTime  = obj.get_bottomDetachmentTime(arraylabel);

            % reduce records with no detected collisionEvent
            bottomDetachmentTime = bottomDetachmentTime(cellfun(@(d) 1 <d, bottomDetachmentTime));

            if ~exist('heightshift','var'), heightshift = 0; end
            if iscell(heightshift), heightshift = cell2mat(heightshift); end

            out = mean( cellfun(@(d,t) real(d.drops(2).centers(t)), obj.DATA(arraylabel),bottomDetachmentTime) - heightshift);
        end

        function out = get_collisionTime(obj,arraylabel)
            %GET_COLLISIONTIME method to get the mean of the collisionTime intervall for each selected record

            out = cellfun(@(d) round(mean(d.Events.collisionTime)),obj.DATA(arraylabel),'UniformOutput',false);
        end
        function out = get_averaged_collisionTime(obj,arraylabel, timeshift)
            %GET_AVERAGED_COLLISIONTIME method to get the averaged collisionTime of selected records

            % get time of collisionEvent for each record
            collisionTime  = cell2mat(obj.get_collisionTime(arraylabel));

            % reduce records with no detected collisionEvent
            collisionTime = collisionTime(1 < collisionTime);

            if ~exist('timeshift','var'), timeshift = 0; end
            if iscell(timeshift), timeshift = cell2mat(timeshift); end

            % shift the colloisionEvent time for each record in the same way as above and
            % get the average of this
            out = mean(collisionTime - timeshift);
        end
        function out = get_averaged_collisionHeight(obj,arraylabel, heightshift)
            %GET_AVERAGED_COLLISIONHEIGHT method to get the averaged collisionHeight of selected records

            % get time of collisionEvent for each record
            collisionTime  = obj.get_collisionTime(arraylabel);

            % reduce records with no detected collisionEvent
            collisionTime = collisionTime(cellfun(@(d) 1 <d, collisionTime));

            if ~exist('heightshift','var'), heightshift = 0; end
            if iscell(heightshift), heightshift = cell2mat(heightshift); end

            out = mean( cellfun(@(d,t) real(d.drops(2).centers(t)), obj.DATA(arraylabel),collisionTime) - heightshift);
        end

%% AVERAGING DATA methods
        function out = get_averaged_height_over_time(obj,arraylabel)
            %GET_AVERAGED_HEIGHT_OVER_TIME method to average height of seleced records for each timestep

            % get the alignment in time with the bottomDetachmentEvent for each selected record
            [t_min, arraylabel] = obj.get_timeshift_bottomDetachmentEvent(arraylabel);

            % get heigth and convert every record to uniform length by adding NaNs
            t_max = max(cell2mat(cellfun(@(d) size(d.drops(2).centers,1),obj.DATA(arraylabel),'UniformOutput',false)));
            height = cell2mat(cellfun(@(d,t)                                ...
                        [ real(d.drops(2).centers(t:end));                  ...
                          NaN(t_max-size(d.drops(2).centers(t:end),1),1)    ...
                        ],                                                  ...
                        obj.DATA(arraylabel),t_min,'UniformOutput',false)');

            % Reduce rows with only NaN
            t = (1:find(sum(~isnan(height),2)> 1,1,'last'))';
            height = height(t,:);

            height_mean = arrayfun(@(t) mean(height(t,~isnan(height(t,:)))),t);
            height_std  = arrayfun(@(t) std (height(t,~isnan(height(t,:)))),t);
            height_max  = arrayfun(@(t) max (height(t,~isnan(height(t,:)))),t);
            height_min  = arrayfun(@(t) min (height(t,~isnan(height(t,:)))),t);

            % get the averaged collisionTime of selected records
            collisionEvent = obj.get_averaged_collisionTime(arraylabel, t_min);

            % get the number of record which are used for averaging of each timestep
            numberOfRecords = sum(~isnan(height),2);

            [px2mm, fps] = obj.get_mean_SIunitsFactors(arraylabel);

            out = struct( ...
                'X' ,            ((1:size(t,1))'-1) / fps, ...
                'Y_mean',         height_mean * px2mm,           ...
                'Y_std',          height_std  * px2mm,           ...
                'Y_max',          height_max  * px2mm,           ...
                'Y_min',          height_min  * px2mm,           ...
                'collisionEvent', (collisionEvent-1) / fps,      ...
                'numberOfRecords',numberOfRecords,               ...
                'arraylabel',     arraylabel                     ...
                        );
        end
        function out = get_averaged_risingVelocity_over_time(obj,arraylabel)
            %GET_AVERAGED_RISINGVELOCITY_OVER_TIME method to average risingVelocity of selected records for each timestep

            % get the alignment in time with the bottomDetachmentEvent for each selected record
            [t_min, arraylabel] = obj.get_timeshift_bottomDetachmentEvent(arraylabel);

            % get velocity and convert every record to uniform length by adding NaNs
            t_max = max(cell2mat(cellfun(@(d) size(d.drops(2).centers,1),obj.DATA(arraylabel),'UniformOutput',false)));
            risingVelocity = cell2mat(cellfun(@(d,t)                                 ...
                                [ real(d.drops(2).velocities(t:end));                ...
                                  NaN(t_max-size(d.drops(2).velocities(t:end),1),1)  ...
                                ],                                                   ...
                                obj.DATA(arraylabel),t_min,'UniformOutput',false)');

            % Reduce rows with only NaN
            t = (1:find(sum(~isnan(risingVelocity),2)> 1,1,'last'))';
            risingVelocity = risingVelocity(t,:);

            risingVelocity_mean = arrayfun(@(t) mean(risingVelocity(t,~isnan(risingVelocity(t,:)))),t);
            risingVelocity_std  = arrayfun(@(t) std (risingVelocity(t,~isnan(risingVelocity(t,:)))),t);
            risingVelocity_max  = arrayfun(@(t) max (risingVelocity(t,~isnan(risingVelocity(t,:)))),t);
            risingVelocity_min  = arrayfun(@(t) min (risingVelocity(t,~isnan(risingVelocity(t,:)))),t);

            % get the averaged collisionTime of selected records
            collisionEvent = obj.get_averaged_collisionTime(arraylabel, t_min);

            % get the number of record which are used for averaging of each timestep
            numberOfRecords = sum(~isnan(risingVelocity),2);

            [px2mm, fps] = obj.get_mean_SIunitsFactors(arraylabel);

            out = struct( ...
                'X' ,             ((1:size(t,1))'-1) / fps, ...
                'Y_mean',         risingVelocity_mean * px2mm * fps,   ...
                'Y_std',          risingVelocity_std  * px2mm * fps,   ...
                'Y_max',          risingVelocity_max  * px2mm * fps,   ...
                'Y_min',          risingVelocity_min  * px2mm * fps,   ...
                'collisionEvent', (collisionEvent-1) / fps,      ...
                'numberOfRecords',numberOfRecords,               ...
                'arraylabel',     arraylabel                     ...
            );
        end

        function [height, risingVelocity, arraylabel] = get_strict_monotone_height(obj,arraylabel)
            %GET_STRICT_MONOTONE_HEIGHT method to force strict monotone decreasing height by croping intervalls

            % get the alignment in time with the bottomDetachmentEvent for each selected record
            [t_min, arraylabel] = obj.get_timeshift_bottomDetachmentEvent(arraylabel);

            % get height and risingVelocity
            height         = cellfun(@(d,t) real(d.drops(2).centers(t:end))   ,obj.DATA(arraylabel),t_min,'UniformOutput',false);
            risingVelocity = cellfun(@(d,t) real(d.drops(2).velocities(t:end)),obj.DATA(arraylabel),t_min,'UniformOutput',false);

            log_index = arraylabel(arraylabel);
            for i = 1:size(height,1)
                % Note: due to the defined coordinate system the droplet height
                % decreases when the droplet in the image rises

                % get index where 0 <= d height ( droplet's center does not rise in image sequence )
                tmp = 0 <=diff(height{i});

                if ~any(tmp), continue; end % height is already strict monotone decreasing

                % get starting indices of intervalls where  0 <= d height
                index_start = get_interval(tmp)*[1;0] ;

                % get the index after each inverall where the height is
                % decreasing again ( droplet's center rises again )
                index_end = arrayfun(@(n)n + find(height{i}((n+1):end) < height{i}(n) ,1,'first'),index_start,'UniformOutput', false);
                index_end(cellfun(@isempty,index_end))={size(height{i},1)+1};

                % get the intervall limits which will be croped
                int = [index_start(:,1)+1,cell2mat(index_end)-1];

                if isempty(int)
                    log_index(i) = false;
                    continue;
                end
                % get rid off nested intervalls
                lin = unique(cell2mat(arrayfun(@(n1,n2) n1:n2, int(:,1),int(:,2),'UniformOutput',false)'));
                tmp = true(size(height{i}));
                tmp(lin) = false;

                height{i}         = height{i}(tmp);
                risingVelocity{i} = risingVelocity{i}(tmp);
            end

            height = height(log_index);
            risingVelocity = risingVelocity(log_index);

            arraylabel(arraylabel) = log_index;
            if ~any(arraylabel), error('view_SEQ:get_strict_monotone_height:noData','Strict monotone decreasing height can not be enforced. You may have to preprocessing this data'); end
        end
        function out = get_averaged_risingVelocity_over_height(obj,arraylabel)
            %GET_AVERAGED_RISINGVELOCITY_OVER_HEIGHT method to average risingVelocity
            %of selected records for each height which were converted to strict monotony

            % force strict monotone decreasing height by cropping intervalls
            [height, risingVelocity, arraylabel] = obj.get_strict_monotone_height(arraylabel);

            % raster data set to the same grid via linear interpolation
            dh_min = mean(cell2mat(cellfun(@(h) min(abs(diff(h))) ,height,'UniformOutput',false)));
            [height, risingVelocity, height_global ] = raster_XY(height, risingVelocity, dh_min, 'linear');

            H = cell2mat(height);
            V = cell2mat(risingVelocity);
            grouped_v = arrayfun(@(h) V(H==h),height_global,'UniformOutput',false);

            % remove empty groups
            tmp = ~cellfun(@isempty,grouped_v);
            height_global = height_global(tmp);
            grouped_v = grouped_v(tmp);

            risingVelocity_mean = cellfun(@mean,grouped_v);
            risingVelocity_std  = cellfun(@std, grouped_v);
            risingVelocity_max  = cellfun(@max, grouped_v);
            risingVelocity_min  = cellfun(@min, grouped_v);

            heightshift = min(height_global);

            % get the averaged bottomDetachmentHeight of selected records
            bottomDetachmentHeight = obj.get_averaged_bottomDetachmentHeight(arraylabel);

            % get the averaged collisionHeight of selected records
            collisionHeight = obj.get_averaged_collisionHeight(arraylabel);

            % get the number of record which are used for averaging of each height
            numberOfRecords = cellfun(@(d) size(d,1), grouped_v);

            [px2mm, fps] = obj.get_mean_SIunitsFactors(arraylabel);

            out = struct( ...
                'X' ,                            (height_global-heightshift) * px2mm, ...
                'Y_mean',                                risingVelocity_mean * px2mm * fps, ...
                'Y_std',                                 risingVelocity_std  * px2mm * fps, ...
                'Y_max',                                 risingVelocity_max  * px2mm * fps, ...
                'Y_min',                                 risingVelocity_min  * px2mm * fps, ...
                'collisionEvent',               (collisionHeight-heightshift)* px2mm, ...
                'bottomDetachmentEvent', (bottomDetachmentHeight-heightshift)* px2mm, ...
                'numberOfRecords',                            numberOfRecords       , ...
                'arraylabel',                                      arraylabel         ...
            );
        end

    end
    methods ( Access = private )
        function out = is_cleanVersions(obj,arraylabel)
            %IS_CLEANVERSION method to check if selected sequences created by commit clean gitVersion

            out = any(cellfun(@(d) isempty(d.commitCLEAN), obj.DATA(arraylabel)));
        end
        function out = is_same_gitVersion(obj,arraylabel)
            %IS_SAME_GITVERSION method to check if selected sequences have the same commitID

            out = size(unique(cellfun(@(d) d.commitID, obj.DATA(arraylabel),'UniformOutput', false)),1) == 1;
        end
        function plot_warnings(obj, arraylabel, varargin)
            %PLOT_WARNINGS method to show warnings when plotting records with different properties

            if iscell(varargin{:})
                varargin = varargin{:}{:};
            end

            if ~ismember('warningOFF', varargin)
                para = obj.label_same_parameterSet(arraylabel);
                if  size(unique(para(para~=0)),1)~=1, warning('plot_warnings:same_parameterSet', 'Selected Records dont have the same parameterSet.'); end
                temp = obj.label_same_temperature (arraylabel);
                if  size(unique(temp(temp~=0)),1)~=1, warning('plot_warnings:same_temperature',  'Selected Records dont have the same temperature.'); end

                if ~ismember('checkVersionsOFF', varargin)
                    if ~obj.is_same_gitVersion(arraylabel), warning('plot_warnings:same_gitVersion', 'Selected Records dont have the same gitVersion.'); end
                    if ~obj.is_cleanVersions(arraylabel)  , warning('plot_warnings:clean_Versions', 'Selected Records evaluated on a unclean gitVersion.'); end
                end
            end
        end
    end

end

