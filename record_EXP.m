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
classdef record_EXP < handle
    %RECORD_EXP class to evaluate single data of experiments

    properties
        reduce_usageOfRam = true; % reduce the usage of RAM by deleting big variables
        path; % directory in which the experimental data with parameter.csv/mat files can be found

        parameters; % struct array of parameters of the test facility
        parameters_csv_file = 'parameters.csv'; % CSV-Name of the parameter's file
        parameters_mat_file = 'parameters.mat'; % MAT-Name of the parameter's file

        imageFileTypes = {'bmp';'jpg'};  % all supported image file types

        STDfactor = 2; % factor of standard deviation, needed for the estimation of the threshold, see get_threshold()
        min_EventDistance = 1; % [px] minimal distance between objects for event detection, see get_bottomDetachmentTime() and get_collisionTime()

        smoothing_range  = 10; % [frames] PostProcessing smoothing range for moving average, see smooth_dropsDynamic()
        cam_syst = { ...  % Definition of the camera systems' specific image noise ratios
                     ...  %ID   maxRelStdNoises [%]
                           1 ,       0.025;...
                           2 ,       0.011;...
                           3 ,       0.03  ...
        };

        frames = struct(... % image data array
                     'filepaths',       {''},       ... % single image paths, string Cell-Array
                     'imgOriginal',     {[]},       ... % original images , int8/int12 Cell-Array
                     'imgFiltered',     {[]},       ... % result of the filtered original images, double Cell-Array
                     'imgBinary',       {[]},       ... % matrix which masks pixels of the objects, result after applied threshold, boolean Cell-Array
                     'numberFrames',    0,          ... % number of image frames of the experiment, double
                     'size',            zeros(2,1), ... % dimension scale of the image frame matrix, double
                     'px2mm',           0,          ... % mm/px-ratio, scale to SI-Units, [mm/px], double
                     'maxRelStdNoises', 0.025       ... % maximum of the relative standard deviation of filtered pixel values, double
        );
        cannulas = struct( ...  % cannulas data array
                     'imgKum' ,     [],           ... % cumulated Image of all Frames , double (image size)
                     'imgBinary',   [],           ... % binarised Images of imgKum; masks cannuala pixels, boolean (image size)
                     'tops_m' ,     zeros(2,1),   ... % coordinates of the midway cannula's top [px]; complex double
                     'tops' ,       zeros(2),     ... % coordinates of the edges of the cannual's top [px]; complex double
                     'roots',       zeros(2),     ... % coordinates of the edges of the cannual's root [px]; complex double
                     'type',        [1; 2],       ... % type of the cannuala [1 -> top, 2 -> bottom]
                     'diameters',   zeros(2,1),   ... % diameter of the cannulas [px], double
                     'length',      zeros(2,1),   ... % length of the cannualas [px], double
                     'diameter_mm', 0.8           ... % diameter of cannulas according to producer [mm], double
        );
        drops = struct( ...     % drops data array for each frame
                     'centers',    {[];[]}, ... % coordinates of drop's center [px], complex double
                     'edgePoints', {{};{}}, ... % coordinates of drop's edge points [px], complex double
                     'mean_radii', {[];[]}, ... % mean drop radius [px], double
                     'min_radii',  {[];[]}, ... % minimal drop radius [px], double
                     'max_radii',  {[];[]}, ... % maximal drop radius [px], double
                     'velocities', {[];[]}, ... % velocity of the drop centres [px/frames], double
                     'validates',  {[];[]}  ... % fulfil drop's criteria [boolean], boolean
        );
        Events  = struct( ... % event data array
                     'bottomDetachmentTime',zeros(1,2),  ... % Time when the bottom drop detachs from the lower cannula [framenumber], double
                     'collisionTime',       zeros(1,2),  ... % Time when both drop's surfaces touch each other [framenumber], double
                     'coalescenceTime',     zeros(1,2),  ... % Time when both drop's coalesce [framenumber], double
                     'repulsionTime',       zeros(1,2)   ... % Time when the drop contact breaks [framenumber], double
        );
        manualValues = struct ( ... % manual evalutated event and preprocessing data
                     'bottomDetachmentTime', [], ... % Time when the bottom drop detaches from the lower cannula [framenumber], double
                     'collisionTime',        [], ... % Time when both drops surfaces touch each other [framenumber], double
                     'coalescenceTime',      [], ... % Time when both drops coalesce [framenumber], double
                     'repulsionTime',        [], ... % Time when the drop contact breaks [framenumber], double
                     'crop_coordinates',     [], ... % diagonal rectangle edges to crop all frames in preprocessing, complex numbers
                     'man_threshold',        []  ... % set threshold by hand for binary frames after filtering, double
                 );
        error = [];            % stores thrown error
        log = {};              % logs methods which are done
        display_logs = false;  % display msg which are added
    end

    methods ( Static )
        function out = create_ParameterMAT(path,framerate)
            %CREATE_PARAMETERMAT create parameterMAT-file from directory name
            %
            % use this method if the parameter.csv-file is not available for an experiment, but be careful to give correct framerates
            % get drop size from the name of the directory

            if ~exist(path,'dir'), error('record_EXP:FalsePath','path does not exist!'); end
            if ~isa(framerate,'double') || isnan(framerate) || isempty(framerate) , error('record_EXP:FramerateNotDouble','Datatype of ''framerate'' is not double(1,1)!'),   end;
            if ispc && path(end) ~='\' % for windows
                    path = [path,'\'];
            elseif ~ispc && path(end) ~='/' % for Linux & Mac
                    path = [path,'/'];
            end
            if ~exist([path,'parameters.csv'],'file')
                para = regexp(path(1:end-1),'(?<d_disp_top_target>\d{1,3}.\d{1,3})_(?<d_disp_bot_target>(\d{1,3}.\d{1,3})|(\d{1,3}))_(?<drop_edge_distance>(\d{1,3}.\d{1,3})|(\d{1,3}))mm_(?<record_num>\d{1,3})$','names');
                if isempty(para)
                    para = regexp(path(1:end-1),'(?<d_disp_top_target>\d{1,3}.\d{1,3})_(?<d_disp_bot_target>(\d{1,3}.\d{1,3})|(\d{1,3}))_(?<drop_edge_distance>(\d{1,3}.\d{1,3})|(\d{1,3}))_(?<record_num>\d{1,3})$','names');
                    if isempty(para)
                        para = regexp(path(1:end-1),'(?<d_disp_top_target>\d{1,3}.\d{1,3})_(?<d_disp_bot_target>(\d{1,3}.\d{1,3})|(\d{1,3}))_(?<drop_edge_distance>(\d{1,3}.\d{1,3})|(\d{1,3}))mm_(?<temperature>\d{1,3})_(?<record_num>\d{1,3})$','names');
                        if isempty(para)
                            para = regexp(path(1:end-1),'(?<d_disp_top_target>\d{1,3}.\d{1,3})_(?<d_disp_bot_target>(\d{1,3}.\d{1,3})|(\d{1,3}))_(?<drop_edge_distance>(\d{1,3}.\d{1,3})|(\d{1,3}))mm(?<record_num>\d{1,3})$','names');
                            if isempty(para)
                                para = regexp(path(1:end-1),'(?<d_disp_top_target>\d{1,3}.\d{1,3})_(?<d_disp_bot_target>(\d{1,3}.\d{1,3})|(\d{1,3}))_(?<drop_edge_distance>(\d{1,3}.\d{1,3})|(\d{1,3}))mm_(?<temperature>\d{1,3})_(?<record_num>\d{1,3})$','names');
                                if isempty(para)
                                    para = regexp(path(1:end-1),'(?<d_disp_top_target>\d{1,3}.\d{1,3})_(?<d_disp_bot_target>(\d{1,3}.\d{1,3})|(\d{1,3}))_(?<drop_edge_distance>(\d{1,3}.\d{1,3})|(\d{1,3}))mm_(?<temperature>\d{1,3})°C_(?<record_num>\d{1,3})$','names');
                                    if isempty(para)
                                        para = regexp(path(1:end-1),'(?<d_disp_top_target>\d{1,3}.\d{1,3})_(?<d_disp_bot_target>(\d{1,3}.\d{1,3})|(\d{1,3}))_(?<drop_edge_distance>(\d{1,3}.\d{1,3})|(\d{1,3}))mm_(?<record_num>\d{1,3})$','names');
                                    end
                                end
                            end
                        end
                    end
                end

                parameters = structfun(@str2double, para, 'UniformOutput', 0); %#ok<*PROP>

                parameters.Framerate = framerate;

                % Set Units
                parameters.units.d_disp_top_target = '[mm]';
                parameters.units.d_disp_bot_target = '[mm]';
                parameters.units.drop_edge_distance = '[mm]';
                parameters.units.Framerate = '[frames/s]';
                if ispc && path(end) ~='\' % for windows
                    path = [path,'\'];
                elseif ~ispc && path(end) ~='/' % for Linux & Mac
                    path = [path,'/'];
                end
                save([path, 'parameters.mat'],'parameters');
                out = parameters;
            else
                out = false;
            end
        end
        function out = create_manualValueFile()
            %CREATE_MANUALVALUEFILE create the empty manualValues file

            global manualValue_file;
            if exist(manualValue_file,'file')
                warning('record_EXP:create_manualValueFile:ExistAlready', 'The manualValueFile exist already Path: %s',manualValue_file);
            else
                path = fileparts(manualValue_file);
                if ~exist(path,'dir'), mkdir(path); end

                fnames = {'current_directory';'bottomDetachmentTime';'collisionTime';'coalescenceTime'; ...
                          'repulsionTime';'hashID';'crop_coordinates';'man_threshold'...
                         };
                manualValues = cell2struct(cell(size(fnames,1),0),fnames);  %#ok<NASGU>

                save(manualValue_file, 'manualValues');
            end
            out = any(exist(manualValue_file,'file'));
        end
    end

    methods(  Access = public )

        function obj = record_EXP(path,CAM_ID)
            %RECORD_EXP contructor of record_EXP class
            %
            %  rec = record_EXP(path,CAM_ID)
            %
            % path - directory in which the experiment data with the parameter.csv/.mat are located
            % CAM  - id of the Camerasystem [1 for the old one, 2 for the new one]

            if ~exist(path,'dir'), error('record_EXP:FalsePath','path does not exist'); end
            config;
            obj.path  = correct_path_ending(path); % correct path ending if necessary

            %% ORGANISATION -- STEP ONE
            obj.load_parameters(); % load parameter from mat-file or generate mat-file form csv-file

            % set camera system specific data
            if exist('CAM_ID','var')
                obj.load_cameraSystem(CAM_ID);
            end
            obj.load_frames();
            obj.set_hashID();               % Generate hashID by hashing the first images-file of the experiment
            obj.load_manualValues();   % Set manual evaluate data if available
            obj.save_parametersMAT();       % save parameterData to a MatLab-File
        end

        function reload_parameters(obj)
            %RELOAD_PARAMETERS reload parameter values

            if ~exist([obj.path,obj.parameters_csv_file],'file')
                error('record_EXP:paraCSVNotFound','paramters.csv file does not exist!\n PFAD: %s',obj.path);
            else
                obj.load_parametersCSV();
                obj.save_parametersMAT();
            end
        end
        function set_error(obj,e)
            %SET_ERROR set the error

            obj.error = e;
        end
        function run_all(obj)
            %RUN_ALL execute all evaluations methods in the right order

            %% IMAGE PROCESSING -- STEP TWO
            obj.filter_frames();      % Filter all images
            obj.detect_imageCorruption();
            obj.binary_frames();      % generate binary images via threshold

            %% GET DATA
            obj.rotate2top_frames();
            obj.detect_cannuals();          % detect cannulas (bottom and top)

            if obj.reduce_usageOfRam    % Reduction unused variable Memory
                obj.frames.imgFiltered = [];
                obj.cannulas.imgKum    = [];
                %obj.cannulas.imgBinary = [];
            end

            obj.set_px2mm();

            obj.reset_DropProperties();
            obj.get_dropCenters(); % detect drops with prediction

            obj.smooth_dropsDynamic();  % smooth drop trajectories; necessary because of the discrete input data
            %plot(real(obj.drops(2).centers));
            %drawnow;
            obj.get_bottomDetachmentTime(); % limits bottomDetachmentTime, see definition
            obj.get_collisionTime();    % limits collsionTime, see definition

            %in progress
            % obj.get_dropCentersAfterCollision();
            %coming soon
            % obj.get_coalescenceTime();
            % obj.get_repulsionTime();
            % obj.get_contactTime();
            % obj.get_drainageTime();
            obj.save_DATA();    % save outcome data
        end

        %% FRAMES
        function out = detect_imageCorruption(obj)
            %DETECT_IMAGECORRUPTION detect image corruption caused by maladjusted ratio of fps and image size

            out = true;
            imgF = obj.frames.imgFiltered{1};

            for i = 2:obj.frames.numberFrames
                imgF = imgF + obj.frames.imgFiltered{i};
            end
            imgF = imgF / obj.frames.numberFrames ;

            imgFx = any(sum(imgF,2) /obj.frames.size(2)/mean(imgF(:)) < 0.25);
            imgFy = any(sum(imgF,1) /obj.frames.size(1)/mean(imgF(:)) < 0.25);
            if imgFx || imgFy
                error('record_EXP:detect_imageCorruption:founded','Detect image corruption!');
            end
        end
        function out = load_frames(obj,i)
            %LOAD_FRAMES load the images data of the experiment to obj.frames.imgOriginal as cell array

            if isempty(obj.frames.filepaths)
                % get all image-files from experiment directory
                obj.frames.filepaths = sort(getAllFiles(obj.path, ['_\d*',strjoin(strcat('(.',obj.imageFileTypes,'$)')','|')] ));

                % allocate memory
                obj.frames.numberFrames = size(obj.frames.filepaths,1);
                obj.frames.imgOriginal = cell(obj.frames.numberFrames,1);

                % Read image files via MatLab Image Processing Toolbox
                for i = 1:obj.frames.numberFrames
                    obj.frames.imgOriginal{i} = imread(obj.frames.filepaths{i});
                    if 3 == size(obj.frames.imgOriginal{i},3)
                        obj.frames.imgOriginal{i} = rgb2gray(obj.frames.imgOriginal{i});
                    else
                        obj.frames.imgOriginal{i} = obj.frames.imgOriginal{i};
                    end
                end
                obj.frames.size = size(obj.frames.imgOriginal{1});
            else exist('i','var')
                out = imread(obj.frames.filepaths{i});
                if 3 == size(out,3)
                    out= rgb2gray(out);
                end
            end

        end
        function out = filter_frames(obj,i)
            %FILTER_FRAMES filter the images to reduce large range gradients along image length

            n = 31;
            gauss = fspecial('gaussian',[n,n],n/3); % two dimensional Gaussian's Bell for a radial weight
            divEdge = conv2(ones(obj.frames.size),gauss,'same');

            if exist('i','var')
                I   = double(obj.load_frames(i));
                out = I./(conv2(I,gauss,'same')./divEdge) ;
            else
                obj.frames.imgFiltered  = cellfun(@(I)  mat2gray(double(I)./(conv2(double(I),gauss,'same')./divEdge)) , obj.frames.imgOriginal,'UniformOutput',false);
            end
%             if obj.reduce_usageOfRam
%                 obj.frames.imgOriginal = {};
%             end

        end
        function binary_frames(obj)
            %BINARY_FRAMES generate image masks to select object pixels of every frame via threshold estimated for each frame

            obj.frames.imgBinary  =  cellfun(@(I) bwareaopen( ... % reduce noised object pixel through
                                                              ... % the criteria of the number of their
                                                              ... % grouped pixels
                                                  I < obj.get_threshold(I) ... % create binary mask for object pixels
                                                  ,6,8 ) ...
                                     , obj.frames.imgFiltered,'UniformOutput',false);

            % choise without noise reduction
            % obj.frames.imgBinary  =  cellfun(@(I) I < obj.get_threshold(I) , obj.frames.imgFiltered,'UniformOutput',false);
        end

        %% OBJECTS
        function detect_cannuals(obj)
            %DETECT_CANNUALS find cannulas through selecting object pixel via specific criterion for cannulas

            obj.set_cannulasImgKum();    % cumulate all filtered Images
            obj.set_cannulasImgBinary(); % create a mask for cannula pixel

            obj.get_cannulasAreas();     % select the cannula objekt from mask
            obj.improve_cannulasShape(); % erase defects on the cannula shape

            obj.set_cannulasExtrema();   % find the four extrema with the far distant of cannula's centroid
            obj.set_cannulasDiameter();  % estimate the diameter of the cannula
            obj.set_cannulasLength();    % estimate the length of the cannula
            obj.set_cannulasTopsM();     % estimate the the middle of the top cannula
        end
        function rotate2top_frames(obj)
            %ROTATE2TOP_FRAMES rotate the images, so that the interested drop is rising.

            n = (1:obj.frames.numberFrames/2)';
            X = 1:obj.frames.size(1); % get x-coordinate image limits
            Y = 1:obj.frames.size(2); % get y-coordinate image limits

            % estimate the centroid of object pixcels
            X_sum_int = cellfun(@(img) sum(img,2)'  ,obj.frames.imgBinary(n),'UniformOutput',false);
            X_centroid = cell2mat(cellfun(@(y_sum) y_sum*X'/sum(y_sum), X_sum_int,'UniformOutput',false));

            % eliminate  frames which are disturbing
            dX = abs(diff(X_centroid));
            int_X = get_interval( dX < 20);
            int_X = int_X(end,:);

            Y_sum_int = cellfun(@(img) sum(img,1) ,obj.frames.imgBinary(n),'UniformOutput',false);
            Y_centroid = cell2mat(cellfun(@(x_sum) x_sum*Y'/sum(x_sum), Y_sum_int,'UniformOutput',false));

            % eliminate frames which are disturbing
            dY = abs(diff(Y_centroid));
            int_Y = get_interval( dY < 20);
            int_Y = int_Y(end,:);

            % linear regression to estimate the global direction of centroid's dynamic
            ax = polyfit((int_X(1):int_X(2))',X_centroid(int_X(1):int_X(2)),1);
            ay = polyfit((int_Y(1):int_Y(2))',Y_centroid(int_Y(1):int_Y(2)),1);

            % reduce to gradient of the first order
            ax = ax(1);
            ay = ay(1);

            % show error if drop dynamic is not detected
            if hypot(ax,ay) < 0.5/obj.parameters.Framerate, error('record_EXP:noDynamicDetected','No Dynamik detected in record!'); end

            % get number of the 90°-rotation which is needed for the correct image orientation
            rot = find((0:4)*pi/2 < (angle((-(ax) + ay*1i)*exp(-3*pi/4*1i))+pi),1,'last')-1;

            % validate the nummber of the rotation
            % check for exiting object sticked on frame edges
            x_ = sum(cell2mat(cellfun(@(x) any(x(1)) & any(x(end)), X_sum_int ,'UniformOutput',false))) >= size(n,1);
            y_ = sum(cell2mat(cellfun(@(y) any(y(1)) & any(y(end)), Y_sum_int ,'UniformOutput',false))) >= size(n,1);
            if xor(x_, y_)
                if (mod(rot,2) == 0 && y_ ) || (mod(rot,2) == 1 && x_)
                    rot = mod(rot+1,4);
                end
            elseif 1 < ceil(log10(abs(ax))) - ceil(log10(abs(ay)))
                error('record_EXP:rotate2top_frames:noRotationDecision','Can''t decide how to rotate the frames.');
            end

            if  rot ~= 0
                obj.frames.imgOriginal = cellfun(@(I) rot90(I,rot),obj.frames.imgOriginal,'UniformOutput', false);
                obj.frames.imgFiltered = cellfun(@(I) rot90(I,rot),obj.frames.imgFiltered,'UniformOutput', false);
                obj.frames.imgBinary   = cellfun(@(I) rot90(I,rot),obj.frames.imgBinary,'UniformOutput', false);

                if mod(rot,2) == 1
                    obj.frames.size = fliplr(obj.frames.size);
                end
            end
        end
        function reset_DropProperties(obj)
            %RESET_DROPPROPERTIES reset and initialize values for drop detection

            % reset validates values
            obj.drops(1).validates = false(obj.frames.numberFrames,1); % top drop
            obj.drops(2).validates = false(obj.frames.numberFrames,1); % bottom drop

            % reset drop velocities
            obj.drops(1).velocities = zeros(obj.frames.numberFrames,1); % top drop
            obj.drops(2).velocities = zeros(obj.frames.numberFrames,1); % bottom drop

            % reset drop edgePoints
            obj.drops(1).edgePoints = cell(obj.frames.numberFrames,1); % top drop
            obj.drops(2).edgePoints = cell(obj.frames.numberFrames,1); % bottom drop

            % reset minDistEdgePoints
            obj.drops(1).minDistEdgePoints = NaN(obj.frames.numberFrames,1); % top drop
            obj.drops(2).minDistEdgePoints = NaN(obj.frames.numberFrames,1); % bottom drop


            % get estimated radii
            radii_top = obj.parameters.d_disp_top_target/obj.frames.px2mm/2;
            radii_bot = obj.parameters.d_disp_bot_target/obj.frames.px2mm/2;

            % set estimated radii as default
            obj.drops(1).min_radii = radii_top*ones(obj.frames.numberFrames,1); % top drop
            obj.drops(2).min_radii = radii_bot*ones(obj.frames.numberFrames,1); % bottom drop

            obj.drops(1).max_radii = radii_top*ones(obj.frames.numberFrames,1); % top drop
            obj.drops(2).max_radii = radii_bot*ones(obj.frames.numberFrames,1); % bottom drop

            obj.drops(1).mean_radii = radii_top*ones(obj.frames.numberFrames,1); % top drop
            obj.drops(2).mean_radii = radii_bot*ones(obj.frames.numberFrames,1); % bottom drop

            % set the estimated drop centers relative to their associated cannula
            obj.drops(1).centers = obj.cannulas.tops_m(1) + obj.drops(1).mean_radii(1)*ones(obj.frames.numberFrames,1); % top drop
            obj.drops(2).centers = obj.cannulas.tops_m(2) - obj.drops(2).mean_radii(1)*ones(obj.frames.numberFrames,1); % bottom drop
        end
        function get_dropCentersAfterCollision(obj)
            %GET_DROPCENTERSAFTERCOLLISION seperate drop pixel into their associated drop types after
            %  collision via moving dividing line constructed through the collsion point

            begin_i = obj.Events.collisionTime(1);%-1+find(obj.drops(2).validates(obj.Events.collisionTime(2):end),1,'first');
            end_i =  begin_i+find(obj.drops(2).validates((begin_i+1):end),1,'last') ;
            [Y, X] = meshgrid(1:obj.frames.size(2),1:obj.frames.size(1));
            XY = X + 1i*Y;
            weight_score = 2/5;  % score_centroid/score_sumed
            for i = begin_i:end_i

                if isnan(obj.drops(1).minDistEdgePoints(i-1))  || isnan(obj.drops(2).minDistEdgePoints(i-1)) || ...
                   isnan(obj.drops(2).velocities(i-1))
                   break;
                end
                % estimate of central edge drop point O
                O = (obj.drops(1).minDistEdgePoints(i-1)+obj.drops(2).minDistEdgePoints(i-1))/2+obj.drops(2).velocities(i-1) ;

                % get drop pixel
                BW = obj.cleared_cannulasImgBinary(i);
                [label,num] = bwlabel(BW,8);
                tmp = label~=0 & abs( XY  - O ) < max(obj.drops(1).mean_radii(i-1),obj.drops(2).mean_radii(i-1));
                tmp1 =  unique(label(tmp));
                label(ismember(label,tmp1)) = tmp1(1);
                num = num - size(tmp1,1) +1;

                % eliminate not interesting objects
                if num ~= 1
                    num = (1:num)';

                    areas     = arrayfun(@(BWl) sum(label(:) == BWl),num);
                    centroids = abs(O - arrayfun(@(BWl) (X(:,1)'*sum(label== BWl,2)+Y(1,:)*sum(label== BWl,1)'*1i)/sum(label(:) == BWl),num));

                    score_centroids = 1-(centroids-min(centroids))/(max(centroids)-min(centroids));
                    score_areas     =   (areas    -min(areas))    /(max(areas)    -min(areas));
                    score = score_centroids*weight_score + score_areas*(1-weight_score);
                    tmp  = score == max(score);
                else
                    tmp = 1;
                end

                BWdrops = label == num(tmp);
                Drops_Points = XY(BWdrops(:)); % coordinates of the drop points

                BWdrops_tmp = BWdrops ;

                % create two global maximum as limits
                x1 = round(real(obj.drops(1).centers(i-1)));
                x2 = round(real(obj.drops(2).centers(i-1)));
                BWdrops_tmp([x1,x2],:) = true;

                % fill drops with minor drop edge circle
                BWdrops_tmp(abs(XY(:)-obj.drops(1).centers(i-1)) < obj.drops(1).min_radii(i-1)-1 ) = true;
                BWdrops_tmp(abs(XY(:)-obj.drops(2).centers(i-1)) < obj.drops(2).min_radii(i-1)-1 ) = true;
                tmp2 = sum(imfill(BWdrops_tmp(x1:x2,:),'holes'),2);

                % find the minimal pixel number in x-coordinated
                [~, id] = min(tmp2);

                % estimate the base of the dividing line
                O = (x1 + id - 1) + 1i*imag(O)  ;

                % construct dividing line
                ang = (pi/2)-angle(obj.drops(1).centers(i-1) - obj.drops(2).centers(i-1)); % angle of the dividing line
                tmp = angle((O-Drops_Points)*exp(1i*ang));  % coordinate translation to the base of the dividing line

                % set types of the pixels through their position relative to the constructed line
                type = ones(size(tmp))*2;
                type(tmp<0) = ones(sum(tmp<0),1)*1;

                obj.drops(1).edgePoints{i} = Drops_Points(type == 1);
                obj.drops(2).edgePoints{i} = Drops_Points(type == 2);

                if isempty(obj.drops(1).edgePoints{i}) || isempty(obj.drops(2).edgePoints{i}), break; end

                for type_ = 1:2 % for each drop type
                    % reduce points to edge points only
                    koords_on_dropCenter = obj.drops(type_).edgePoints{i}-obj.drops(type_).centers(i-1);
                    pixel_props =  [ obj.drops(type_).edgePoints{i},    ... % index 1
                                     koords_on_dropCenter,              ... % index 2
                                     abs(koords_on_dropCenter),         ... % index 3
                                     angle(koords_on_dropCenter)        ... % index 4
                                   ];
                    % reduce pixels via twice the standard deviation of the distance between drop pixel and estimated center
                    tmp = ~( mean(pixel_props(:,3)) + 2*std(pixel_props(:,3)) < pixel_props(:,3) );
                    pixel_props = pixel_props(tmp,:);

                    % sort pixel by angle
                    pixel_props  = sortrows([pixel_props, angle(pixel_props(:,2))],4);
                    pixel_props = [flipud(pixel_props(real(pixel_props(:,4)) < 0,:)); pixel_props( 0 <= real(pixel_props(:,4)),:)]; % negativ values are not considered by the sortrows function

                    % periodic continuation of the pixels' list
                    pixel_props = [pixel_props; pixel_props(1,:)];
                    pixel_props(end,4) = pixel_props(end,4)+2*pi;

                    % local major search with defined angle range
                    [~, index] = findpeaks_distance(pixel_props(:,4), pixel_props(:,3), acos(1-2/obj.drops(type_).mean_radii(i-1)^2));
                    pixel_props = pixel_props(index,:); % reduce to point with maximal radius inside the angle range

                    obj.drops(type_).edgePoints{i} = pixel_props(:,1);

                    %% adjustment of the estimated drop center by the detected drop edge through weights centroid determination
                      % get angular difference for each edge points between their two adjacent edge points
                      diffAngle2neighbours = diff([pixel_props(end,4)-2*pi; pixel_props(:,4)])  ...
                                            +diff([pixel_props(:,4); pixel_props(1,4)+2*pi]);

                      % adjust the estimated drop center
                      dS = diffAngle2neighbours'*pixel_props(:,2)/(4*pi); % Korrektur des Tropfenmittelpunktes
                      obj.drops(type_).centers(i) = obj.drops(type_).centers(i-1) + dS;

                    if i ~= 1
                        obj.drops(type_).velocities(i) = obj.drops(type_).centers(i)-obj.drops(type_).centers(i-1);
                    end

                    % adjust radius of the drop
                    radiiEdgePoints = abs(pixel_props(:,3) - dS);
                    obj.drops(type_).mean_radii(i) = mean(radiiEdgePoints);
                    obj.drops(type_).min_radii(i)  = min(radiiEdgePoints);
                    obj.drops(type_).max_radii(i)  = max(radiiEdgePoints);
                    obj.drops(type_).validates(i)  = true;
                end

                % select drop edge points inside lower quadrantal for the top drop
                lowerRAD_drop1 = obj.drops(1).edgePoints{i}(angle((obj.drops(1).edgePoints{i}-obj.drops(1).centers(i))*exp(+1i*pi/2)) > 0);
                % select drop edge points inside upper quadrantal for the bottom drop
                upperRAD_drop2 = obj.drops(2).edgePoints{i}(angle((obj.drops(2).edgePoints{i}-obj.drops(2).centers(i))*exp(-1i*pi/2)) > 0);

                if isempty(lowerRAD_drop1) || isempty(upperRAD_drop2), continue;  end

                % estimate distance between each selected drop edge point with different type
                kronLowerRAD_drop1 = kron(ones(size(upperRAD_drop2)),lowerRAD_drop1.');
                kronUpperRAD_drop2 = kron(ones(size(lowerRAD_drop1))',upperRAD_drop2 );
                distantDropEdgePoints = abs(kronLowerRAD_drop1-kronUpperRAD_drop2);

                % select drop edge points with minimal distance to each other
                minDIST = distantDropEdgePoints == min(distantDropEdgePoints(:));
                kronLowerRAD_drop1 = kronLowerRAD_drop1(minDIST);
                kronUpperRAD_drop2 = kronUpperRAD_drop2(minDIST);

                % select drop edge points with minimal angular difference to perpendicular
                 P = mean([kronLowerRAD_drop1,kronUpperRAD_drop2],2);
                 [~,id] = min(abs(  angle(P-obj.drops(1).centers(i)) ...
                                     - angle(obj.drops(2).centers(i)- obj.drops(1).centers(i)) ...
                                    ) ...
                                );

                obj.drops(1).minDistEdgePoints(i) = kronLowerRAD_drop1(id(1),:);
                obj.drops(2).minDistEdgePoints(i) = kronUpperRAD_drop2(id(1),:);

%                 vec = view_EXP(obj);
%                 vec.fh=1;
%                     figure(vec.fh);
%                     vec.show_imgBinary(i);
%                     hold on;
%                     n = exp(1i*(-ang));
%                     plot(imag(O),real(O),'+b');
%                     X_ = [-50;50]*n+O;
%                     plot(imag(X_),real(X_),'-r')
%                     drawnow;
%                     hold off;
            end
        end

        %% EVENTS
        function get_bottomDetachmentTime(obj)
            %GET_BOTTOMDETACHMENTTIME restricts the time interval of the bottomDetachmentTime

            %% FIRST RESTRICTION LEVEL (manual restriction)
            int = 1:ceil(obj.frames.numberFrames/3); % gives a loose manual interval

            %% SECOND RESTRICTION LEVEL (drop approximately descripted by a circular ring)
            distance_bottomDropCannual = abs(obj.cannulas.tops_m(2)- obj.drops(2).centers(int)); % distance between drop center and top of the bottom cannula
            EdgeDistance_max =  distance_bottomDropCannual - obj.drops(2).max_radii(int);
            EdgeDistance_min =  distance_bottomDropCannual - obj.drops(2).min_radii(int);

            max_ = get_interval(obj.min_EventDistance < EdgeDistance_max)  ;
            max_ = max_(find( 10 < diff(max_,1,2) ,1,'first'),1);

            min_ = get_interval( obj.min_EventDistance < EdgeDistance_min);
            min_ = min_(find( 10 < diff(min_,1,2) ,1,'first'),1);

            % SET the limits of the restricted intervall
            if any(min_)
                obj.Events.bottomDetachmentTime(1) = min_;
            else
                obj.Events.bottomDetachmentTime(1) = int(1); % fallback to the first restriction level
            end
            if any(max_)
                obj.Events.bottomDetachmentTime(2) = max_;
            else
                obj.Events.bottomDetachmentTime(2) = int(end); % fallback to the first restriction level
            end

            %% THIRD RESTRICTION LEVEL (check distance between each drop edge pixels and pixels of top's bottom cannual)
            %
            int = obj.Events.bottomDetachmentTime(1):obj.Events.bottomDetachmentTime(2);
            if size(int,2)~=1
                minDIST = NaN(size(int));
                for i = int
                    lowerRAD_drop2 = obj.drops(2).edgePoints{i}(angle((obj.drops(2).edgePoints{i}-obj.drops(2).centers(i))*exp(+1i*pi/2)) > 0);
                    if isempty(lowerRAD_drop2), continue; end
                    dS = abs(diff(obj.cannulas.tops(2,:)));
                    upperRAD_cannula2 = obj.cannulas.tops(2,1) + (0:ceil(dS))/ceil(dS).*diff(obj.cannulas.tops(2,:));

                    kronLowerRAD_drop2    = kron(ones(size(upperRAD_cannula2)),lowerRAD_drop2);
                    kronUpperRAD_cannula2 = kron(ones(size(lowerRAD_drop2)),upperRAD_cannula2 );

                    distantDropEdgePointCannulaTop = abs(kronUpperRAD_cannula2-kronLowerRAD_drop2);

                    leftLowerDropEdgePoints  = imag(lowerRAD_drop2)-imag(obj.drops(2).centers(i)) < 0;
                    rightLowerDropEdgePoints = imag(lowerRAD_drop2)-imag(obj.drops(2).centers(i)) > 0;

                    if ~any(leftLowerDropEdgePoints(:)) || ~any(rightLowerDropEdgePoints(:)), continue; end

                    minDistLeft = distantDropEdgePointCannulaTop == min(min(distantDropEdgePointCannulaTop(leftLowerDropEdgePoints,:)));
                    minDistLeft(~leftLowerDropEdgePoints,:) = false;
                    id_minDistLeft = find(minDistLeft(:),1,'first');
                    minDistLeftPoint = kronLowerRAD_drop2(id_minDistLeft)-obj.drops(2).centers(i);

                    minDistRight = distantDropEdgePointCannulaTop == min(min(distantDropEdgePointCannulaTop(rightLowerDropEdgePoints,:)));
                    minDistRight(~rightLowerDropEdgePoints,:) = false;
                    id_minDistRight = find(minDistRight,1,'first');
                    minDistRightPoint = kronLowerRAD_drop2(id_minDistRight)-obj.drops(2).centers(i);

                    angle_LR =  angle(minDistLeftPoint)/angle(minDistRightPoint);
                    R_0 = ( abs(minDistLeftPoint) - abs(minDistLeftPoint)*angle_LR ) / ( 1-angle_LR );

                    min_DropEdgePoint = R_0 + obj.drops(2).centers(i);

                    minDIST(i) = min(abs(upperRAD_cannula2 - min_DropEdgePoint));
                end

                bDTime = find(minDIST<1.5,1,'last');
                if ~isempty(bDTime)
                     obj.Events.bottomDetachmentTime = ones(1,2)*bDTime;
                end
            end

            %% limit time interval
            obj.Events.bottomDetachmentTime(obj.Events.bottomDetachmentTime < 0 ) = 0;
            obj.Events.bottomDetachmentTime(obj.Events.bottomDetachmentTime > obj.frames.numberFrames) = obj.frames.numberFrames;

        end
        function get_collisionTime(obj)
            %GET_COLLISIONTIME restricts the interval of the collisionTime

            %% FIRST RESTRICTION LEVEL (manual restriction)
            %int = ceil(obj.frames.numberFrames/3):obj.frames.numberFrames ; % gives a loose manual interval
            int = (obj.Events.bottomDetachmentTime(2)+1):obj.frames.numberFrames ; % gives a loose manual interval

            %% SECOND RESTRICTION LEVEL (drop approximately descripted by a circular ring)
            distance_centers = abs(obj.drops(1).centers(int) - obj.drops(2).centers(int)); % distance between drops in a single image frame
            EdgeDistance_max = distance_centers - abs(obj.drops(1).max_radii(int) + obj.drops(2).max_radii(int));
            EdgeDistance_min = distance_centers - abs(obj.drops(1).min_radii(int) + obj.drops(2).min_radii(int));

            min_ = obj.Events.bottomDetachmentTime(2)+find( EdgeDistance_max < obj.min_EventDistance ,1,'first');
            max_ = obj.Events.bottomDetachmentTime(2)+find( EdgeDistance_min < obj.min_EventDistance ,1,'first');

            if any(min_)
                obj.Events.collisionTime(1) = min_;
            else
                obj.Events.collisionTime(1) = int(1); % fallback to the first restriction level
            end
            if any(max_)
                tmp = obj.Events.collisionTime(1)-1+find(~obj.drops(2).validates(obj.Events.collisionTime(1):end),1,'first');
                if any(tmp)
                    obj.Events.collisionTime(2) = min(max_,tmp);
                else
                    obj.Events.collisionTime(2) = max_;
                end
            else
                obj.Events.collisionTime(2) = int(end); % fallback to the first restriction level;
            end

            %% THIRD RESTRICTION LEVEL (check distance between each drop edge pixels in single image frame)
            int = obj.Events.collisionTime(1):obj.Events.collisionTime(2);
            if size(int,2)~=1
                for i = int
                    lowerRAD_drop1 = obj.drops(1).edgePoints{i}(angle((obj.drops(1).edgePoints{i}-obj.drops(1).centers(i))*exp(+1i*pi/2)) > 0);
                    upperRAD_drop2 = obj.drops(2).edgePoints{i}(angle((obj.drops(2).edgePoints{i}-obj.drops(2).centers(i))*exp(-1i*pi/2)) > 0);

                    if isempty(lowerRAD_drop1) || isempty(upperRAD_drop2), continue;  end

                    kronLowerRAD_drop1 = kron(ones(size(upperRAD_drop2)),lowerRAD_drop1.');
                    kronUpperRAD_drop2 = kron(ones(size(lowerRAD_drop1))',upperRAD_drop2 );
                    distantDropEdgePoints = abs(kronLowerRAD_drop1-kronUpperRAD_drop2);

                    % selection by minimal distance
                    minDIST = distantDropEdgePoints == min(distantDropEdgePoints(:));
                    kronLowerRAD_drop1 = kronLowerRAD_drop1(minDIST);
                    kronUpperRAD_drop2 = kronUpperRAD_drop2(minDIST);

                    % selection by minimal distance to centroid axis
                    P = mean([kronLowerRAD_drop1;kronUpperRAD_drop2],1);
                    [~,id] = min(abs(  angle(P-obj.drops(1).centers(i)) ...
                                     - angle(obj.drops(2).centers(i)- obj.drops(1).centers(i)) ...
                                    ) ...
                                );

                    obj.drops(1).minDistEdgePoints(i) = kronLowerRAD_drop1(id(1));
                    obj.drops(2).minDistEdgePoints(i) = kronUpperRAD_drop2(id(1));
                end
                dist  = abs(obj.drops(1).minDistEdgePoints(int) - obj.drops(2).minDistEdgePoints(int));
                index = int(1)-1+find(dist == 2 ,1,'last')+1;
                if ~isempty( index)
                    obj.Events.collisionTime(1:2) = index;
                end
            end

            %% limit time interval
            obj.Events.collisionTime(obj.Events.collisionTime < 0 ) = 0;
            obj.Events.collisionTime(obj.Events.collisionTime > obj.frames.numberFrames) = obj.frames.numberFrames;

        end
        function get_coalescenceTime(obj)
            %GET_COALESCENCETIME restricts the interval of the coalescenceTime

             % coming soon
            warning('record_EXP:get_coalescenceTime:TODO', 'Function ''obj.get_coalescenceTime()'' is in development.')
            obj.Events.coalescenceTime;
        end
        function get_repulsionTime(obj)
            %GET_REPULSIONTIME restricts the interval of the repulsionTime

             % coming soon
             warning('record_EXP:get_coalescenceTime:TODO', 'Function ''obj.get_coalescenceTime()'' is in development.')
            obj.Events.repulsionTime;
        end

        function get_contactTime(obj)
            %GET_CONTACTTIME restricts the interval of the contactTime

             % coming soon
            warning('record_EXP:get_contactTime:TODO', 'Function ''obj.get_contactTime()'' is in development.')
            obj.Events.contactTime;
        end
        function get_drainageTime(obj)
            %GET_DRAINAGETIME restricts the interval of the drainageTime

             % coming soon
            warning('record_EXP:get_drainageTime:TODO', 'Function ''obj.get_drainageTime()'' is in development.')
            obj.Events.drainageTime;
        end

        %% MISC
        function out = crop_frames(obj,coords)
            %CROP_FRAMES crop frame by two given coordinates

            out = false;
            if exist('coords','var')
                obj.manualValues.crop_coordinates = coords;
                obj.add_log('Coordinates for cropping images are overwritten.');
            end

            if ~isempty(obj.manualValues.crop_coordinates) && ~any(size(obj.manualValues.crop_coordinates) ~= [2,1])
                if  any(obj.manualValues.crop_coordinates < kron(1+1i,ones(2,1))) || any(kron(obj.frames.size(1)+obj.frames.size(2)*1i,ones(2,1)) <= obj.manualValues.crop_coordinates ), error('record_EXP:crope_frames:coordinatesNotInBorder','Given cooridnates are not defined within frames size.'); end;

                obj.frames.imgOriginal = cellfun(@(I)I(real(obj.manualValues.crop_coordinates(1)):real(obj.manualValues.crop_coordinates(2)), ...
                                                           imag(obj.manualValues.crop_coordinates(1)):imag(obj.manualValues.crop_coordinates(2))), ...
                                                           obj.frames.imgOriginal,'UniformOutput',false);
                obj.frames.size = size(obj.frames.imgOriginal{1});
                obj.add_log('Images are cropped by hand.');
                out = true;
            end
        end
        function save_DATA(obj)
            %SAVE_DATA saves archived data of the record

            % GIT STATUS
            global GIT_PATH;
            if exist(GIT_PATH, 'file')
                [~, commitID]    = system([GIT_PATH,' rev-parse HEAD']);
                [~, commitCLEAN] = system([GIT_PATH,' status -s']);
            else
                commitID = '';
                commitCLEAN = [];
            end

            data = struct( 'parameters' ,obj.parameters,      ...
                           'drops', obj.drops,                ...
                           'cannulas', obj.cannulas,          ...
                           'Events', obj.Events,              ...
                           'manualValues' , obj.manualValues, ...
                           'dateTime', datestr(now),          ...
                           'commitID', strtrim(commitID),     ...
                           'commitCLEAN',commitCLEAN,         ...
                           'log', {obj.log},                  ...
                           'frames', rmfield(obj.frames, {'imgOriginal','imgFiltered','imgBinary'}), ...
                           'error', obj.error                 ...
                       ); %#ok<NASGU>
            save([obj.path,'exp_data.mat'],'data');
        end
        function save_manualValuesFiles(obj)
            %SAVE_MANUALVALUESFILES save values set by hand

            global manualValue_file;
            load(manualValue_file);
            if ~exist('manualValues','var'), error('record_EXP:load_manualValues:noDATA','Loading manualValues from file ''%s'' failed.', manualValue_file); end

            % add fields if they are not existent
            fnames = {'current_directory';'bottomDetachmentTime';'collisionTime';'coalescenceTime'; ...
                      'repulsionTime';'hashID';'crop_coordinates';'man_threshold'...
                     };
            tmp = ismember(fnames, fieldnames(manualValues)); %#ok<NODEF>
            if any(~tmp)
                fnames = fnames(tmp);
                for i = 1:size(fnames,1)
                    manualValues(:).(fnames{i}) = {};
                end
            end

            tmp = ismember({manualValues(:).hashID},obj.parameters.hashID);
            if sum(tmp) == 1
                manualValues(tmp).crop_coordinates =  obj.manualValues.crop_coordinates;
                manualValues(tmp).man_threshold    =  obj.manualValues.man_threshold; %#ok<NASGU>
            elseif ~any(tmp)
                add = struct('current_directory',    obj.path, ...
                             'bottomDetachmentTime', [] , ...
                             'collisionTime',        [] , ...
                             'coalescenceTime',      [] , ...
                             'repulsionTime',        [] , ...
                             'hashID',               obj.parameters.hashID,  ...
                             'crop_coordinates',     obj.manualValues.crop_coordinates, ...
                             'man_threshold',        obj.manualValues.man_threshold     ...
                       );
                manualValues = [manualValues; add];  %#ok<NASGU>
            end
            save(manualValue_file,'manualValues');
            obj.add_log('manualValues saved.')
        end
        function disp_fullLog(obj)
            %DISP_FULLLOG diplay a list of messages created by method used with unusual adjustments for replicability

            disp(strjoin([{[obj.path,':']},obj.log],sprintf('\n')));
        end
    end
    methods (Access = private)

        %% FRAMES
        function out = get_threshold (obj, I)
            %GET_THRESHOLD estimate threshold for given image based on properties of the camera system which is used

            if isempty(obj.manualValues.man_threshold)
                % decrease average instensity by the static standard deviation of the camera system
                out = mean(I(:))*(1-obj.STDfactor*obj.frames.maxRelStdNoises);
            else
                out = obj.manualValues.man_threshold;
                obj.add_log('Threshold value is selected by hand.');
            end
        end

        %% CANNUALS
        function set_cannulasImgKum(obj)
            %SET_CANNULASIMGKUM generate average image by arithmetic averaging pixel values with same image position

            obj.cannulas.imgKum = zeros(obj.frames.size);
            for n = 1:obj.frames.numberFrames
                obj.cannulas.imgKum = obj.cannulas.imgKum +  obj.frames.imgFiltered{n};
            end
            obj.cannulas.imgKum = obj.cannulas.imgKum/obj.frames.numberFrames; % normalize
        end
        function set_cannulasImgBinary(obj)
            %SET_CANNULASIMGBINARY generate image masks to select object pixels of cannula imgKum via threshold

            obj.cannulas.imgBinary = bwareaopen(obj.cannulas.imgKum < obj.get_threshold(obj.cannulas.imgKum),4,8);
        end
        function get_cannulasAreas(obj)
            %GET_CANNULASAREAS select and label cannulas pixel

            imgLabeled = bwlabel(obj.cannulas.imgBinary,4); % label linked object pixels

            % eliminate objects which have  no pixel on top or bottom frame edges
            index_topEdge = unique(imgLabeled(1,:));          % ( 0 is the background)
            index_topEdge = index_topEdge(index_topEdge ~=0); % erase index 0
            index_bottomEdge = unique(imgLabeled(end,:));
            index_bottomEdge = index_bottomEdge(index_bottomEdge ~=0); % erase index 0

            if isempty(index_topEdge) && isempty(index_bottomEdge), error('cannulas:find:canNOTonEdge','No cannulas found at image frame edge!'); end

            %% Calculate the area and the object centroids and select one of each edge
            %  with central centroid in y-coordinate and the biggest area as cannula

            Y = (1:obj.frames.size(2))' - obj.frames.size(2)/2; % shift coordinate system so that the image y-center is at 0
            weight_score = 2/5;  % score_centroid/score_sumed

            % TOP FRAME EDGE
            if 1 < size(index_topEdge,2)

                % calculate criteria for Area and Centroid
                centroids= zeros(size(index_topEdge));
                areas_     = zeros(size(index_topEdge));
                for i = 1:size(index_topEdge,2)     % for each object on the top frame edge
                    sum_x = sum(imgLabeled==index_topEdge(i),1);
                    areas_(i) = sum(sum_x);
                    centroids(i) = abs(sum_x*Y/areas_(i));
                end
                % Estimate scores for each criteria
                score_centroids = 1-(centroids-min(centroids))/(max(centroids)-min(centroids));
                score_areas     =   (areas_    -min(areas_))    /(max(areas_)    -min(areas_));

                score = score_centroids*weight_score + score_areas*(1-weight_score);
                index_top = index_topEdge(score ==max(score));
            else
                index_top = index_topEdge;
            end

            % BOTTOM FRAME EDGE
            if 1 < size(index_bottomEdge,2)

                % calculate criteria for Area and Centroid
                centroids = zeros(size(index_bottomEdge));
                areas_     = zeros(size(index_bottomEdge));
                for i = 1:size(index_bottomEdge,2) % for each object on the bottom frame edge
                    sum_x = sum(imgLabeled==index_bottomEdge(i),1);
                    areas_(i) = sum(sum_x);
                    centroids(i) = abs(sum_x*Y/areas_(i));
                end
                % Estimate scores for each criteria
                score_centroids = 1-(centroids-min(centroids))/(max(centroids)-min(centroids));
                score_areas     =   (areas_    -min(areas_))    /(max(areas_)    -min(areas_));

                score = score_centroids*weight_score + score_areas*(1-weight_score);
                index_bottom = index_bottomEdge(score ==max(score));
            else
                index_bottom = index_bottomEdge;
            end

            if isempty(index_top) || isempty(index_bottom), error('cannulas:find:canNOT','No cannula found at bottom frame edge!'); end

            obj.cannulas.imgBinary = {imgLabeled == index_top;imgLabeled == index_bottom};
        end
        function improve_cannulasShape(obj)
            %IMPROVE_CANNULASSHAPE eliminate defects on the cannula shape

            for i = 1:2 % for each cannula
                if isempty(obj.cannulas.imgBinary{i}), continue; end

                %% generate rectangular mask for the cannula
                X = sum(obj.cannulas.imgBinary{i},2);
                Y = sum(obj.cannulas.imgBinary{i},1);
                int_x = get_interval(X > mean(X));
                int_y = get_interval(Y > mean(Y));

                X = false(size(X));
                Y = false(size(Y));
                if int_x(1) == 1
                    % exclude all objects without pixels at upper frame edge.
                    tmp_edge = arrayfun(@(y1,y2) any(obj.cannulas.imgBinary{i}(1,y1:y2)),int_y(:,1),int_y(:,2));
                    X(1:int_x(1,2)) = true;
                else
                    % exclude all objects without pixels at lower frame edge.
                    tmp_edge = arrayfun(@(y1,y2) any(obj.cannulas.imgBinary{i}(end,y1:y2)),int_y(:,1),int_y(:,2));
                    X(int_x(1,1):end) = true;
                end
                int_y = int_y(tmp_edge,:); % exclusion

                % merge remaining intervalls
                int_y = [int_y(1,1),int_y(end,2)];
                Y(int_y(1,1):int_y(1,2)) = true;

                % generate mask finally
                cannulas_mask  = false(obj.frames.size);
                cannulas_mask(X,Y) = true;

                % apply mask to binary image of the cannula
                obj.cannulas.imgBinary{i} = obj.cannulas.imgBinary{i} & cannulas_mask;
            end
        end
        function set_cannulasExtrema(obj)
            %SET_CANNULASEXTREMA set extremal quadrangle edge points of each cannula

            obj.cannulas.roots = zeros(2);
            obj.cannulas.tops  = zeros(2);
            for i= 1:2 % for each cannula type
                if isempty(obj.cannulas.imgBinary{i}), continue; end

                % estimate cannula's centroid
                [Y, X ] = meshgrid(1:obj.frames.size(2),1:obj.frames.size(1));
                y_mid = (obj.cannulas.imgBinary{i}(:)'*Y(:))/sum(obj.cannulas.imgBinary{i}(:));
                x_mid = (obj.cannulas.imgBinary{i}(:)'*X(:))/sum(obj.cannulas.imgBinary{i}(:));

                koord = complex(X(obj.cannulas.imgBinary{i}(:)),Y(obj.cannulas.imgBinary{i}(:)));

                % sort pixel coordinates by their distance to the cannula centroid
                tmp = sortrows([koord,                                         ... % absolute pixel coordinate of the cannula
                                abs(koord-(x_mid+y_mid*1i)),                   ... % distance between cannula pixel and  centroid
                                ceil(angle(koord-(x_mid+y_mid*1i))/(pi/2)+2)], ... % number of the quadrant relatived to the cannula centroid
                               2);
                tmp1 = arrayfun(@(i) find(tmp(:,3)==i,1,'last'),(1:4)','UniformOutput',false); % find extrem cannula points in each quadrant
                if any(cellfun(@isempty,tmp1)), error('record_EXP:set_cannulasExtrema:noDetected','Cannula type %i could not be detected !',i); end
                tmp = tmp(cell2mat(tmp1),1); % reduce to extrem

                % set extremal points as cannula properties
                if i == 1
                    obj.cannulas.tops(i,1:2) = tmp(2:3,1).';
                    obj.cannulas.roots(i,1:2)  = tmp([4,1],1).';
                else
                    obj.cannulas.tops(i,1:2) = tmp([4,1],1).';
                    obj.cannulas.roots(i,1:2)  = tmp(2:3,1).';
                end
            end
        end
        function set_cannulasDiameter(obj)
            %SET_CANNULASDIAMETER estimate cannulas' diameter

            obj.cannulas.diameter = zeros(2,1);
            for i =1:2
                d = sum(obj.cannulas.imgBinary{i},2);
                obj.cannulas.diameters(i) = mean(d(any(d,2)));
            end
        end
        function set_cannulasLength(obj)
            %SET_CANNULASLENGTH estimate cannulas' length

            obj.cannulas.diameter = zeros(2,1);
            for i =1:2
                obj.cannulas.length(i) = max(sum(obj.cannulas.imgBinary{i},1));
            end
        end
        function set_cannulasTopsM(obj)
            %SET_CANNULASTOPSM estimate cannula mean top

            obj.cannulas.tops_m = mean(obj.cannulas.tops,2);
        end
        function set_px2mm(obj)
            %SET_PX2MM estimate the ratio of millimeter per pixel length

            obj.frames.px2mm = obj.cannulas.diameter_mm/mean(obj.cannulas.diameters);
        end

        %% DROPS
        function get_dropCenters(obj)
            %GET_DROPCENTERS estimate and adjust drop center through drop edge

            [Y, X] = meshgrid(1:obj.frames.size(2),1:obj.frames.size(1));
            koords = X + 1i*Y;
            i_0 = obj.get_firstDropFrame();
            for i = i_0:obj.frames.numberFrames
                BW = obj.cleared_cannulasImgBinary(i);
                label = bwlabel(BW);
                for type = 1:2
                    if i ~= i_0
                        % estimate position of the drop center based on the previous detected drop
                        obj.drops(type).centers(i)    = obj.drops(type).centers(i-1) ;%+ obj.drops(type).velocities(i-1) ;
                        if obj.drops(type).validates(i-1)
                            % use radii values from pervious detected drop
                            obj.drops(type).mean_radii(i) = obj.drops(type).mean_radii(i-1);
                            obj.drops(type).min_radii(i)  = obj.drops(type).min_radii(i-1);
                            obj.drops(type).max_radii(i)  = obj.drops(type).max_radii(i-1);
                        end
                    end
                    % reduce points to edge points only
                    koords_on_dropCenter = koords(BW(:))- obj.drops(type).centers(i);
                    pixel_props = [  koords(BW(:)),                 ... % index 1
                                     label(BW(:)) ,                 ... % index 2
                                     koords_on_dropCenter,          ... % index 3
                                     abs(koords_on_dropCenter)      ... % index 4
                                     ... angle(koords_on_dropCenter)... % index 5
                                  ];
                    obj.drops(type).edgePoints{i} = pixel_props(:,1);

                    % eliminate grouped pixels which are not completely inside 1.5 of then maximal radii
                    tmp = ~ismember(pixel_props(:,2), unique(pixel_props(obj.drops(type).max_radii(i)*1.5<pixel_props(:,4),2)));
                    pixel_props = pixel_props(tmp,:);

                    if any(tmp)

                        % reduce pixels via twice the standard deviation of the distance between drop pixel and estimated center
                        tmp = ~( mean(pixel_props(:,4))+3*std(pixel_props(:,4)) < pixel_props(:,4) );
                        pixel_props = pixel_props(tmp,:);

                        % sort pixel by angle
                        pixel_props  = sortrows([pixel_props, angle(pixel_props(:,3))],5);
                        pixel_props = [flipud(pixel_props(real(pixel_props(:,5)) < 0,:)); pixel_props( 0 <= real(pixel_props(:,5)),:)]; % negativ values are not considered by the sortrows function

                        % periodic continuation of the pixels' list
                        pixel_props = [pixel_props; pixel_props(1,:)];
                        pixel_props(end,5) = pixel_props(end,5)+2*pi;

                        % local major search with defined angle range
                        [~, index] = findpeaks_distance(pixel_props(:,5), pixel_props(:,4), acos(1-2/obj.drops(type).mean_radii(i)^2));
                        pixel_props = pixel_props(index,:); % reduce to point with maximal radius inside the angle range

                        obj.drops(type).edgePoints{i} = pixel_props(:,1);
                        % close all, h = figure(); imshow(obj.frames.imgFiltered{i},[]),hold on; plot(imag(obj.drops(type).edgePoints{i}),real(obj.drops(type).edgePoints{i}),'.m'), uiwait(h);

                        % adjustment of the estimated drop center by the detected drop edge through weights centroid determination
                        % get angular difference for each edge points between their two adjacent edge points
                        diffAngle2neighbours =    diff([pixel_props(end,5)-2*pi; pixel_props(:,5)])  ...
                                                + diff([pixel_props(:,5); pixel_props(1,5)+2*pi]);

                        % adjust the estimated center drop
                        dS = diffAngle2neighbours'*pixel_props(:,3)/(4*pi); % Korrektur des Tropfenmittelpunktes
                        obj.drops(type).centers(i) = obj.drops(type).centers(i) + dS;

                        if i ~= 1
                            obj.drops(type).velocities(i) = obj.drops(type).centers(i)-obj.drops(type).centers(i-1);
                        end

                        % adjust radius of the drop
                        radiiEdgePoints = abs(pixel_props(:,3) - dS);
                        obj.drops(type).mean_radii(i) = mean(radiiEdgePoints);
                        obj.drops(type).min_radii(i)  = min(radiiEdgePoints);
                        obj.drops(type).max_radii(i)  = max(radiiEdgePoints);
                        obj.drops(type).validates(i)  = true;
                    end
                    if i == i_0+1
                        % copy drop data to the begin of the record
                        obj.drops(type).centers(1:i_0) = obj.drops(type).centers(i);
                        obj.drops(type).mean_radii(1:i_0) = obj.drops(type).mean_radii(i);
                        obj.drops(type).min_radii(1:i_0)  = obj.drops(type).min_radii(i);
                        obj.drops(type).max_radii(1:i_0)  = obj.drops(type).max_radii(i);
                    end
                end
            end
        end
        function get_DropsEdgePoints(obj)
            %GET_DROPSEDGEPOINTS

            warning('record_EXP:get_DropsEdgePoints:TODO', 'Function ''obj.get_DropsEdgePoints()'' is in development.')
        end
        function out = get_firstDropFrame(obj)
            %GET_FIRSTDROPFRAME get the frame number when the lower drop has been generated

            out = 1;
            numberOFpixel_ROI = cell(2,1);
            types = 2;
            for type = types
                int = ( -ceil(obj.drops(type).max_radii(1)):ceil(obj.drops(type).max_radii(1)) );
                koord_x = round(real(obj.drops(type).centers(1))) - int ;
                koord_y = round(imag(obj.drops(type).centers(1))) - int ;
                numberOFpixel_ROI{type,1} = cellfun(@(img) sum(sum(img( koord_x , koord_y ) ) ) ,  obj.frames.imgBinary );
                numberOFpixel_ROI{type,1} = numberOFpixel_ROI{type,1} < mean(numberOFpixel_ROI{type,1});
            end

            int = get_interval(numberOFpixel_ROI{2,1});

            if int(1,1) == 1
                out = int(1,2)+1;
                obj.add_log('Time shifting is needed because the lower drop is not shown in the first frame.');
            end
        end
        function out = cleared_cannulasImgBinary(obj,i)
            %CLEARED_CANNULASIMGBINARY eliminate cannulas and disturbing artefacts from binary image

            out = obj.frames.imgBinary{i};
            lin = 1:obj.frames.size(1);

            tmp = lin < real(obj.cannulas.tops_m(1))+1 | real(obj.cannulas.tops_m(2))-1 <= lin ;
            out(tmp,:) = false;

            % eliminate disturbing artefacts at left and right frame edge (rinsing cannula)
            label = bwlabel(out,8);
            tmp2 = unique([label(:,1);label(:,end)]);
            out = ~ismember(label,tmp2);
        end
        function smooth_dropsDynamic(obj)
            %SMOOTH_DROPSDYNAMIC smooth drop trajectory by moving average, see smoothing_range property for adjustments

            tmp = ones(obj.smoothing_range*2+1,1);
            for type = 1:2
                int = get_interval(obj.drops(type).validates);
                if isempty(int), continue; end
                int = int(1,1) - 1 + get_interval(abs(diff(obj.drops(type).centers(int(1,1):int(1,2)),2,1)) < 10); %
                if isempty(int), continue; end
                int = int(obj.smoothing_range*2+1 <= diff(int,1,2),:);
                for i = 1:size(int,1)
                    obj.drops(type).centers(int(i,1):int(i,2)) = conv2(obj.drops(type).centers(int(i,1):int(i,2)), tmp, 'same') ...
                                                               ./conv2(ones(int(i,2)-int(i,1)+1,1), tmp, 'same');
                end
            end

            % estimate new drop velocity
            obj.get_dropsVelocities();

            % estimate new drop radii
            obj.get_dropsRadii();
        end
        function get_dropsVelocities(obj)
            %GET_DROPSVELOCITIES estimate drop velocities via central difference quotient

            tmp = [1;0;-1];
            obj.drops(1).velocities = conv2(obj.drops(1).centers, tmp, 'same')/2;
            obj.drops(2).velocities = conv2(obj.drops(2).centers, tmp, 'same')/2;

            % adjustements of the limits
            for i = 1:2
                obj.drops(i).velocities(1)   =  obj.drops(i).centers(1) - obj.drops(i).centers(2);
                obj.drops(i).velocities(end) =  obj.drops(i).centers(end-1) - obj.drops(i).centers(end);
            end
        end
        function get_dropsRadii(obj)
            %GET_DROPSRADII estimate drop radii from given drop centers and drop edge

            % get raddii to the drop based on the arccording drop centers
            radii_center_1 = cellfun(@(p,c) abs(p - c) ,obj.drops(1).edgePoints,num2cell(obj.drops(1).centers),'UniformOutput',false);
            radii_center_2 = cellfun(@(p,c) abs(p - c) ,obj.drops(2).edgePoints,num2cell(obj.drops(2).centers),'UniformOutput',false);

            % if coordinates are empty values are set to 0
            radii_center_1(cellfun(@isempty, radii_center_1)) = {0};
            radii_center_2(cellfun(@isempty, radii_center_2)) = {0};

            % get statistic radii value as drop properties
            obj.drops(1).min_radii = cellfun(@min,radii_center_1);
            obj.drops(2).min_radii = cellfun(@min,radii_center_2);

            obj.drops(1).max_radii = cellfun(@max,radii_center_1);
            obj.drops(2).max_radii = cellfun(@max,radii_center_2);

            obj.drops(1).mean_radii = cellfun(@mean,radii_center_1);
            obj.drops(2).mean_radii = cellfun(@mean,radii_center_2);
        end

        %% MISC
        function add_log(obj,msg)
            %ADD_LOG add message to log for replicability when methods are used with unusual adjustments

            if ~ischar(msg), error('record_EXP:add_log:MSGisntchar','Wrong type of given msg, required type: char, not %s!',class(msg)); end
            obj.log = [obj.log;{msg}];
            if obj.display_logs,  disp(msg); end
        end
        function load_cameraSystem(obj, CAM_ID)
            %LOAD_CAMERASYSTEM load properties of the used camera system. Properties can be added as

            if exist('CAM_ID','var')
                % set the max. relative standard noise of the camera system
                obj.frames.maxRelStdNoises = obj.cam_syst{CAM_ID,2};
            end
        end
        function load_parameters(obj)
            %LOAD_PARAMETERS load parameter values

            if ~exist([obj.path,obj.parameters_mat_file],'file')
                if ~exist([obj.path,obj.parameters_csv_file],'file')
                    error('record_EXP:paraCSVNotFound','paramters.csv file does not exist!\n PFAD: %s',obj.path);
                else
                    obj.load_parametersCSV();
                end
            else
                obj.load_parametersMAT();
            end
        end
        function load_parametersCSV(obj)
            %LOAD_PARAMETERSCSV load parameter values from given CSV-file

            fid = fopen([obj.path,obj.parameters_csv_file],'r','n','windows-1252');
               str = textscan(fid,'%s','delimiter','\n');
            fclose(fid);
               str = regexprep (str{1},'(?<=\d)\,(?=\d)', '.');
               str = regexprep (str,' & ', '_');
               str = regexprep (str,'  ', ' ');
               str = regexprep (str,'\*','');
               str = regexprep (str,'?', '');
               str = regexprep (str,'(?<=[a-zA-Z])[\, ](?=[a-zA-Z])', '_');

               % eliminate headlines
               str = str(cellfun(@isempty,regexp(str,'--')) & cellfun(@isempty,regexp(str,'D:\\')));

               % split parameters data into names and according values
               str = regexp(str(1:end),';', 'split');
               str = [cellfun(@(a) a{:,1},str,'UniformOutput',false), cellfun(@(a) a{:,2},str,'UniformOutput',false)];

               % get the units of the values
               c = cell(size(str(:,1)));
               tmp = ~cellfun(@isempty,regexp(str(:,1),'['));
               c(~tmp) = cellfun(@(a) [a, {''}], str(~tmp,1),'UniformOutput',false);
               c(tmp) = cellfun(@(c, i) { strtrim(c(1:i-1)) ,strtrim(c(i:end))} , str(tmp,1), regexp(str(tmp,1),'['),'UniformOutput',false );
               num = cellfun(@str2double,str(:,2));
               str(~isnan(num),2)= num2cell(num(~isnan(num)));
               fin = [cellfun(@(a) regexprep(a{:,1},'\.',''),c,'UniformOutput',false), cellfun(@(a) a{:,2},c,'UniformOutput',false),str(:,2)];

               % rename variable (coused by changing parameter.csv)
               rename_list = {% new                  old               ...
                               'd_disp_top_target', 'd_oil_top';       ...
                               'd_disp_bot_target', 'd_oil_bottom'     ...
                             };
               for i = 1:size(rename_list,1)
                  tmp = ismember(fin(:,1),rename_list{i,2});
                  if 1 == sum(tmp)
                      fin(tmp,1) = rename_list(i,1);
                  end
               end

               % convert into struct-array
               obj.parameters = cell2struct(fin(:,3),fin(:,1),1);
               obj.parameters.units = cell2struct(fin(:,2),fin(:,1),1);
        end
        function save_parametersMAT(obj)
            %SAVE_PARAMETERSMAT save parameter values to MAT-file

            if ~isempty(obj.parameters)
                parameters = obj.parameters;  %#ok<NASGU>
                save([obj.path,obj.parameters_mat_file],'parameters');
            else
                warning('record_EXP:saveParametersMat:isempty','parameters is empty');
            end
        end
        function load_parametersMAT(obj)
            %LOAD_PARAMETERSMAT load paramter values from converted MAT-file

            load([obj.path,obj.parameters_mat_file]);
            if exist('parameters','var')
                obj.parameters = parameters;  %#ok<CPROP>
            elseif exist('para','var') % old variable name
                obj.parameters = para;
            else
                error('record_EXP:loadParametersMat:noPara','In parameters.mat ''parameters'' is not defined!');
            end
        end
        function set_hashID(obj)
            %SET_HASHID generate hashID from the first image of the record as unique identifier

            if ~isfield(obj.parameters,'hashID')
                mddigest   = java.security.MessageDigest.getInstance('MD5');
                filestream = java.io.FileInputStream(java.io.File(obj.frames.filepaths{1}));
                digestream = java.security.DigestInputStream(filestream,mddigest);
                while(digestream.read() ~= -1), end
                obj.parameters.hashID = reshape(dec2hex(typecast(mddigest.digest,'uint8'))',1,[]);
            end
        end
        function load_manualValues(obj)
            %LOAD_MANUALVALUES load the manual values for this record based on given hashID
            global manualValue_file;
            if ~exist(manualValue_file,'file'), error('record_EXP:load_manualValues:FileNotExist','manualValue_file does not exist\n %s.\n If necessary you can create it with record_EXP.create_manualValuesFile()',manualValue_file); end

            load(manualValue_file);
            if ~exist('manualValues','var'), error('record_EXP:load_manualValues:noDATA','manualValues does not exit in file:\n %s.',manualValue_file); end
            tmp = ismember({manualValues(:).hashID},obj.parameters.hashID);  %#ok<NODEF>
            if sum(tmp) == 1
                obj.manualValues.bottomDetachmentTime = manualValues(tmp).bottomDetachmentTime;
                obj.manualValues.collisionTime        = manualValues(tmp).collisionTime;
                obj.manualValues.coalescenceTime      = manualValues(tmp).coalescenceTime;
                obj.manualValues.repulsionTime        = manualValues(tmp).repulsionTime;
                obj.manualValues.crop_coordinates     = manualValues(tmp).crop_coordinates;
                obj.manualValues.man_threshold        = manualValues(tmp).man_threshold;
            end
        end
    end
end