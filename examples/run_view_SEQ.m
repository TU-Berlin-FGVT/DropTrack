%% Getting Started with view_SEQ
%
% Because you may have not a parameterSets, we just generate them before
% (Preparation Block, see below). Then you are able to use this class.
% The final usage you find below.
%
% Usually you have already evaluated your records with record_EXP class and
% saved the exp_data.mat in the records directory.
%
%
%% Preparation
if any(regexp(cd(),'examples$','once'))
    cd('..');
end
config;

% generate exp_data.mat
path = osPATH('.\examples\paraSet01\');
mkdir(path)
mkdir([path,'fig']);
path_exp_data = strcat(path, num2cell(num2str((0:9)')),'exp_data.mat');

parameters.d_disp_top_target = 1;
parameters.d_disp_bot_target = 23;
parameters.drop_edge_distance = 1;
parameters.Framerate = 10000;

data = struct( 'parameters' , parameters,      ...
               'dateTime', datestr(now),          ...
               'commitID', 'abcdefghijk',         ...
               'commitCLEAN','',                  ...
               'frames', struct('px2mm',0.08),    ...
               'log', 'GENERATED DATA FILE'       ...
             );
T = 400;
t = (1:T)';
for i = 1:size(path_exp_data,1)
    dt_rand = ceil(rand(1)*10);
    dh_rand = rand(1,1)*2;
    dv_rand = 0.5*rand(1,1)*sin(30*2*pi*t/T*rand(1,1)).*exp(-t/T);
    heights = sqrt(T)+2-(sqrt(t-dt_rand)-dh_rand);
    velocities = [diff(heights);heights(end)-heights(end-1)]-dv_rand;
    data.drops = [struct('centers', [],     'velocities', []);          ...
                  struct('centers', heights, 'velocities', velocities)  ...
                 ];
    data.Events= struct('bottomDetachmentTime', dt_rand*ones(1,1), ...
                        'collisionTime', round((T*3/4)-4*rand(1))*ones(1,1));
    if mod(i,2) == 0
        data.parameters.d_oil_top = 12;
        data.parameters.d_oil_bottom = 10;
        data.parameters.drop_edge_distance = 1;
        data.parameters.Framerate = 10000;
    else
        data.parameters = parameters;
    end
    save(path_exp_data{i},'data');
end
clear all;

%% Example Code
config;

path = osPATH('.\examples\paraSet01\');
path_output = osPATH('.\examples\paraSet01\fig\');

seq = view_SEQ(path);

% load all exp_data.mat-files in path recursive
seq.load_data();

% label all records with the same parameterSet (defined by drops' diameter
% and their distance)
parameterSets = seq.label_same_parameterSet();

% select one parameter set
selected_parameterset = parameterSets == 1;

% plot averaged drop trajectory with three graphs
[fh, data] = seq.subplot_averaged(selected_parameterset);

savefig(fh, [path_output,'seq_001.fig']);
save([path_output,'seq_001.mat'],'data');
