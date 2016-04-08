%% preparation
if any(regexp(cd(),'examples$','once'))
    cd('..');
end
config;
unzip(osPATH([ROOT,'examples\record01.zip']),osPATH([ROOT,'examples\']));
record_EXP.create_manualValueFile(); % create an empty manualValuesFile

clear all;

%% example code
config;

% process the record
path = osPATH([ROOT,'examples\record01\']); % full path to the record
rec  = record_EXP(path);
rec.run_all();

% access to processed data of the record
BW = rec.frames.imgBinary{1}; % get the first binary image frame
rec.drops(2).centers; % get all drop centers of the lower drop
rec.Events.bottomDetachmentTime; % get the restricted time interval of the detachtment event of the bottom drop

% access to processed data of saved data
load(osPATH([ROOT,'examples\record01\exp_data.mat'])); % load processed data in `data`
data.drops(1).centers; % get all drop centers of the lower drop

%% plot drop trajectory
bottDetachTIME = mean(rec.Events.bottomDetachmentTime);
collTIME = mean(rec.Events.collisionTime);

figure();
imshow(rec.frames.imgBinary{collTIME})
hold on;
plot( imag(rec.drops(2).centers) , real(rec.drops(2).centers) ,'-b')
hold off;

%% plot over time
figure();
int = (1:rec.frames.numberFrames)';
subplot(2,1,1)
hold on;
l1 = plot( bottDetachTIME*ones(2,1) , [max(real(rec.drops(2).centers));min(real(rec.drops(2).centers))],'-.r');
l2 = plot( collTIME*ones(2,1)       , [max(real(rec.drops(2).centers));min(real(rec.drops(2).centers))],'-.m');
plot( int , real(rec.drops(2).centers) ,'-b')
legend([l1,l2],{'bottom detachment event','collision event'},'southoutside');
title('height/frame of the bottom drop')
xlabel('[frame]');
ylabel('[px]');
hold off;

subplot(2,1,2)
hold on;
plot( bottDetachTIME*ones(2,1) , [max(real(rec.drops(2).velocities));min(real(rec.drops(2).velocities))],'-.r')
plot( collTIME*ones(2,1)       , [max(real(rec.drops(2).velocities));min(real(rec.drops(2).velocities))],'-.m')
plot( int , real(rec.drops(2).velocities) ,'-b')
title('rising velocity / frame of the bottom drop')
xlabel('[frame]');
ylabel('[px/frame]');
hold off;
