clear; close all; clc;
global target;
global safe_limit;
global base_init;
global MARKER;
global LINE;


%% Plot step-by-step prediction and past positions of needle tip
%% Plot also the calculated and performed control outputs for each step

stop_step = false; % true for stopping at each insertion step (hit keyboard)

%% Load Dataset
trial = 04;
folder = 'placement_accuracy_invitro';
name = 'phantomA_closedloop_';
load(strcat(folder,'/',name,num2str(trial,'%2.2d'),'.mat'));

%% Configure simulation
S = size(up,2);  % Number of steps
safe_limit = 6;
base_init = base(:,1);

% Configure plot
INSERTION_STEP = -5;
MARKER = 10;
LINE = 1.5;
LEGLOC = 'northeast';

%% Loop key
k_key = find(key); % Samples when key was pressed
ns = length(k_key);
T = zeros(1,ns);

%Invert negative values
target = [-target(1); -target(2); target(3)];
base_init = [-base_init(1); -base_init(2); base_init(3)];

for i=1:ns
    k = k_key(i);   %sample
    Jc_step{i} = J{k}(1:3,:);
    Z_step(i,:) = Z(:,k);
    X_step(i,:) = X(:,k);
    cmd_step(i,:) = cmd(:, k);
    
    % Invert negative values
    Z_step(i,1) = -Z_step(i,1);
    Z_step(i,2) = -Z_step(i,2);
    X_step(i,1) = -X_step(i,1);
    X_step(i,2) = -X_step(i,2);
    cmd_step(i,1) = -cmd_step(i,1);
    cmd_step(i,2) = -cmd_step(i,2);
    if i<ns
        yp{i}(:,1) = -yp{i}(:,1);
        yp{i}(:,2) = -yp{i}(:,2);
        up{i}(:,1) = -up{i}(:,1);
    end
    
    if size(yp{1},2) == 3
        err_step(i,:) = Z_step(i,1:3)'-target;
    else
        err_step(i,:) = Z_step(i,1:5)'-[target;0;0];
    end
end

%% Saturation percentage
% Get saturation limits
x_lim = [base_init(1)+safe_limit base_init(1), base_init(1)-safe_limit base_init(1)];
z_lim = [base_init(3)+safe_limit base_init(3), base_init(3)-safe_limit base_init(3)];
x_max = max(x_lim)*ones(ns,1);
x_min = min(x_lim)*ones(ns,1);
z_max = max(z_lim)*ones(ns,1);
z_min = min(z_lim)*ones(ns,1);

% Calculate saturation
sat_x = or(X_step(:,1)>=x_max, X_step(:,1)<=x_min);
sat_z = or(X_step(:,3)>=z_max, X_step(:,3)<=z_min);
sat = or(sat_x,sat_z);
sat_x_perc = sum(sat_x)/length(sat_x);
sat_z_perc = sum(sat_z)/length(sat_z);
sat_perc = sum(sat)/length(sat);

%% Calculate final results
final_results = sprintf('*** RESULTS ***\n\n');
final_results = [final_results, sprintf('Error X [mm] = %0.2f\n', abs(err_step(end,1)))];
final_results = [final_results, sprintf('Error Z [mm] = %0.2f\n', abs(err_step(end,3)))];
final_results = [final_results, sprintf('Trajectory error [mm] = %0.2f\n', sqrt(err_step(end,1)^2+err_step(end,3)^2))];
final_results = [final_results, sprintf(' \n')];
final_results = [final_results, sprintf('Saturation X (horizontal) [0-1] = %0.2f\n', sat_x_perc)];
final_results = [final_results, sprintf('Saturation Y (vertical) [0-1] = %0.2f\n', sat_z_perc)];
final_results = [final_results, sprintf('Control saturation [0-1] = %0.2f\n', sat_x_perc)];
fprintf(final_results);

depth = [0:5:100];

% First insertion step
f = figure;
f.Position = [0 0 650 600];
f.PaperOrientation = 'landscape'; 

subplot(4,1,1);
plot(depth(1:2), Z_step(1:2,1), '.-m','LineWidth',LINE, 'MarkerSize',MARKER);
plot_target_X('--k');
title('Needle tip trajectory');
ylabel('Tip X [mm]'); legend('executed', 'target', 'Location', LEGLOC);
xlim([0 100]);

subplot(4,1,3);
plot(depth(1:2), X_step(1:2,1), '.-m','LineWidth',LINE, 'MarkerSize',MARKER);
plot_safe_limit_X('k');
title('Needle guide displacement');
ylabel('Guide X [mm]'); legend('executed', 'safe limit', 'Location', LEGLOC);
xlim([0 100]);

subplot(4,1,2);
plot(depth(1:2), Z_step(1:2,3), '.-m','LineWidth',LINE, 'MarkerSize',MARKER);
plot_target_Z('--k');
ylabel('Tip Z [mm]'); legend('executed', 'target', 'Location', LEGLOC);
xlim([0 100]);

subplot(4,1,4);
plot(depth(1:2), X_step(1:2,3), '.-m','LineWidth',LINE, 'MarkerSize',MARKER);
plot_safe_limit_Z('k');
xlabel('Length [mm]'),ylabel('Guide Z [mm]'); legend('executed', 'safe limit','Location', LEGLOC);
xlim([0 100]);
if stop_step
    fprintf('*** Press keyboard for next insertion step***\n');
    pause;
else
    pause(0.01);
end

past_u = [X_step(1:2,1) X_step(1:2,3)];
past_y = Z_step(1:2,1:3);

future_u = [past_u(2,:); up{1}];
future_y = [past_y(2,:); yp{1}(:,1:3)];

%% Loop prediction steps
for i=2:S    
    if (i>=ns/2)
        LEGLOC = 'northwest';
    end
    H = min(S-i+1,size(up{i},1));
    subplot(4,1,1);
    plot(depth(1:i), past_y(:,1), '.-m', depth(i:H+i), future_y(:,1), '.-b', 'LineWidth',LINE, 'MarkerSize',MARKER);
    plot_target_X('--k');
    title('Needle tip trajectory');
    ylabel('Tip X [mm]'); legend('executed','prediction',  'target', 'Location', LEGLOC);
    xlim([0 100]);
    
    subplot(4,1,3);
    plot(depth(1:i), past_u(:,1), '.-m', depth(i:H+i), future_u(:,1), '.-b', 'LineWidth',LINE, 'MarkerSize',MARKER);
    plot_safe_limit_X('k');
    title('Needle guide displacement');
    ylabel('Guide X [mm]'); legend('executed','prediction',  'safe limit', 'Location', LEGLOC);
    xlim([0 100]);
    
    subplot(4,1,2);
    plot(depth(1:i), past_y(:,3), '.-m', depth(i:H+i), future_y(:,3), '.-b', 'LineWidth',LINE, 'MarkerSize',MARKER);
    plot_target_Z('--k');
    ylabel('Tip Z [mm]'); legend('executed','prediction',  'target', 'Location', LEGLOC);
    xlim([0 100]);
    
    subplot(4,1,4);
    plot(depth(1:i), past_u(:,2), '.-m', depth(i:H+i), future_u(:,2), '.-b', 'LineWidth',LINE, 'MarkerSize',MARKER);
    plot_safe_limit_Z('k');
    xlabel('Length [mm]'),ylabel('Guide Z [mm]'); legend('executed', 'prediction', 'safe limit', 'Location', LEGLOC);
    xlim([0 100]);

    if (i~=S)
        futH = min(S-i ,size(up{i},1));
        future_u = [X_step(i+1,1), X_step(i+1,3); up{i}(1:futH,:)];
        future_y = [Z_step(i+1,1:3); yp{i}(1:futH,1:3)];
    end
    past_u = [past_u; [X_step(i+1,1) X_step(i+1,3)]];
    past_y = [past_y; Z_step(i+1,1:3)];

    if stop_step
        fprintf('*** Press keyboard for next insertion step***\n');
        pause;
    else
        pause(0.01);
    end
   
end

subplot(4,1,1);
plot(depth, Z_step(:,1), '.-m', 'LineWidth',LINE, 'MarkerSize',MARKER);
plot_target_X('--k');
title('Needle tip trajectory');
ylabel('Tip X [mm]'); legend('executed', 'target', 'Location', LEGLOC);
xlim([0 100]);

subplot(4,1,3);
plot(depth, X_step(:,1), '.-m', 'LineWidth',LINE, 'MarkerSize',MARKER);
plot_safe_limit_X('k');
title('Needle guide displacement');
ylabel('Guide X [mm]'); legend('executed', 'safe limit', 'Location', LEGLOC);
xlim([0 100]);

subplot(4,1,2);
plot(depth, Z_step(:,3), '.-m', 'LineWidth',LINE, 'MarkerSize',MARKER);
plot_target_Z('--k');
ylabel('Tip Z [mm]'); legend('executed', 'target', 'Location', LEGLOC);
xlim([0 100]);

subplot(4,1,4);
plot(depth, X_step(:,3), '.-m', 'LineWidth',LINE, 'MarkerSize',MARKER);
plot_safe_limit_Z('k');
xlabel('Length [mm]'),ylabel('Guide Z [mm]'); legend('executed', 'safe limit', 'Location', LEGLOC);
xlim([0 100]);

 
function plot_target(line)
    global target;
    plot_target_X(line);
    plot_target_Z(line);
end

function plot_target_X(line)
    global target;
    global LINE;
    yline(target(1), line, 'LineWidth', LINE);
    ylim([target(1)-15, target(1)+15])
end

function plot_target_Z(line)
    global target;
    global LINE;
    yline(target(3), line, 'LineWidth', LINE);
    ylim([target(3)-15, target(3)+15])
end

function plot_safe_limit(line)
    global LINE;
    plot_safe_limit_X(line, 'LineWidth', LINE);
    plot_safe_limit_Z(line, 'LineWidth', LINE);
end

function plot_safe_limit_X(line)
    global base_init;
    global safe_limit;
    global LINE;
    yline(base_init(1)+safe_limit,line, 'LineWidth', LINE);
    yline(base_init(1)-safe_limit,line, 'LineWidth', LINE);
    ylim([base_init(1)-safe_limit-3, base_init(1)+safe_limit+3])

end

function plot_safe_limit_Z(line)
    global base_init;
    global safe_limit;
    global LINE;
    yline(base_init(3)+safe_limit,line, 'LineWidth', LINE);
    yline(base_init(3)-safe_limit,line, 'LineWidth', LINE);
    ylim([base_init(3)-safe_limit-3, base_init(3)+safe_limit+3])
end