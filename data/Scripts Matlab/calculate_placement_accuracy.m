clear all;
close all;
clc;

%% Calculate mean and std deviation values for all insertions from a given trial
%% Select appropriate dataset (with manual selection of intervals for trial and extra)
%% Select if no_compensation (00) or compensated trials
%% Select if MPC3 algorithm (with angles) or previous version

folder = 'placement_accuracy_invitro';
trial = 1:5;

% Compare open-loop and closed-loop for phantom A
get_trial_errors(folder, 'A', trial, true);
get_trial_errors(folder, 'A', trial, false);

% Compare with closed-loop for phantom A without skin
get_trial_errors(folder, 'A_noskin', trial, false);



function get_trial_errors(folder, phantom, trial, open_loop)

    MARKER = 10;
    LINE = 1.5;

    %% Select Dataset

    name = ['phantom',phantom,'_'];
    mpc3 = true;    %True = Used MPC version with angles

    safe_limit = 5.9999;

    %% Set figure
    f1 = figure();
    f1.Position = [0 0 650 146];
    f1.PaperOrientation = 'landscape'; 

    %% Initialize vectors
    if open_loop == true
        name = [name, 'openloop_'];
    else
        name = [name, 'closedloop_'];
    end
    N = length(trial);
    if mpc3 == true
        err_step = zeros(5,N);
        [err_angle_v, err_angle_h] = deal(zeros(N,1));
    else
        err_step = zeros(3,N);
    end
    [err_x, err_y, err_z, err_2d, err_3d] = deal(zeros(N,1));   

    fprintf('**********\nName: %s\n', name);

    %% Get error for each trial
    for j=1:N
        %% Load data
        load(strcat(folder,'/',name,num2str(trial(j),'%2.2d'),'.mat'));

        %% Loop all steps
        k_key = find(key); % Samples when key was pressed
        ns = length(k_key);
        depth = 5*(0:ns-1);
        % Get saturation limits
        x_lim = [base_init(1)+safe_limit base_init(1), base_init(1)-safe_limit base_init(1)];
        z_lim = [base_init(3)+safe_limit base_init(3), base_init(3)-safe_limit base_init(3)];
        x_max = max(x_lim)*ones(ns,1);
        x_min = min(x_lim)*ones(ns,1);
        z_max = max(z_lim)*ones(ns,1);
        z_min = min(z_lim)*ones(ns,1);

        X_step = zeros(ns,3);
        for i=1:ns
            k = k_key(i);   %sample
            X_step(i,:) = X(:,k);
            Z_step(i,:) = Z(:,k);
            err_traj(i,:) = sqrt((Z(1,k)-Z(1,1))^2 + (Z(3,k)-Z(3,1))^2);
        end

        % Plot figures

        figure(f1.Number)
        ax = gca;
        switch phantom
            case 'b3'
                pair = 0;
            case 'b2'
                pair = 1;
            case 'b1'
                pair = 2;
            otherwise 
                pair = 3;
        end
        if(open_loop)
            ax.ColorOrderIndex = 1+2*pair;
        else
            ax.ColorOrderIndex = 2+2*pair;
        end
             
        plot(depth, err_traj, '.-','LineWidth',LINE, 'MarkerSize',MARKER);
        hold on
        title('Trajectory error'),xlabel('Insertion depth [mm]'),ylabel('Error [mm]');
        

        %% Saturation percentage
        sat_x = or(X_step(:,1)>=x_max, X_step(:,1)<=x_min);
        sat_z = or(X_step(:,3)>=z_max, X_step(:,3)<=z_min);
        sat = or(sat_x,sat_z);
        sat_x_perc(j) = sum(sat_x)/length(sat_x);
        sat_z_perc(j) = sum(sat_z)/length(sat_z);
        sat_perc(j) = sum(sat)/length(sat);

        %% RMS error
        rms_err_2d(j) = rms(err_traj);

        %% Get final error
        if size(yp{1},2) == 3
            err_step = Z(1:3,k_key(ns)) - target;
        else
            err_step = Z(1:5,k_key(ns)) - [target;0;0];
        end    
        err_x(j) = abs(err_step(1,end));
        err_y(j) = abs(err_step(2,end));
        err_z(j) = abs(err_step(3,end));

        err_2d(j) = sqrt(err_step(1,end)^2+err_step(3,end)^2);
        err_3d(j) = sqrt(err_step(1,end)^2+err_step(2,end)^2+err_step(3,end)^2);

        if size(err_step,2) == 5
            err_angle_v(j) = abs(err_step(4,end));
            err_angle_h(j) = abs(err_step(5,end));
        end

        if open_loop == true
            fprintf('Trial %2.2d \t Err 2D[mm] = %0.4f \t RMS Err 2D[mm] = %0.4f\n', trial(j), err_2d(j), rms_err_2d(j));
        else
            fprintf('Trial %2.2d \t Err 2D[mm] = %0.4f \t Sat[0-1] = %0.4f \t RMS Err 2D[mm] = %0.4f\n', trial(j), err_2d(j), sat_perc(j), rms_err_2d(j));
        end
    end


    fprintf('\nTotal error 2D[mm] = %0.4f +- %0.4f\n', mean(err_2d), std(err_2d));
    fprintf('Total RMS error 2D[mm] = %0.4f +- %0.4f\n', mean(rms_err_2d), std(rms_err_2d));
    if open_loop == false
        fprintf('Controller saturation[0-1] = %0.4f +- %0.4f\n', mean(sat_perc), std(sat_perc));
    end
end
