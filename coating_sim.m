%% ELEC 7310 - Advanced Electrodynamics I - Project 2
% Jonathan Gunasekaran
% 2026/11/15

clear all; close all;

%% Simulation Setup
% constants
epsilon_0 = 8.854e-12;
mu_0 = pi*4e-7;
c = 1/sqrt(mu_0*epsilon_0);

% set simulation frequency - eventually setup as input for frequency sweep
% analysis
f0 = 20e6;
lambda_0 = c/f0;


% define 1-D spatial array 
delta_x = 0.1;         % spatial step in meters
max_x = 1000;            % max bound on simulation space
x = 0:delta_x:max_x;    % spatial array

% define temporal array based on spatial step value to ensure numerical
% stability
tmax = 4e-6;             % simulation end time
S = 0.7;                % Courant stability factor (want less than 1 in all media)
delta_t = S*delta_x/c;  % time step
t = 0;                  % current simulation time (start at 0)

% define multi-layer stackup material parameters
load_eps_r = 5;
layer_eps_r = [2, 3, 2, 3];
num_layers = length(layer_eps_r);

% relative permittivity spatial array
eps = ones(size(x)) .* epsilon_0; 

layer_start_x = 0.6*max_x;  % start point for multi-layer interface in x

for ii=1:num_layers
    lambda = lambda_0/sqrt(layer_eps_r(ii));
    layer_end_x = round(layer_start_x+lambda/4, 2);
    
    layer_x_ind = find(x >= layer_start_x & x < layer_end_x);
    eps(layer_x_ind) = eps(layer_x_ind) .* layer_eps_r(ii);
    
    layer_start_x = layer_end_x;
end

load_x_ind = find(x >= layer_start_x);
eps(load_x_ind) = eps(load_x_ind) .* load_eps_r;

% define locations for the virtual electric field probes 
boundary_x_ind = find(diff(eps));
if isempty(boundary_x_ind)
    boundary_x_ind(1) = layer_start_x/delta_x;
end

x_dist_to_probe = max_x/10*1.5;
probe1_x_ind = boundary_x_ind(1) - x_dist_to_probe/delta_x;
probe2_x_ind = boundary_x_ind(1) + x_dist_to_probe/delta_x;

t_vec = t:delta_t:tmax;


%% Set up frequency sweep
bw = f0;  % look at band from f0-f0/2 to f0+f0/2
num_f_sweep_pts = 21;   % needs to be odd to include center freq
f = linspace(f0-bw/2, f0+bw/2, num_f_sweep_pts);
f = linspace(0, 2*f0, num_f_sweep_pts); f = f(2:end);
% f = s0e6;
% f = [10, 30] * 1e6;

refl_coef = zeros(size(f));

plot_movie = 0;
plot_probe_meas = 0;


for f_ind = 1:length(f)
    f_oper = f(f_ind);
    
    % re-initialize electric field vector - store values for past, present,
    % and future time step
    E = zeros(length(x), 3);

    % reset time-series probe measurements to 0
    E_probe1 = zeros(size(t_vec));
    E_probe2 = zeros(size(t_vec));


    %% 1-D FDTD Loop
    if plot_movie
        fdtd_movie = figure;
        frame_decimation_rate = 500;
    end

    tic
    loop_iter = 0;
    while t <= tmax
        loop_iter = loop_iter + 1;
        
        % Add source (incoming E-field) - simulate for num_cycles
        num_cycles = 8;
        if t <= num_cycles/f_oper
            E(1,3) = sin(2*pi*f_oper*t);
        else
            E(1,3) = 0;
        end
    
        % 1-D FDTD update equation
        E(2:end-1,3) = (delta_t/delta_x)^2 ./ (mu_0.*eps(2:end-1).') ...
                        .* (E(3:end,2) - 2*E(2:end-1,2) + E(1:end-2,2)) ... 
                        + 2*E(2:end-1,2) ...
                        - E(2:end-1, 1);
        
        % update stored E-field history with calculated values for the current time step
        E(:,1) = E(:,2);
        E(:,2) = E(:,3);
    
        % save E-field value at probe positions
        t_ind = floor(t/delta_t)+1;
        E_probe1(t_ind) = E(probe1_x_ind,2);
        E_probe2(t_ind) = E(probe2_x_ind,2);
    
        % increment time
        t = t + delta_t;
    
    
        % 1-D FDTD sim plot movie
        if plot_movie && mod(loop_iter, frame_decimation_rate) == 0
            plot(x, E(:,2));
            yline(0);
            xline(x(boundary_x_ind));
            xline(x([probe1_x_ind, probe2_x_ind]), '--r');
            title(["t = " + num2str(t) + "s"]);
            xlabel("x (m)");
            ylabel("E (V/m)");
            ylim([-2.5, 2.5]);
            xlim([x(1), x(end)]);
            drawnow limitrate;
        end
    end
    toc
    
    t = 0;
    
    %% Plotting of Reflection Coefficient
    sec2usec = @(x) x.*1e6;
    
    if plot_probe_meas
        figure;
        subplot(2,1,1);
        plot(sec2usec(t_vec), E_probe1, 'LineWidth',2);
        hold on;
        % yline(0);
        ylim([-1.5 1.5]);
        title("Incident + Reflected Wave")
        
        subplot(2,1,2);
        plot(sec2usec(t_vec), E_probe2, 'LineWidth',2);
        hold on;
        % yline(0);
        ylim([-1.5 1.5]);
        title("Transmitted Wave");
        
        xlabel("t (us)");
        ylabel("E (V/m)");
        sgtitle("E-Field Probe Measurements")
    end
    
    
    % filter out undesired frequency components to determine reflection
    % coefficient at current operating frequency
    
    % set window-size to match wavelength at given frequency
    lambda_oper_0 = c/f_oper;
    
    % find time difference equal to one wavelength
    wavelength_time_diff = lambda_oper_0/c;  % in free space (for reflected wave)
    
    % convert to start/end index
    t_ind_diff = floor(wavelength_time_diff/delta_t);
    
    % % naive approach to filter out signal components not at frequency of interest
    % % TODO: better approach would probably involve frequency-domain filtering
    % E0_f_oper = zeros(size(t_vec));
    % for ii=t_ind_diff+1:length(t_vec)-t_ind_diff
    %     % for each time index, only keep the E-field amplitude value if it is
    %     % equal to the E-field amplitude values at the time-steps corresponding
    %     % to +/-1 wavelength
    %     if abs(E_probe1(ii)) >= 1e-6 && ...
    %         (abs(E_probe1(ii) - E_probe1(ii+t_ind_diff)) < 1e-2 && abs(E_probe1(ii) - E_probe1(ii-t_ind_diff)) < 1e-2)
    % 
    %         E0_f_oper(ii) = E_probe1(ii);
    %     end
    % end
    % % figure; 
    % % plot(E0_f_oper);
    % % hold on;
    % 
    % % slide a window across the filtered time-series data and calculate the max 
    % % value across the window to determine the amplitude of the electric field 
    % % seen by the probe at a given operating frequency 
    % % then, isolate the non-zero stable regions to isolate the incident and
    % % reflected waves and find the magnitude of each
    % window_size = 2*t_ind_diff;  % look at time indices over 2 wavelengths
    % 
    % window_max_E0_f_oper = zeros(size(E_probe1));
    % for ii=1:length(t_vec)-window_size
    %     window_max_E0_f_oper(ii+window_size/2) = max(abs(E0_f_oper(ii+1:ii+window_size)));
    % end
    % plot(t_vec,window_max_E0_f_oper);
    % 
    % nonzeroE_t_ind = find(window_max_E0_f_oper >= 1e-3);
    % disjoint_marker_ind = find(diff(nonzeroE_t_ind) ~= 1);
    % demarcation_t_ind = nonzeroE_t_ind(disjoint_marker_ind);
    % 
    % E0_r = max(window_max_E0_f_oper(demarcation_t_ind:end));
    % E0_i = max(window_max_E0_f_oper(1:demarcation_t_ind));

    [~,refl_wave_sample_t_ind] = min(abs(t_vec - 2.65e-6));
    refl_wave_sample = E_probe1(refl_wave_sample_t_ind:refl_wave_sample_t_ind+t_ind_diff);

    E0_r = max(refl_wave_sample);
    E0_i = 1;

    refl_coef(f_ind) = E0_r/E0_i;
    disp(refl_coef(f_ind));
end

% Plot simulated reflection coefficient as a function of frequency
figure;
plot(f/1e6, refl_coef, 'LineWidth',2);
xlabel("Frequency (MHz)")
ylabel("Input Reflection Coefficient")
title("Input Reflection Coefficient vs. Frequency")

