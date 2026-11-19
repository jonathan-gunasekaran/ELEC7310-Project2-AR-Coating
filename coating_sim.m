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
delta_x = 0.01;         % spatial step in meters
max_x = 500;            % max bound on simulation space
x = 0:delta_x:max_x;    % spatial array

% define temporal array based on spatial step value to ensure numerical
% stability
tmax = 1.2e-6;             % simulation end time
S = 0.8;                % Courant stability factor (want less than 1 in all media)
delta_t = S*delta_x/c;  % time step
t = 0;                  % current simulation time (start at 0)

% define multi-layer stackup material parameters
load_eps_r = 5;
layer_eps_r = [2, 3];
num_layers = length(layer_eps_r);

% relative permittivity spatial array
eps = ones(size(x)) .* epsilon_0; 

layer_start_x = 0.4*max_x;  % start point for multi-layer interface in x

for ii=1:num_layers
    layer_end_x = round(layer_start_x+lambda_0/4, 2);
    
    layer_x_ind = find(x >= layer_start_x & x < layer_end_x);
    eps(layer_x_ind) = eps(layer_x_ind) .* layer_eps_r(ii);
    
    layer_start_x = layer_end_x;
end

load_x_ind = find(x >= layer_start_x);
eps(load_x_ind) = eps(load_x_ind) .* load_eps_r;


% initialize electric field vector - store values for past, present, and
% future time step
E = zeros(length(x), 3);


% define locations for the virtual electric field probes 
boundary_x_ind = find(diff(eps));

x_dist_to_probe = max_x/10;
probe1_x_ind = boundary_x_ind(1) - x_dist_to_probe/delta_x;
probe2_x_ind = boundary_x_ind(1) + x_dist_to_probe/delta_x;


%% 1-D FDTD Loop
fdtd_movie = figure;
frame_decimation_rate = 500;

tic
loop_iter = 0;
while t <= tmax
    loop_iter = loop_iter + 1;
    
    % Add source (incoming E-field) - simulate for 5 cycles
    if t <= 5/f0
        E(1,3) = sin(2*pi*f0*t);
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

    % increment time
    t = t + delta_t;


    % 1-D FDTD sim plot movie
    if mod(loop_iter, frame_decimation_rate) == 0
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



%% Plotting

