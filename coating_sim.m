%% ELEC 7310 - Advanced Electrodynamics I - Project 2
% Jonathan Gunasekaran
% 2026/11/15


%% Simulation Setup
% constants
epsilon_0 = 8.854e-12;
mu_0 = pi*4e-7;
c = 1/sqrt(mu_0*epsilon_0);

% set simulation frequency - eventually setup as input for frequency sweep
% analysis
f = 20e6;
lambda_0 = c/f;


% define 1-D spatial array 
delta_x = 0.05;         % spatial step in meters
max_x = 100;            % max bound on simulation space
x = 0:delta_x:max_x;    % spatial array

% define temporal array based on spatial step value to ensure numerical
% stability
tmax = 100e-6;             % simulation end time
S = 0.8;                % Courant stability factor (want less than 1 in all media)
delta_t = S*delta_x/c;  % time step
t = 0;                  % current simulation time (start at 0)

% define multi-layer stackup material parameters
num_layers = 2;
layer_eps_r = [2, 3];
eps = ones(size(x)) .* epsilon_0; 

layer_start_x = 40;     % start layers at 40 m

for ii=1:num_layers
    layer_end_x = round(layer_start_x+lambda_0/4, 2);
    
    layer_x_ind = find(x >= layer_start_x & x < layer_end_x);
    eps(layer_x_ind) = eps(layer_x_ind) .* layer_eps_r(ii);
    
    layer_start_x = layer_end_x;
end


% initialize electric field vector - store values for past, present, and
% future time step
E = zeros(length(x), 3);


% define locations for the virtual electric field probes 
boundary_x_ind = find(diff(eps));
ind_layer_width = unique(diff(boundary_x_ind));
probe1_x_ind = boundary_x_ind(1) - ind_layer_width;
probe2_x_ind = boundary_x_ind(end) + ind_layer_width;


%% 1-D FDTD Loop
tic
while t <= tmax
    % disp(t)
    
    E(2:end-1,3) = (delta_t/delta_x)^2 ./ (mu_0.*eps(2:end-1).') ...
                    .* (E(3:end,2) - 2*E(2:end-1,2) + E(1:end-2,2)) ... 
                    + 2*E(2:end-1,2) ...
                    - E(2:end-1, 1);
    
    E(2:end-1,1) = E(2:end-1,2);
    E(2:end-1,2) = E(2:end-1,3);
    t = t + delta_t;

end
toc



%% Plotting


