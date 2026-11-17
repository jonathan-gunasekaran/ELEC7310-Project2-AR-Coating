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
delta_x = 0.01;         % spatial step in meters
max_x = 100;            % max bound on simulation space
x = 0:delta_x:max_x;    % spatial array

% define temporal array based on spatial step value to ensure numerical
% stability
tmax = 100;
S = 0.8;        % Courant stability factor (want less than 1 in all media)


% define multi-layer stackup material parameters
layer_delx = [];
layer_eps_r = [];

eps = ones(size(x)) .* epsilon_0; 

for ii = 1:length(layer_delx)
    layer_x_ind = layer_delx(ii,:);
    eps(layer_x_ind) = eps(layer_x_ind) .* layer_eps_r(ii);
end


% initialize electric field vector - store values for past, present, and
% future time step
E = zeros(length(x), 3);


%% 1-D FDTD Loop





%% Plotting


