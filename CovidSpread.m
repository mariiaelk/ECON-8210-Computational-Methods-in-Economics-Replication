% Replication of
% The Puzzling Behavior of Spreads During Covid
% by Stelios Fourakis, Loukas Karabarbounis

%----------------------------------------------------------------
% 0. Housekeeping
%----------------------------------------------------------------

clc
clear
close all

%----------------------------------------------------------------
% 1. Generate results
%----------------------------------------------------------------

params.betaG = 0.96;    % Discount factor of the government
params.sigma = 1.16;       % Risk aversion
params.frisch = 0.56; % Frisch elasticity
params.chiH = 0.46; % Labour level parameter of hand-to-mouth
params.thetaH = 0.36; % Labour efficiency of hand-to-mouth
params.rho = 0.84;       % Autoregressive coefficient
params.sd_shock = 0.04; % Standard deviation of shocks
params.zbar = 0.32; % Mean of technology process
params.r = 0.02;       % Risk-free interest rate
params.mu = 0.05; % Cost of default parameter: level
params.muZ = 0.09; % Cost of default parameter: elasticity
params.lambdaP = 0.1; % Maturity of private debt: good rating
params.kappaP = 0.04; % Coupon of private debt: good rating
params.lambdaD = 0.06; % Maturity of private debt: bad rating
params.kappaD = 0.065; % Coupon of private debt: bad rating
params.lambdaG = 0.04; % Maturity of official debt
params.eta = 0.47;      % Haircut
params.psi = 0.14; % Probability of good credit rating conditional on default
params.issuancecost = 0.8; % Debt issuance cost parameter
params.pi_delta = 0.25; % Probability of decrease in promised official loans
params.pi_delta_hat = 0.74; % Probability of extra official loans in case of default
params.delta = [0,0.35,0.65]; % State space for official loans

% Probability matrices for official lending in case of restructure
params.deltahat_matrix = [1 - params.pi_delta_hat params.pi_delta_hat/2 params.pi_delta_hat/2; ...
                          0 1 0;
                          0 0 1];
params.delta_matrix = [1 0 0; ...
                       params.pi_delta 1-params.pi_delta 0;
                       0 params.pi_delta 1-params.pi_delta];

params.ny = 3;         % Size of the grid for y
params.nB = 15;        % Size of the grid for bonds
params.nF = 5;        % Size of the grid for official loans
params.nNB = 100;     % Size of the grid for revenue from debt issuance
params.alpha_upd = 0.3; % Weight on old guess when updating prices
params.alpha_upd2 = 0; % Weight on old guess when updating VF and probabilities

% Utility function
u = @(c,l) (c.^(1 - params.sigma)-1) / (1 - params.sigma) - ...
    params.chiH*(l.^(1 + 1/params.frisch)) / (1 + 1/params.frisch);

% Obtain value function and policy functions
model = setUpModel(params);
model = vfi(model, u);

%----------------------------------------------------------------
% 2. Make plots
%----------------------------------------------------------------

% Options
% Choose indices for variables which are not varied in plots
choose_f = 2;
choose_delta_hat = 2;
choose_z = 2;

% Choose indices for variables which are varied in plots

% Technology: high versus low
iy_high = 3;
iy_low = 1;

% New official loans: high versus low
ideltahat_high = 3;
ideltahat_low = 1;

% Official loans: high versus low
if_high = 4;
if_low = 2;

% Plot comparison with not in default?
comp = false;

%----------------------------------------------------------------
% 2(a-1). Bond price schedule: technology
%----------------------------------------------------------------

% Define colors for the dashed lines
cmap = colororder();

% Initialize arrays for the plot
x = [];
q_low_nodef = [];
q_high_nodef = [];
q_low_def = [];
q_high_def = [];

% Extract a suitable plot grid
for i = 1:model.nB
    b = model.Bgrid(i);
    if 0 <= b && b <= 20 % Limit values of b_t+1
        x = [x, b]; % Append to x
        q_low_nodef = [q_low_nodef, model.q_nodef(i+(choose_f-1)*model.nB,choose_delta_hat+...
                            (iy_low-1)*size(model.delta,2))]; % Append to q_low
        q_low_def = [q_low_def, model.q_def(i+(choose_f-1)*model.nB,choose_delta_hat+...
                            (iy_low-1)*size(model.delta,2))]; % Append to q_low
        q_high_nodef = [q_high_nodef, model.q_nodef(i+(choose_f-1)*model.nB,choose_delta_hat+...
                            (iy_high-1)*size(model.delta,2))]; % Append to q_high
        q_high_def = [q_high_def, model.q_def(i+(choose_f-1)*model.nB,choose_delta_hat+...
                            (iy_high-1)*size(model.delta,2))]; % Append to q_high
    end
end

% Generate the plot
figure;
plot(x, q_low_nodef, 'DisplayName', 'Low, not in default');
hold on;
plot(x, q_high_nodef, 'DisplayName', 'High, not in default');
if comp 
    plot(x, q_low_def, '--', 'Color', cmap(1,:), 'DisplayName', 'Low, in default');
    plot(x, q_high_def, '--', 'Color', cmap(2,:), 'DisplayName', 'High, in default');
end
title('Bond price schedule $q(b^\prime,f,z,\hat{\delta})$ for different level of technology shock', 'Interpreter', 'latex');
xlabel('$B^\prime$', 'Interpreter', 'latex');
ylabel('$q$', 'Interpreter', 'latex');
legend_handle = legend('show', 'Location', 'northeast'); % Capture the legend handle
legend_title = 'y'; % Legend title text
title(legend_handle, legend_title); % Add title to the legend
hold off;

%----------------------------------------------------------------
% 2(a-2). Bond price schedule: new official loans
%----------------------------------------------------------------

% Define colors for the dashed lines
cmap = colororder();

% Initialize arrays for the plot
x = [];
q_low_nodef = [];
q_high_nodef = [];
q_low_def = [];
q_high_def = [];

% Extract a suitable plot grid
for i = 1:model.nB
    b = model.Bgrid(i);
    if 0 <= b && b <= 20 % Limit values of b_t+1
        x = [x, b]; % Append to x
        q_low_nodef = [q_low_nodef, model.q_nodef(i+(choose_f-1)*model.nB,ideltahat_low+...
                            (choose_z-1)*size(model.delta,2))]; % Append to q_low
        q_low_def = [q_low_def, model.q_def(i+(choose_f-1)*model.nB,ideltahat_low+...
                            (choose_z-1)*size(model.delta,2))]; % Append to q_low
        q_high_nodef = [q_high_nodef, model.q_nodef(i+(choose_f-1)*model.nB,ideltahat_high+...
                            (choose_z-1)*size(model.delta,2))]; % Append to q_high
        q_high_def = [q_high_def, model.q_def(i+(choose_f-1)*model.nB,ideltahat_high+...
                            (choose_z-1)*size(model.delta,2))]; % Append to q_high
    end
end

% Generate the plot
figure;
plot(x, q_low_nodef, 'DisplayName', 'Low, not in default');
hold on;
plot(x, q_high_nodef, 'DisplayName', 'High, not in default');
if comp
    plot(x, q_low_def, '--', 'Color', cmap(1,:), 'DisplayName', 'Low, in default');
    plot(x, q_high_def, '--', 'Color', cmap(2,:), 'DisplayName', 'High, in default');
end
title('Bond price schedule $q(b^\prime,f,z,\hat{\delta})$ for different value of new official loans', 'Interpreter', 'latex');
xlabel('$B^\prime$', 'Interpreter', 'latex');
ylabel('$q$', 'Interpreter', 'latex');
legend_handle = legend('show', 'Location', 'northeast'); % Capture the legend handle
legend_title = 'y'; % Legend title text
title(legend_handle, legend_title); % Add title to the legend
hold off;

%----------------------------------------------------------------
% 2(a-3). Bond price schedule: level of official loans
%----------------------------------------------------------------

% Define colors for the dashed lines
cmap = colororder();

% Initialize arrays for the plot
x = [];
q_low_nodef = [];
q_high_nodef = [];
q_low_def = [];
q_high_def = [];

% Extract a suitable plot grid
for i = 1:model.nB
    b = model.Bgrid(i);
    if 0 <= b && b <= 20 % Limit values of b_t+1
        x = [x, b]; % Append to x
        q_low_nodef = [q_low_nodef, model.q_nodef(i+(if_low-1)*model.nB,choose_delta_hat+...
                            (choose_z-1)*size(model.delta,2))]; % Append to q_low
        q_low_def = [q_low_def, model.q_def(i+(if_low-1)*model.nB,choose_delta_hat+...
                            (choose_z-1)*size(model.delta,2))]; % Append to q_low
        q_high_nodef = [q_high_nodef, model.q_nodef(i+(if_high-1)*model.nB,choose_delta_hat+...
                            (choose_z-1)*size(model.delta,2))]; % Append to q_high
        q_high_def = [q_high_def, model.q_def(i+(if_high-1)*model.nB,choose_delta_hat+...
                            (choose_z-1)*size(model.delta,2))]; % Append to q_high
    end
end

% Generate the plot
figure;
plot(x, q_low_nodef, 'DisplayName', 'Low, not in default');
hold on;
plot(x, q_high_nodef, 'DisplayName', 'High, not in default');
if comp
    plot(x, q_low_def, '--', 'Color', cmap(1,:), 'DisplayName', 'Low, in default');
    plot(x, q_high_def, '--', 'Color', cmap(2,:), 'DisplayName', 'High, in default');
end
title('Bond price schedule $q(b^\prime,f,z,\hat{\delta})$ for different level of accumulated official loans', 'Interpreter', 'latex');
xlabel('$B^\prime$', 'Interpreter', 'latex');
ylabel('$q$', 'Interpreter', 'latex');
legend_handle = legend('show', 'Location', 'northeast'); % Capture the legend handle
legend_title = 'y'; % Legend title text
title(legend_handle, legend_title); % Add title to the legend
hold off;

%----------------------------------------------------------------
% 2(b-1). Value functions: technology
%----------------------------------------------------------------

% Plot value functions for low and high states: technology
ind_start = 1+(choose_f-1)*model.nB;
ind_end = model.nB+(choose_f-1)*model.nB;
figure;
plot(model.Bgrid, model.vf_gc(ind_start:ind_end,choose_delta_hat+...
                            (iy_low-1)*size(model.delta,2)), 'DisplayName', 'Low, good credit status'); % Plot for low state
hold on;
plot(model.Bgrid, model.vf_gc(ind_start:ind_end,choose_delta_hat+...
                            (iy_high-1)*size(model.delta,2)), 'DisplayName', 'High, good credit status'); % Plot for high state
if true
    plot(model.Bgrid, model.vf_bc(ind_start:ind_end,choose_delta_hat+...
                                (iy_low-1)*size(model.delta,2)), '--', 'Color', cmap(1,:), 'DisplayName', 'Low, bad credit status'); % Plot for low state
    plot(model.Bgrid, model.vf_bc(ind_start:ind_end,choose_delta_hat+...
                                (iy_high-1)*size(model.delta,2)), '--', 'Color', cmap(2,:), 'DisplayName', 'High, bad credit status'); % Plot for high state
end

% Add labels, title, and legend
xlabel('$B$', 'Interpreter', 'latex');
ylabel('$V(y,B)$', 'Interpreter', 'latex');
title('Value functions: different levels of technology shock', 'Interpreter', 'latex');
legend_handle = legend('show', 'Location', 'southwest');
title(legend_handle, '$y$', 'Interpreter', 'latex'); % Add title to the legend
hold off;

%----------------------------------------------------------------
% 2(b-2). Value functions: new official loans
%----------------------------------------------------------------

% Plot value functions for low and high states: new official loans
ind_start = 1 + (choose_f - 1) * model.nB;
ind_end = model.nB + (choose_f - 1) * model.nB;

figure;
plot(model.Bgrid, model.vf_gc(ind_start:ind_end, ideltahat_low + ...
    (choose_z - 1) * size(model.delta, 2)), 'DisplayName', 'Low, good credit status'); % Plot for low state
hold on;
plot(model.Bgrid, model.vf_gc(ind_start:ind_end, ideltahat_high + ...
    (choose_z - 1) * size(model.delta, 2)), 'DisplayName', 'High, good credit status'); % Plot for high state
if true
    plot(model.Bgrid, model.vf_bc(ind_start:ind_end, ideltahat_low + ...
        (choose_z - 1) * size(model.delta, 2)), '--', 'Color', cmap(1, :), 'DisplayName', 'Low, bad credit status'); % Plot for low state
    plot(model.Bgrid, model.vf_bc(ind_start:ind_end, ideltahat_high + ...
        (choose_z - 1) * size(model.delta, 2)), '--', 'Color', cmap(2, :), 'DisplayName', 'High, bad credit status'); % Plot for high state
end
% Add labels, title, and legend
xlabel('$B$', 'Interpreter', 'latex');
ylabel('$V(y,B)$', 'Interpreter', 'latex');
title('Value functions: different levels of new official loans', 'Interpreter', 'latex');
legend_handle = legend('show', 'Location', 'southwest');
title(legend_handle, '$y$', 'Interpreter', 'latex'); % Add title to the legend
hold off;

%----------------------------------------------------------------
% 2(b-3). Value functions: accumulated official loans
%----------------------------------------------------------------

% Plot value functions for low and high states: accumulated official loans
ind_start_high = 1+(if_high-1)*model.nB;
ind_end_high = model.nB+(if_high-1)*model.nB;
ind_start_low = 1+(if_low-1)*model.nB;
ind_end_low = model.nB+(if_low-1)*model.nB;

figure;
plot(model.Bgrid, model.vf_gc(ind_start_low:ind_end_low,choose_delta_hat+...
                            (choose_z-1)*size(model.delta,2)), 'DisplayName', 'Low, good credit status'); % Plot for low state
hold on;
plot(model.Bgrid, model.vf_gc(ind_start_high:ind_end_high,choose_delta_hat+...
                            (choose_z-1)*size(model.delta,2)), 'DisplayName', 'High, good credit status'); % Plot for high state
if true
    plot(model.Bgrid, model.vf_bc(ind_start_low:ind_end_low,choose_delta_hat+...
                                (choose_z-1)*size(model.delta,2)), '--', 'Color', cmap(1,:), 'DisplayName', 'Low, bad credit status'); % Plot for low state
    plot(model.Bgrid, model.vf_bc(ind_start_high:ind_end_high,choose_delta_hat+...
                                (choose_z-1)*size(model.delta,2)), '--', 'Color', cmap(2,:), 'DisplayName', 'High, bad credit status'); % Plot for high state
end

% Add labels, title, and legend
xlabel('$B$', 'Interpreter', 'latex');
ylabel('$V(y,B)$', 'Interpreter', 'latex');
title('Value functions: different levels of accumulated official loans', 'Interpreter', 'latex');
legend_handle = legend('show', 'Location', 'southwest');
title(legend_handle, '$y$', 'Interpreter', 'latex'); % Add title to the legend
hold off;

%----------------------------------------------------------------
% 2(b-4). Policy functions: technology
%----------------------------------------------------------------

% Define titles for the plot
title_policy = {'Policy function $B^\prime(b, f, z, \hat{\delta})$: different technology shocks'};

% Grid for B
B_grid = model.Bgrid;

% Policy function as a function of B
x_B = B_grid; % Use B as the x-axis
policy_low_B = zeros(size(B_grid));
policy_high_B = zeros(size(B_grid));

% Compute policy for B'
for i = 1:length(B_grid)

    % Transform indices to actual values using the B grid
    policy_low_B(i) = model.pol_b(i+(choose_f-1)*model.nB,choose_delta_hat+...
                            (iy_low-1)*size(model.delta,2)); % Low state
    policy_high_B(i) = model.pol_b(i+(choose_f-1)*model.nB,choose_delta_hat+...
                            (iy_high-1)*size(model.delta,2)); % High state
end

% Plot for B
figure;
plot(x_B, policy_low_B, 'DisplayName', 'Low');
hold on;
plot(x_B, policy_high_B, 'DisplayName', 'High');
xlabel('$B$', 'Interpreter', 'latex');
ylabel('$B^\prime$', 'Interpreter', 'latex');
title(title_policy{1}, 'Interpreter', 'latex');
legend_handle = legend('show', 'Location', 'northwest');
title(legend_handle, '$y$', 'Interpreter', 'latex'); % Add title to the legend
hold off;

%----------------------------------------------------------------
% 2(b-5). Policy functions: new official loans
%----------------------------------------------------------------

% Define titles for the plot
title_policy = {'Policy function $B^\prime(b, f, z, \hat{\delta})$: different new official loans'};

% Grid for B
B_grid = model.Bgrid;

% Policy function as a function of B
x_B = B_grid; % Use B as the x-axis
policy_low_B = zeros(size(B_grid));
policy_high_B = zeros(size(B_grid));

% Compute policy for B'
for i = 1:length(B_grid)

    % Transform indices to actual values using the B grid
    policy_low_B(i) = model.pol_b(i+(choose_f-1)*model.nB,ideltahat_low+...
                            (choose_z-1)*size(model.delta,2)); % Low state
    policy_high_B(i) = model.pol_b(i+(choose_f-1)*model.nB,ideltahat_high+...
                            (choose_z-1)*size(model.delta,2)); % High state
end

% Plot for B
figure;
plot(x_B, policy_low_B, 'DisplayName', 'Low');
hold on;
plot(x_B, policy_high_B, 'DisplayName', 'High');
xlabel('$B$', 'Interpreter', 'latex');
ylabel('$B^\prime$', 'Interpreter', 'latex');
title(title_policy{1}, 'Interpreter', 'latex');
legend_handle = legend('show', 'Location', 'northwest');
title(legend_handle, '$y$', 'Interpreter', 'latex'); % Add title to the legend
hold off;

%----------------------------------------------------------------
% 2(b-6). Policy functions: accumulated official loans
%----------------------------------------------------------------

% Define titles for the plot
title_policy = {'Policy function $B^\prime(b, f, z, \hat{\delta})$: different level of accumulated official loans'};

% Grid for B
B_grid = model.Bgrid;

% Policy function as a function of B
x_B = B_grid; % Use B as the x-axis
policy_low_B = zeros(size(B_grid));
policy_high_B = zeros(size(B_grid));

% Compute policy for B'
for i = 1:length(B_grid)

    % Transform indices to actual values using the B grid
    policy_low_B(i) = model.pol_b(i+(if_low-1)*model.nB,choose_delta_hat+...
                            (choose_z-1)*size(model.delta,2)); % Low state
    policy_high_B(i) = model.pol_b(i+(if_high-1)*model.nB,choose_delta_hat+...
                            (choose_z-1)*size(model.delta,2)); % High state
end

% Plot for B
figure;
plot(x_B, policy_low_B, 'DisplayName', 'Low');
hold on;
plot(x_B, policy_high_B, 'DisplayName', 'High');
xlabel('$B$', 'Interpreter', 'latex');
ylabel('$B^\prime$', 'Interpreter', 'latex');
title(title_policy{1}, 'Interpreter', 'latex');
legend_handle = legend('show', 'Location', 'northwest');
title(legend_handle, '$y$', 'Interpreter', 'latex'); % Add title to the legend
hold off;


%----------------------------------------------------------------
% 2(c-1). Default probability: technology
%----------------------------------------------------------------

% Bonds and technology, not in default

% Define the grid for the heatmap
ind_start = 1+(choose_f-1)*model.nB;
ind_end = model.nB+(choose_f-1)*model.nB;

% Extract data for the heatmap
Bgrid = model.Bgrid; % Grid for rows
ygrid = model.ygrid; % Grid for columns
restrprob = model.restrprob_nodef(ind_start:ind_end, ...
    choose_delta_hat:size(model.delta, 2):size(model.delta, 2)*model.ny);

% Generate the heatmap
figure;
imagesc(ygrid, Bgrid, restrprob); % Create the heatmap with ygrid on x-axis and Bgrid on y-axis
colorbar; % Add colorbar to show the scale
set(gca, 'YDir', 'normal'); % Flip Y-axis to match standard orientation

% Add labels and title
ylabel('$B^\prime$', 'Interpreter', 'latex');
xlabel('$y$', 'Interpreter', 'latex');
title('Probability of restructuring: not in default', 'Interpreter', 'latex');

% Optionally, customize the color scale
colormap('parula'); % Use the 'hot' colormap for a visually appealing heatmap

% Bonds and technology, in default

% Define the grid for the heatmap
ind_start = 1+(choose_f-1)*model.nB;
ind_end = model.nB+(choose_f-1)*model.nB;

% Extract data for the heatmap
Bgrid = model.Bgrid; % Grid for rows
ygrid = model.ygrid; % Grid for columns
restrprob = model.restrprob_def(ind_start:ind_end, ...
    choose_delta_hat:size(model.delta, 2):size(model.delta, 2)*model.ny);

% Generate the heatmap
figure;
imagesc(ygrid, Bgrid, restrprob); % Create the heatmap with ygrid on x-axis and Bgrid on y-axis
colorbar; % Add colorbar to show the scale
set(gca, 'YDir', 'normal'); % Flip Y-axis to match standard orientation

% Add labels and title
ylabel('$B^\prime$', 'Interpreter', 'latex');
xlabel('$y$', 'Interpreter', 'latex');
title('Probability of restructuring: in default', 'Interpreter', 'latex');

% Optionally, customize the color scale
colormap('parula'); % Use the 'hot' colormap for a visually appealing heatmap

%----------------------------------------------------------------
% 2(c-2). Default probability: new official loans
%----------------------------------------------------------------

% Bonds and new official loans, not in default

% Define the grid for the heatmap
ind_start = 1+(choose_f-1)*model.nB;
ind_end = model.nB+(choose_f-1)*model.nB;

ind_start_col = 1 + (choose_z-1)*model.ny;
ind_end_col = size(model.delta, 2)+(choose_z-1)*model.ny;

% Extract data for the heatmap
Bgrid = model.Bgrid; % Grid for rows
deltahat_grid = model.delta; % Grid for columns
restrprob = model.restrprob_nodef(ind_start:ind_end, ...
    ind_start_col:ind_end_col);

% Generate the heatmap
figure;
imagesc(deltahat_grid, Bgrid, restrprob); % Create the heatmap with ygrid on x-axis and Bgrid on y-axis
colorbar; % Add colorbar to show the scale
set(gca, 'YDir', 'normal'); % Flip Y-axis to match standard orientation

% Add labels and title
ylabel('$B^\prime$', 'Interpreter', 'latex');
xlabel('$\hat{\delta}$', 'Interpreter', 'latex');
title('Probability of restructuring: not in default', 'Interpreter', 'latex');

% Optionally, customize the color scale
colormap('parula'); % Use the 'hot' colormap for a visually appealing heatmap

% Bonds and new official loans, in default

% Define the grid for the heatmap
ind_start = 1+(choose_f-1)*model.nB;
ind_end = model.nB+(choose_f-1)*model.nB;

ind_start_col = 1 + (choose_z-1)*model.ny;
ind_end_col = size(model.delta, 2)+(choose_z-1)*model.ny;

% Extract data for the heatmap
Bgrid = model.Bgrid; % Grid for rows
deltahat_grid = model.delta; % Grid for columns
restrprob = model.restrprob_def(ind_start:ind_end, ...
    ind_start_col:ind_end_col);

% Generate the heatmap
figure;
imagesc(deltahat_grid, Bgrid, restrprob); % Create the heatmap with ygrid on x-axis and Bgrid on y-axis
colorbar; % Add colorbar to show the scale
set(gca, 'YDir', 'normal'); % Flip Y-axis to match standard orientation

% Add labels and title
ylabel('$B^\prime$', 'Interpreter', 'latex');
xlabel('$\hat{\delta}$', 'Interpreter', 'latex');
title('Probability of restructuring: in default', 'Interpreter', 'latex');

% Optionally, customize the color scale
colormap('parula'); % Use the 'hot' colormap for a visually appealing heatmap

%----------------------------------------------------------------
% 2(c-3). Default probability: accumulated official loans
%----------------------------------------------------------------

% Bonds and official loans, not in default
% Extract data for the heatmap
Bgrid = model.Bgrid; % Grid for rows
Fgrid = model.Fgrid; % Grid for columns
restrprob = model.restrprob_nodef(:, ...
    choose_delta_hat+(choose_z-1)*size(model.delta,2));
restrprob = reshape(restrprob,model.nB,model.nF);

% Generate the heatmap
figure;
imagesc(Fgrid, Bgrid, restrprob); % Create the heatmap with ygrid on x-axis and Bgrid on y-axis
colorbar; % Add colorbar to show the scale
set(gca, 'YDir', 'normal'); % Flip Y-axis to match standard orientation

% Add labels and title
ylabel('$B^\prime$', 'Interpreter', 'latex');
xlabel('$F^\prime$', 'Interpreter', 'latex');
title('Probability of restructuring: not in default', 'Interpreter', 'latex');

% Optionally, customize the color scale
colormap('parula'); % Use the 'hot' colormap for a visually appealing heatmap

% Bonds and official loans, in default
% Extract data for the heatmap
Bgrid = model.Bgrid; % Grid for rows
Fgrid = model.Fgrid; % Grid for columns
restrprob = model.restrprob_def(:, ...
    choose_delta_hat+(choose_z-1)*size(model.delta,2));
restrprob = reshape(restrprob,model.nB,model.nF);

% Generate the heatmap
figure;
imagesc(Fgrid, Bgrid, restrprob); % Create the heatmap with ygrid on x-axis and Bgrid on y-axis
colorbar; % Add colorbar to show the scale
set(gca, 'YDir', 'normal'); % Flip Y-axis to match standard orientation

% Add labels and title
ylabel('$B^\prime$', 'Interpreter', 'latex');
xlabel('$F^\prime$', 'Interpreter', 'latex');
title('Probability of restructuring: not in default', 'Interpreter', 'latex');

% Optionally, customize the color scale
colormap('parula'); % Use the 'hot' colormap for a visually appealing heatmap

%----------------------------------------------------------------
% 3. Functions
%----------------------------------------------------------------

function resid = hand_to_mouth(c, transfers, technology, params)

    l = ((technology*params.thetaH/params.chiH)*c^(-params.sigma))^params.frisch;
    resid = c - technology*params.thetaH*l - transfers;

end

function model = setUpModel(params)

    options = optimset('Display','off','TolFun',10^(-8),'TolX',10^(-8));

    % Grids
    % Private bonds grid
    %params.Bgrid = linspace(0, 20, params.nB);  % Private bonds grid
    fine_fraction = 0.5; % Fraction of points between 0 and 2
    % Generate a linearly spaced grid from 0 to 1
    lin_grid = linspace(0, 1, params.nB + 1);

    % Apply a nonlinear transformation to distribute points
    transformed_grid = lin_grid.^2; % Use a power function to create a fine grid near 0

    % Scale the grid so that 50% of points lie between 0 and 2
    fine_points = floor((params.nB + 1) * fine_fraction);
    coarse_points = params.nB + 1 - fine_points;

    % Define the breakpoints and scaling
    fine_scale = linspace(0, 2, fine_points);
    coarse_scale = linspace(2, 20, coarse_points);

    % Combine the two parts of the grid
    params.Bgrid = [fine_scale, coarse_scale(2:end)]; % Avoid duplicate point at 2

    % Other grids
    params.Fgrid = linspace(0, 40, params.nF);     % Official loans grid
    params.NBgrid = linspace(-10, 15, params.nNB);     % Revenue from debt issuance grid
    [params.trans_mat, grid_logy] = tauchen(params.ny, params.rho, ...
        params.sd_shock, params.zbar); % Transition matrix and states for income
    params.trans_mat_transpose = params.trans_mat';
    params.ygrid = exp(grid_logy); % Income grid
    [params.BGgrid_F, params.BGgrid_B] = meshgrid(params.Fgrid,params.Bgrid);

    % Create a grid with rows (b,f) and columns (z,delta)
    % Matrix with values of f
    params.Fgrid_matrix = []; % Initialize the final matrix
    for i = 1:params.nF
        % Create a params.nF by params.ny*size(params.delta,2) matrix filled with the current value of kk_grid
        temp_matrix = params.Fgrid(i) * ones(params.nB, params.ny*size(params.delta,2));
        
        % Stack this matrix to the final matrix
        params.Fgrid_matrix = [params.Fgrid_matrix; temp_matrix];
    end
    % Matrix with values of b
    params.Bgrid_matrix = repmat(params.Bgrid', params.nF, params.ny*size(params.delta,2));
    % Matrix with values of z
    params.ygrid_matrix = []; % Initialize the final matrix
    for i = 1:params.ny
        % Create a params.nF*params.nB by size(params.delta,2) matrix filled with the current value of kk_grid
        temp_matrix = params.ygrid(i) * ones(params.nF*params.nB, size(params.delta,2));
        
        % Stack this matrix to the final matrix
        params.ygrid_matrix = [params.ygrid_matrix temp_matrix];
    end
    % Matrix with values of delta
    params.deltagrid_matrix = repmat(params.delta, params.nB*params.nF, params.ny);

    % Find optimal consumption and labor conditional on NB (rows) and z
    % (columns)
    params.consH = ones(params.nNB,params.ny)*(-999);
    params.laborH = zeros(params.nNB,params.ny);
    [params.NBgrid_j, params.ygrid_j] = meshgrid(params.NBgrid, params.ygrid);
    params.NBgrid_j = params.NBgrid_j';
    params.ygrid_j = params.ygrid_j';

    for i = 1:params.nNB
        for j = 1:params.ny
            params.consH(i,j) = fsolve(@(x) hand_to_mouth(x, params.NBgrid(i), ...
                params.ygrid(j), params), 1, options);
        end
    end
    params.laborH = ((params.ygrid_j*params.thetaH/params.chiH).*...
        params.consH.^(-params.sigma)).^params.frisch;

    % Define value functions
    params.evf_d = zeros(size(params.deltagrid_matrix)); % Expected value function in default
    params.evf_nd = zeros(size(params.deltagrid_matrix)); % Expected value function conditional on no default
    params.vf_d_nr = zeros(size(params.deltagrid_matrix)); % Value function conditional on default and no restructuring: choose level of debt
    params.vf_d_r = zeros(size(params.deltagrid_matrix)); % Value function conditional on default and restructuring: choose level of debt
    params.vf_n = zeros(size(params.deltagrid_matrix)); % Value function conditional on no default and no restructuring: choose level of debt
    params.vf_gc = zeros(size(params.deltagrid_matrix)); % Value function conditional on good credit standing: choose whether to restructure
    params.vf_bc = zeros(size(params.deltagrid_matrix)); % Value function conditional on bad credit standing: choose whether to restructure
    
    params.pol_eta_gc = zeros(size(params.deltagrid_matrix)); % Policy function: whether to default conditional on good credit standing
    params.pol_eta_bc = zeros(size(params.deltagrid_matrix)); % Policy function: whether to default conditional on bad credit standing
    params.pol_b = zeros(size(params.deltagrid_matrix)); % Policy function: for debt conditional on good credit standing and not restructuring
    params.pol_labor = zeros(size(params.deltagrid_matrix)); % Policy function: for labor conditional on good credit standing and not restructuring
    params.pol_c = zeros(size(params.deltagrid_matrix)); % Policy function: for consumption conditional on good credit standing and not restructuring
    
    params.q_nodef = ones(size(params.deltagrid_matrix)) .* (1 / (1 + params.r)); % Pricing function conditional on no default
    params.q_def = ones(size(params.deltagrid_matrix)) .* (1 / (1 + params.r)); % Pricing function conditional on default
    params.restrprob_nodef = zeros(size(params.deltagrid_matrix)); % Restructure probabilities conditional on no default
    params.restrprob_def = zeros(size(params.deltagrid_matrix)); % Restructure probabilities conditional on default
    
    % Store results
    model = params;

end

% This code updates all value functions and policy functions given
% price schedule and restructure probabilities

function model = one_step_update(model, u)
    
    % Loop over possible realizations of income
    for ind_y = 1:model.ny
        
        % Loop over possible realizations of new official loans
        for ind_delta = 1:size(model.delta,2)

                % Loop over possible values of official loans
                for ind_f = 1:model.nF

                    % Loop over possible values of debt
                    for ind_b = 1:model.nB
        
                        y = model.ygrid(ind_y);
                        f = model.Fgrid(ind_f);
                        b = model.Bgrid(ind_b);
                        cur_delta_hat = model.delta(ind_delta);

                        % Value in bad credit standing or current
                        % restructuring
                        nb_val = cur_delta_hat - (model.lambdaD + ...
                            model.kappaD)*b - (model.lambdaG + ...
                            model.r)*f;
                        fp = (1-model.lambdaG)*f + cur_delta_hat;

                        % If fp lies outside of the grid, adjust it 
                        if fp > model.Fgrid(end)
                            fp = model.Fgrid(end);
                            cur_delta_hat = fp - (1-model.lambdaG)*f;
                            nb_val = cur_delta_hat - (model.lambdaD + ...
                                model.kappaD)*b - (model.lambdaG + ...
                                model.r)*f;
                        end

                        % Current utility
                        c = interp1(model.NBgrid_j(:,ind_y),model.consH(:,ind_y),...
                            nb_val,'linear','extrap');
                        l = interp1(model.NBgrid_j(:,ind_y),model.laborH(:,ind_y),...
                            nb_val,'linear','extrap');
                        cur_util = u(c,l);

                        % New debt in case of bad credit standing, no restructuring
                        bp_bc_nr = (1-model.lambdaD)*b;

                        % New debt in case of bad credit standing, restructuring
                        bp_bc_r = (1-model.lambdaD)*(1-model.eta)*b;

                        % Value of bad credit standing, no restructure
                        mmatrix = model.evf_d(:,ind_delta+...
                            (ind_y-1)*size(model.delta,2));
                        % Reshape the data
                        mmatrix = reshape(mmatrix, model.nB, model.nF);
                        exp_util = interp2(model.BGgrid_F, model.BGgrid_B, ...
                            mmatrix, fp, bp_bc_nr, 'spline');
                        model.vf_d_nr(ind_b+(ind_f-1)*model.nB,ind_delta+...
                            (ind_y-1)*size(model.delta,2)) = cur_util + ...
                            model.betaG*exp_util;

                        % Value of bad credit standing, restructure
                        exp_util = interp2(model.BGgrid_F, model.BGgrid_B, ...
                            mmatrix, fp, bp_bc_r, 'spline');
                        model.vf_d_r(ind_b+(ind_f-1)*model.nB,ind_delta+...
                            (ind_y-1)*size(model.delta,2)) = cur_util - ...
                            - model.mu - model.muZ * log(model.ygrid(ind_y)) ...
                            + model.betaG*exp_util;

                        % Value of good credit standing
                        % Loop over possible next period values of debt
                        current_max = -1e14;
                        pol_ind = 0;
                        for ind_bp = 1:model.nB
                
                            bp = model.Bgrid(ind_bp);
                            % NB if a certain level of debt is chosen
                            price = model.q_nodef(ind_bp+(ind_f-1)*model.nB,...
                                ind_delta+(ind_y-1)*size(model.delta,2));
                            %mmatrix = model.q_nodef(:,ind_delta+...
                            %    (ind_y-1)*size(model.delta,2));
                            % Reshape the data
                            %mmatrix = reshape(mmatrix, model.nB, model.nF);
                            %price = interp2(model.BGgrid_F, model.BGgrid_B, ...
                            %    mmatrix, f, bp, 'spline');
                            restr_def = interp1(model.Fgrid,...
                                model.restrprob_nodef(find(model.Bgrid_matrix(:,1) == ...
                                model.Bgrid(ind_bp)),ind_delta+...
                                (ind_y-1)*size(model.delta,2)),...
                                fp,'linear','extrap');
                            %mmatrix = model.restrprob_nodef(:,ind_delta+...
                            %    (ind_y-1)*size(model.delta,2));
                            % Reshape the data
                            %mmatrix = reshape(mmatrix, model.nB, model.nF);
                            %prob_def = interp2(model.BGgrid_F, model.BGgrid_B, ...
                            %    mmatrix, fp, bp, 'spline');
                            issuance_cost = price*(bp - (1-model.lambdaP)*b)*...
                                0.5*(1+sin(((restr_def-model.issuancecost)/...
                                (1-model.issuancecost) - 0.5)*pi))*...
                                (bp - (1-model.lambdaP)*b > 0);
                            nb_val = price*(bp - (1-model.lambdaP)*b) + ...
                                    cur_delta_hat - (model.lambdaP + ...
                                    model.kappaP)*b - (model.lambdaG + ...
                                    model.r)*f - issuance_cost; 

                            % Current utility
                            c = interp1(model.NBgrid_j(:,ind_y),model.consH(:,ind_y),...
                                nb_val,'linear','extrap');
                            l = interp1(model.NBgrid_j(:,ind_y),model.laborH(:,ind_y),...
                                nb_val,'linear','extrap');
                            cur_util = u(c,l);
                            % Expected utility
                            exp_util = interp1(model.Fgrid,...
                                model.evf_nd(find(model.Bgrid_matrix(:,1) == ...
                                model.Bgrid(ind_bp)),ind_delta+...
                                (ind_y-1)*size(model.delta,2)),...
                                fp,'linear','extrap');
                            %mmatrix = model.evf_nd(:,ind_delta+...
                            %    (ind_y-1)*size(model.delta,2));
                            % Reshape the data
                            %mmatrix = reshape(mmatrix, model.nB, model.nF);
                            %exp_util = interp2(model.BGgrid_F, model.BGgrid_B, ...
                            %mmatrix, fp, bp, 'spline');
                            % Value if this level of debt is chosen
                            v_val = cur_util + model.betaG*exp_util;
                            if v_val > current_max
                                current_max = v_val;
                                pol_ind = ind_bp;
                                lopt = l;
                                copt = c;
                            end
    
                        end

                        % Update value and policy functions
                         model.vf_n(ind_b+(ind_f-1)*model.nB,ind_delta+...
                                (ind_y-1)*size(model.delta,2)) = current_max;
                         model.pol_b(ind_b+(ind_f-1)*model.nB,ind_delta+...
                                (ind_y-1)*size(model.delta,2)) = model.Bgrid(pol_ind);
                         model.pol_labor(ind_b+(ind_f-1)*model.nB,ind_delta+...
                                (ind_y-1)*size(model.delta,2)) = lopt;
                         model.pol_c(ind_b+(ind_f-1)*model.nB,ind_delta+...
                                (ind_y-1)*size(model.delta,2)) = copt;

                         % Now determine whether to restructure
                         % Now we have delta, not delta_hat as state
                         % If government begins period in bad credit
                         % standing
                         % Expected value if the government restructures
                         start_ind = 1+(ind_y-1)*size(model.delta,2);
                         val_d_r_bar = model.vf_d_r(ind_b+(ind_f-1)*model.nB,...
                            start_ind:(start_ind+size(model.delta,2)-1))*...
                            model.deltahat_matrix(ind_delta,:)';
                         % Value if the government does not restructure
                         val_d_nr = model.vf_d_nr(ind_b+(ind_f-1)*model.nB,ind_delta+...
                                (ind_y-1)*size(model.delta,2));
                         model.vf_bc(ind_b+(ind_f-1)*model.nB,ind_delta+...
                                (ind_y-1)*size(model.delta,2)) = max(val_d_r_bar,val_d_nr);
                         % If equal, suppose we choose not to restructure
                         model.pol_eta_bc(ind_b+(ind_f-1)*model.nB,ind_delta+...
                                (ind_y-1)*size(model.delta,2)) = val_d_r_bar > val_d_nr;

                         % If government begins period in good credit
                         % standing
                         % Value if the government does not restructure
                         val_n = model.vf_n(ind_b+(ind_f-1)*model.nB,ind_delta+...
                                (ind_y-1)*size(model.delta,2));
                         model.vf_gc(ind_b+(ind_f-1)*model.nB,ind_delta+...
                                (ind_y-1)*size(model.delta,2)) = max(val_d_r_bar,val_n);
                         % If equal, suppose we choose not to restructure
                         model.pol_eta_gc(ind_b+(ind_f-1)*model.nB,ind_delta+...
                                (ind_y-1)*size(model.delta,2)) = val_d_r_bar > val_n;
            
                    end
        
                end
        
        end

    end

    % Update expected value function
    % Loop over possible realizations of income: y
    for ind_y = 1:model.ny
        
        % Loop over possible realizations of new official loans: hat delta
        for ind_delta = 1:size(model.delta,2)

            % Update expected value function: no prior default
            prob_val = repmat(model.trans_mat(ind_y,:),model.ny,1);
            prob_val = prob_val(:)';
            prob_val = prob_val .* repmat(model.delta_matrix(ind_delta,:),...
                1,size(model.delta,2));
            model.evf_nd(:,ind_delta+(ind_y-1)*size(model.delta,2)) = ...
                model.vf_gc*prob_val';
    
            % Update expected value function: prior default
            model.evf_d(:,ind_delta+(ind_y-1)*size(model.delta,2)) = ...
                model.psi*model.vf_gc*prob_val' + ...
                (1-model.psi)*model.vf_bc*prob_val';

            % Update restructure probabilities: no prior default
            model.restrprob_nodef(:,ind_delta+(ind_y-1)*size(model.delta,2)) = ...
                model.pol_eta_gc*prob_val';
    
            % Update restructure probabilities: prior default
            model.restrprob_def(:,ind_delta+(ind_y-1)*size(model.delta,2)) = ...
                model.psi*model.pol_eta_gc*prob_val' + ...
                (1-model.psi)*model.pol_eta_bc*prob_val';

        end

    end

end       

% This function updates prices and restructure probabilities given
% value functions and policy functions

function model = compute_prices(model)

    % Loop over possible realizations of income: y
    for ind_y = 1:model.ny
        
        % Loop over possible realizations of new official loans: delta hat
        for ind_delta = 1:size(model.delta,2)

                % Update restructure probabilities
                prob_val = repmat(model.trans_mat(ind_y,:),model.ny,1);
                prob_val = prob_val(:)';
                prob_val = prob_val .* repmat(model.delta_matrix(ind_delta,:),...
                    1,size(model.delta,2));

                % Loop over possible values of official loans: f
                for ind_f = 1:model.nF

                    % Loop over possible values of debt: bp
                    for ind_b = 1:model.nB

                    y = model.ygrid(ind_y);
                    f = model.Fgrid(ind_f);
                    bp = model.Bgrid(ind_b);
                    cur_delta_hat = model.delta(ind_delta);
                    % f' can be off the grid
                    fp = (1-model.lambdaG)*f + model.delta(ind_delta);

                    if fp > model.Fgrid(end)
                        fp = model.Fgrid(end);
                        cur_delta_hat = fp - (1-model.lambdaG)*f;
                    end

                    % Optimal b'': conditional on b',f',z',delta hat' and good credit
                    % standing next period
                    b_indices = find(model.Bgrid_matrix(:,1) == model.Bgrid(ind_b));
                    mmatrix = model.pol_b(b_indices,:);
                    bpp = interp1(model.Fgrid',mmatrix,...
                            fp,'linear','extrap');
                    bpp = min(bpp, model.Bgrid(end));
                    bpp = max(bpp,model.Bgrid(1));
                    mmatrix = model.pol_eta_gc(b_indices,:);
                    % Optimal eta' conditional on b',f',z',delta' and good
                    % credit standing next period
                    etapp_gc = interp1(model.Fgrid',mmatrix,...
                            fp,'linear','extrap');
                    mmatrix = model.pol_eta_bc(b_indices,:);
                    etapp_bc = interp1(model.Fgrid',mmatrix,...
                            fp,'linear','extrap');

                    % bpp is possibly off the grid

                    % Determine q(b'', f',z',delta hat',default' = 0) 
                    % and q(b'', f',z',delta hat',default' = 1)
                    q_bpp_nodef = zeros(1,size(model.q_nodef,2));
                    q_bpp_def = zeros(1,size(model.q_nodef,2));
                    for jj = 1:size(model.q_nodef,2)
                        mmatrix = model.q_nodef(:,jj);
                        mmatrix = reshape(mmatrix, model.nB, model.nF);
                        q_bpp_nodef(jj) = interp2(model.BGgrid_F, model.BGgrid_B, ...
                                    mmatrix, fp, bpp(jj), 'spline');

                        mmatrix = model.q_def(:,jj);
                        mmatrix = reshape(mmatrix, model.nB, model.nF);
                        q_bpp_def(jj) = interp2(model.BGgrid_F, model.BGgrid_B, ...
                                    mmatrix, fp, bpp(jj), 'spline');
                    end

                    % Update price for the case when default = 0
                    % Part when gov't chooses not to restructure
                    p1 = (1-etapp_gc).*(model.lambdaP+model.kappaP+(1-model.lambdaP)*...
                        q_bpp_nodef); 
                    p1 = p1*prob_val';
                    % Part when gov't chooses to restructure
                    % In this case we need to restructure delta' and
                    % delta_hat'
                    % Probabilities for 9 combinations of (z,delta') +
                    % delta_hat' = 0; then 9 combinations with delta_hat'=
                    % deltaL; then 9 combinations with delta_hat'=deltaH
                    prob_val_full = repmat(prob_val,1,size(model.delta,2));
                    etapp_gc_full = repmat(etapp_gc,1,size(model.delta,2));
                    prob_val_full_hat = [];
                    for jj = 1:size(model.delta,2)
                        prob_val_full_hat = [prob_val_full_hat; ...
                            repmat(model.deltahat_matrix(:,jj),size(model.delta,2),1)];
                    end
                    prob_val_full = prob_val_full.*prob_val_full_hat';
                    q_bpp_def_full = reshape(q_bpp_def,model.ny,size(model.delta,2))';
                    q_bpp_def_full = q_bpp_def_full(:);
                    q_bpp_def_full = repelem(q_bpp_def_full,3);

                    p2 = etapp_gc_full.*(model.lambdaD+model.kappaD+...
                        (1-model.eta)*(1-model.lambdaD)*q_bpp_def_full');
                    p2 = p2*prob_val_full';
                    model.q_nodef(ind_b+(ind_f-1)*model.nB,ind_delta+...
                                (ind_y-1)*size(model.delta,2)) = ...
                                (p1+p2)/(1+model.r);

                    % Update price for the case when default = 1
                    etapp_bc_full = repmat(etapp_bc,1,size(model.delta,2));
                    p3 = model.lambdaD+model.kappaD+(1-etapp_bc_full*model.eta)*...
                        (1-model.lambdaD).*q_bpp_def_full';
                    p3 = p3*prob_val_full';
                    model.q_def(ind_b+(ind_f-1)*model.nB,ind_delta+...
                                (ind_y-1)*size(model.delta,2)) = ...
                                (model.psi*(p1+p2)+(1-model.psi)*p3)/(1+model.r);

                    end
                end
        end
    end    

end

function model = vfi(model, u, tol, maxit)

    % Default values
    if nargin < 3
        tol2 = 0.01; % For VF
        tol1 = 0.03; % For prices
    end

    if nargin < 4
        maxit = 10000;
    end

    % Iterations for prices
    it_price = 0;
    dist_price = 100;

    while dist_price > tol1 && it_price < maxit

        it_price = it_price + 1;
        q_nodef_enter = model.q_nodef;
        q_def_enter = model.q_def;

        % Iteration stuff - VF
        it = 0;
        dist = 100;
    
        % Allocate memory for update
        %V_enter = zeros(size(ae.vf), class(ae.vf));
    
        while dist > tol2 && it < maxit
    
            it = it + 1;
    
            % Compute expectations for this iterations
            % (we need Pi' because of order value function dimensions)
            EVd_enter = model.evf_d;
            EVnd_enter = model.evf_nd;
            restrprob_def_enter = model.restrprob_def;
            restrprob_nodef_enter = model.restrprob_nodef;
    
            % Update value function
            model = one_step_update(model, u);
            
            % Check convergence of the value function
            dist = max([max(abs(EVd_enter - model.evf_d), [], 'all'),...
                       max(abs(EVnd_enter - model.evf_nd), [], 'all')]);
                       %max(abs(restrprob_def_enter - model.restrprob_def), [], 'all'),...
                       %max(abs(restrprob_nodef_enter - model.restrprob_nodef), [], 'all')]);


            % Update value functions and probabilities if no convergence
            if dist > tol2
                model.evf_d = (1-model.alpha_upd2)*model.evf_d + model.alpha_upd2*EVd_enter;
                model.evf_nd = (1-model.alpha_upd2)*model.evf_nd + model.alpha_upd2*EVnd_enter;
                model.restrprob_def = (1-model.alpha_upd2)*model.restrprob_def + model.alpha_upd2*restrprob_def_enter;
                model.restrprob_nodef = (1-model.alpha_upd2)*model.restrprob_nodef + model.alpha_upd2*restrprob_nodef_enter;
            end
    
            if mod(it, 25) == 0
                fprintf('Finished iteration %d for given prices with dist of %f\n', it, dist);
            end
    
        end

        % Update prices
        model = compute_prices(model);

        % Check convergence of prices
        dist_price = max([max(abs(q_nodef_enter - model.q_nodef), [], 'all'),...
                         max(abs(q_def_enter - model.q_def), [], 'all')]);
                         %max(abs(restrprob_def_enter - model.restrprob_def), [], 'all'),...
                         %max(abs(restrprob_nodef_enter - model.restrprob_nodef), [], 'all')]);

        if mod(it, 1) == 0
                fprintf('Price interation %d with dist of %f\n', it_price, dist_price);
        end

        % Update price if no convergence
        if dist_price > tol1
            model.q_nodef = (1-model.alpha_upd)*model.q_nodef + model.alpha_upd*q_nodef_enter;
            model.q_def = (1-model.alpha_upd)*model.q_def + model.alpha_upd*q_def_enter;
            %model.restrprob_def = (1-model.alpha_upd)*model.restrprob_def + model.alpha_upd*restrprob_def_enter;
            %model.restrprob_nodef = (1-model.alpha_upd)*model.restrprob_nodef + model.alpha_upd*restrprob_nodef_enter;
        end

    end

end

