clear; close all;

%% 0. Initialize Parameters
L = 1200;               % Length of bridge
n = L;               % Discretize into 1 mm seg.
P = 400;                % Total weight of train [N]
x = linspace(0, L, n); % x-axis

% Note: the position is 1mm greater then the physical location from the left
% edge of the bridge; A_pos = 1 means the location 0mm from the left of the
% bridge.
A_pos = 1;            % position of first support from left edge
B_pos = A_pos + 1200; % position of second support from left edge


%% 1. SFD, BMD under train loading
x_train = [52 228 392 568 732 908]; % Train Load Locations
P_train = [1 1 1 1 1 1] * P/6;

% distance_from_train_left_axel_to_right_axel
dis_between_l_r_axels = max(x_train) - min(x_train)

n_train = 900;                 % num of train locations, as in how many different simulations of train locations we're simulating
n_positions = L + 1 % We want information from pos=1 (0mm) to pos=1201 (1200mm) 

% Create an array of positions where the right axel of the train can be
% located, from pos=1 (0mm from left of bridge) to pos=1201 + distance
% between wheels. This will make the left axel be at pos=1201 (1200mm) at
% the last position.
train_space = floor(linspace(1, n_positions + dis_between_l_r_axels, n_train));
% For the following 3 matrices, the information for each simulation is
% stored in the row and there can be information at each 1mm segment across
% the bridge. 
loadsi = zeros(n_train, n_positions);    % 1 load plot for each train loc.
SFDi = zeros(n_train, n_positions);     % 1 SFD for each train loc.
BMDi = zeros(n_train, n_positions);     % 1 BMD for each train loc.

% Solve for SFD and BMD with the train at different locations
for J = 1:n_train
    i = train_space(J); % how much we're incrementing the position of the train, the position of the right axel of the train
    % start location of train
    ttrain_pos_i = x_train - max(x_train) + i; % train position 
    ttrain_pos_i(ttrain_pos_i <= 0) = NaN; % Get positions that are on the bridge only (take out positions before 0 position of the bridge)
    % any train position greater than 1201 means it's off the bridge since
    % pos 1201 = 1200mm from the left edge
    ttrain_pos_i(ttrain_pos_i > n_positions) = NaN; % Get positions that are on the bridge only (take out positions after 1200 position of the bridge)

    p_ts = P_train(~isnan(ttrain_pos_i)); % get loads at the positions that are on the bridge only
    pos_ti = ttrain_pos_i(~isnan(ttrain_pos_i)); % get positions of the train that are on the bridge only

    % sum of moments at A eqn
    rxn_B = sum((A_pos - pos_ti) .* p_ts) / -L;

    % sum of Fy eqn
    rxn_A = sum(p_ts) - rxn_B;
    % construct applied loads
    % w(x)
    loads = zeros(1, n_positions);
    % set forces from the train (negative loads)
    loads(pos_ti) = -1 .* p_ts;
    % set reaction forces at their positions with their respective forces,
    % you may need to add these forces because if the wheel is at the
    % support there is no net load.
    loads(A_pos) = loads(A_pos) + rxn_A;
    loads(B_pos) = loads(B_pos) + rxn_B;


    %loadi stuff
    loadsi(i,:) = loads; % place loads at each position [0, 0, F_A, . . . , -P/6. . ., 0 . . . ]
    % SFD = num. integral(w)
    SFDi(i,:) = cumsum(loads); % up-down method for SFD 
    % BMD = num. integral(SFD)
    BMDi(i,:) = cumtrapz(SFDi(i,:)); % integral of BMD from SFD
end

% this is the SFD and the positions of the train that gives the max SFD at that point
[SFD, max_loads] = max(abs(SFDi));   % SFD envelope
BMD = max(BMDi);        % BMD envelope

max(SFD)
max(BMD)
figure
plot(BMD)


% to use max_loads, plot it
% x-axis is position along the bridge, y-axis is
% the position of the rightmost axel relative to x=0
% you can click on the plot to get the point
% compare to the SFE and BME to get what the
% shear/moment is at that point
