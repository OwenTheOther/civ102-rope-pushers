clear; close all;

%% 0. Initialize Parameters
L = 1200;               % Length of bridge
n = L;                  % Discretize into 1 mm seg.
P = 400;                % Total weight of train [N]
x = 1:1:L;              % x-axis

% Note: the position is 1mm greater then the physical location from the left
% edge of the bridge; A_pos = 1 means the location 0mm from the left of the
% bridge.
A_pos = 1;            % position of first support from left edge
B_pos = 1200; % position of second support from left edge


%% 1. SFD, BMD under train loading
x_train = [52 228 392 568 732 908]; % Train Load Locations
P_train = [1 1 1 1 1 1] * P/6;

n_train = 349;                 % num of train locations
train_space = floor(linspace(1, L + max(x_train), n_train));
loadsi = zeros(n_train, n);    % 1 load plot for each train loc.
SFDi = zeros(n_train, n);     % 1 SFD for each train loc.
BMDi = zeros(n_train, n);     % 1 BMD for each train loc.

% Solve for SFD and BMD with the train at different locations
for J = 1:n_train
    i = train_space(J);
    % start location of train
    ttrain_pos_i = x_train - max(x_train) + i;
    ttrain_pos_i(ttrain_pos_i <= 0) = NaN;
    ttrain_pos_i(ttrain_pos_i >= L) = NaN;

    p_ts = P_train(~isnan(ttrain_pos_i));
    pos_ti = ttrain_pos_i(~isnan(ttrain_pos_i));

    % sum of moments at A eqn
    rxn_B = sum((A_pos - pos_ti) .* p_ts) / -1200;

    % sum of Fy eqn
    rxn_A = sum(p_ts) - rxn_B;
    % construct applied loads
    % w(x)
    loads = zeros(1, L);
    loads(pos_ti) = -1 .* p_ts;
    loads([A_pos+1 B_pos]) = [rxn_A rxn_B];

    %loadi stuff
    loadsi(i,:) = loads;
    % SFD = num. integral(w)
    SFDi(i,:) = cumsum(loads);
    % BMD = num. integral(SFD)
    BMDi(i,:) = cumtrapz(SFDi(i,:));

end

[SFD, max_loadsSFD] = max(abs(SFDi));   % SFD envelope
[BMD, max_loadsBMD] = max(BMDi);        % BMD envelope

% to use max_loads, plot it
% x-axis is position along the bridge, y-axis is
% the position of the rightmost axel relative to x=0
% you can click on the plot to get the point
% compare to the SFE and BME to get what the
% shear/moment is at that point
figure(1)
plot(SFD)

figure(2)
plot(BMD)

figure(3)
plot(max_loadsSFD)

figure(4)
plot(max_loadsBMD)

%max shear moment (using 600 because bridge is symmetrical and the code gets weird towards right end)
max_shear = max(SFD(1:600))
location_on_bridge = find(round(SFD) == round(max_shear))
train_rightmost_axel = max_loadsSFD(location_on_bridge)



%max bending moment
max_moment = max(BMD(1:600))
location_on_bridge = find(round(BMD) == round(max_moment))
train_rightmost_axel = max_loadsBMD(location_on_bridge)
