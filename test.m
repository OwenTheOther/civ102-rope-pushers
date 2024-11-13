clear; close all;

%% 0. Initialize Parameters
L = 1260;               % Length of bridge
n = L;               % Discretize into 1 mm seg.
P = 400;                % Total weight of train [N]
x = linspace(0, L, n); % x-axis
A_pos = 40;            % position of first support from left edge
B_pos = A_pos + 1200;  % position of second support from left edge


%% 1. SFD, BMD under train loading
x_train = [52 228 392 568 732 908]; % Train Load Locations
P_train = [1 1 1 1 1 1] * P/6;

n_train = 900;                 % num of train locations
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
    loads([A_pos B_pos]) = [rxn_A rxn_B];

    %loadi stuff
    loadsi(i,:) = loads;
    % SFD = num. integral(w)
    SFDi(i,:) = cumsum(loads);
    % BMD = num. integral(SFD)
    BMDi(i,:) = cumtrapz(SFDi(i,:));

end

SFD = max(abs(SFDi));   % SFD envelope
BMD = max(BMDi);        % BMD envelope

plot(BMD)