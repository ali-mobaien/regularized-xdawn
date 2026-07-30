%% Example: Regularized xDAWN for P300-based BCI data
%
% Expected data dimensions:
%
% Signal       : [Nt × Ns × N_char]
% StimulusType : [N_char × Nt]
% Cv           : [Ns × Ns]
%
% Nt      = number of temporal samples per character
% Ns      = number of EEG channels
% N_char  = number of characters or recording sequences
%
% The example assumes that the training data and Cv have already
% been loaded into the MATLAB workspace.

clc;
clear;
close all;

%% Parameters

Fs = 240;               % Sampling frequency in Hz
P300Length = 0.8;       % Analyzed ERP duration in seconds
alpha = 0.05;           % RxDAWN regularization parameter
numberOfFilters = 4;    % Number of spatial filters to retain

%% Load data
%
% Replace this section with the location of your own data.
%
% The loaded file should contain:
%   Signal       : EEG data [Nt × Ns × N_char]
%   StimulusType : target/non-target labels [N_char × Nt]
%
% load('Subject_B_Train.mat');
%
% If the original Signal dimensions are:
%   [N_char × Nt × Ns]
%
% convert them using:
%
% Signal = permute(Signal, [2 3 1]);

Signal = double(Signal);

[Nt, Ns, N_char] = size(Signal);

%% Optional band-pass filtering

applyBandpassFilter = true;
filterOrder = 4;
frequencyBand = [0.1 15];

if applyBandpassFilter

    [b, a] = butter( ...
        filterOrder, frequencyBand / (Fs / 2), 'bandpass');

    % filtfilt operates along the first dimension, which corresponds
    % to time in Signal [Nt × Ns × N_char].
    Signal = filtfilt(b, a, Signal);

end

%% Zero-mean normalization
%
% Remove the temporal mean independently from every channel and
% every character recording.

channelMean = mean(Signal, 1);

Signal = Signal - repmat(channelMean, [Nt 1 1]);

%% Unit-variance normalization
%
% Normalize each EEG channel independently within each character
% recording.

for characterIndex = 1:N_char

    channelStd = std(Signal(:, :, characterIndex), 0, 1);

    % Prevent division by zero for constant or inactive channels.
    channelStd(channelStd < eps) = 1;

    Signal(:, :, characterIndex) = ...
        Signal(:, :, characterIndex) ./ channelStd;

end

%% Define the VEP penalty matrix
%
% Cv must be an Ns-by-Ns symmetric matrix representing the spatial
% contribution of stimulus-related VEP activity.
%
% Replace this example with the Cv construction used in your study.
%
% Example:
%
% load('VEP_energy.mat', 'energy_vep');
% Cv = diag(energy_vep);

if ~exist('Cv', 'var')
    error(['Cv is not available. Load or calculate an ', ...
           'Ns-by-Ns VEP penalty matrix before running this example.']);
end

if ~isequal(size(Cv), [Ns Ns])
    error('Cv must have dimensions Ns-by-Ns.');
end

%% Concatenate all character recordings
%
% Original:
%   Signal [Nt × Ns × N_char]
%
% After permutation:
%   Signal [Nt × N_char × Ns]
%
% After reshaping:
%   X [Nt*N_char × Ns]

X = permute(Signal, [1 3 2]);
X = reshape(X, Nt * N_char, Ns);

%% Construct target-stimulus onset vector
%
% StimulusType is expected to have dimensions:
%   [N_char × Nt]
%
% The following operation identifies transitions corresponding to
% target-stimulus onsets.

targetOnsetMask = ...
    conv2(1, [0 1 -1], 2 * (StimulusType == 1), 'same') > 1;

% Arrange the target mask in the same temporal order used when the
% EEG recordings were concatenated.
targetOnsetMask = targetOnsetMask';
targetOnsetMask = targetOnsetMask(:);

% RxDAWN expects sample indices rather than a logical indicator vector.
targetStimulusOnsets = find(targetOnsetMask);

%% Calculate RxDAWN spatial filters

[U, A_hat, eigenvalues] = RxDAWN( ...
    X, ...
    targetStimulusOnsets, ...
    Fs, ...
    P300Length, ...
    Cv, ...
    alpha);

%% Select the leading spatial filters

numberOfFilters = min(numberOfFilters, size(U, 2));

selectedFilters = U(:, 1:numberOfFilters);

%% Project the concatenated EEG data

filteredX = X * selectedFilters;

%% Restore the original character-based organization
%
% Output dimensions:
%   filteredSignal [Nt × numberOfFilters × N_char]

filteredSignal = reshape( ...
    filteredX, Nt, N_char, numberOfFilters);

filteredSignal = permute(filteredSignal, [1 3 2]);

%% Display results

fprintf('RxDAWN completed successfully.\n');
fprintf('Input channels: %d\n', Ns);
fprintf('Selected spatial filters: %d\n', numberOfFilters);
fprintf('Filtered signal dimensions: [%d × %d × %d]\n', ...
    size(filteredSignal, 1), ...
    size(filteredSignal, 2), ...
    size(filteredSignal, 3));

disp('Leading generalized eigenvalues:');
disp(eigenvalues(1:numberOfFilters));

%% Optional visualization of estimated ERP components

timeVector = (0:size(A_hat, 1)-1) / Fs;

estimatedERPComponents = A_hat * selectedFilters;

figure;

plot(timeVector, estimatedERPComponents, 'LineWidth', 1.2);

xlabel('Time (s)');
ylabel('Amplitude (normalized units)');
title('Estimated ERP Components after RxDAWN');
grid on;
