function [U, A_hat, eigenvalues] = RxDAWN( ...
    X, stimulusOnsets, Fs, P300Length, Cv, alpha)
%RXDAWN Compute regularized xDAWN spatial filters.
%
%   [U, A_HAT, EIGENVALUES] = RXDAWN(X, STIMULUSONSETS, FS, ...
%       P300LENGTH, CV, ALPHA) estimates spatial filters that enhance
%   target P300 responses while penalizing stimulus-related VEP activity.
%
%   Inputs
%   ------
%   X
%       EEG data matrix of size Nt-by-Ns, where Nt is the number of
%       temporal samples and Ns is the number of EEG channels.
%
%       Each channel should be zero-mean before calling this function.
%
%   stimulusOnsets
%       Vector containing the sample indices of target-stimulus onsets.
%
%   Fs
%       Sampling frequency in Hz.
%
%   P300Length
%       Length of the modeled P300 response in seconds, typically
%       between 0.6 and 1.0 seconds.
%
%   Cv
%       VEP covariance or scatter matrix of size Ns-by-Ns.
%
%   alpha
%       Regularization weight in the interval [0, 1].
%
%       alpha = 0 produces the conventional xDAWN formulation.
%       Increasing alpha places greater emphasis on suppressing VEPs.
%
%   Outputs
%   -------
%   U
%       Spatial-filter matrix of size Ns-by-Ns. Columns are ordered by
%       decreasing generalized eigenvalue. The first columns represent
%       the spatial filters with the highest RxDAWN objective values.
%
%   A_hat
%       Estimated target ERP response of size Ne-by-Ns, where Ne is the
%       number of samples in the modeled P300 response.
%
%   eigenvalues
%       Generalized eigenvalues sorted in descending order.
%
%   Reference
%   ---------
%   A. Mobaien, R. Boostani, and S. Sanei,
%   "Improving the performance of P300-based BCIs by mitigating the
%   effects of stimuli-related evoked potentials through regularized
%   spatial filtering," Journal of Neural Engineering, vol. 21, no. 1,
%   016023, 2024.
%
%   Copyright (c) 2026 Ali Mobaien
%   Released under the MIT License.

    %% Validate inputs

    if nargin ~= 6
        error('RxDAWN requires exactly six input arguments.');
    end

    if ~isnumeric(X) || ~ismatrix(X) || isempty(X)
        error('X must be a nonempty numeric matrix.');
    end

    if any(~isfinite(X(:)))
        error('X must contain only finite values.');
    end

    [Nt, Ns] = size(X);

    if ~isnumeric(stimulusOnsets) || ~isvector(stimulusOnsets) || ...
            isempty(stimulusOnsets)
        error('stimulusOnsets must be a nonempty numeric vector.');
    end

    stimulusOnsets = stimulusOnsets(:);

    if any(~isfinite(stimulusOnsets)) || ...
            any(stimulusOnsets ~= round(stimulusOnsets))
        error('stimulusOnsets must contain integer sample indices.');
    end

    if any(stimulusOnsets < 1) || any(stimulusOnsets > Nt)
        error('All stimulus onsets must be between 1 and size(X,1).');
    end

    if ~isscalar(Fs) || ~isfinite(Fs) || Fs <= 0
        error('Fs must be a positive finite scalar.');
    end

    if ~isscalar(P300Length) || ~isfinite(P300Length) || ...
            P300Length <= 0
        error('P300Length must be a positive finite scalar.');
    end

    if ~isnumeric(Cv) || ~isequal(size(Cv), [Ns, Ns])
        error('Cv must be an Ns-by-Ns numeric matrix.');
    end

    if any(~isfinite(Cv(:)))
        error('Cv must contain only finite values.');
    end

    if ~isscalar(alpha) || ~isfinite(alpha) || ...
            alpha < 0 || alpha > 1
        error('alpha must be a scalar in the interval [0, 1].');
    end

    %% Construct the target-stimulus Toeplitz matrix

    Ne = round(P300Length * Fs);

    if Ne < 1
        error('P300Length and Fs must produce at least one ERP sample.');
    end

    firstColumn = zeros(Nt, 1);
    firstColumn(stimulusOnsets) = 1;

    firstRow = zeros(1, Ne);
    firstRow(1) = firstColumn(1);

    D = toeplitz(firstColumn, firstRow);

    %% Estimate the target ERP response

    DtD = D' * D;
    DtX = D' * X;

    % Least-squares estimate:
    % A_hat = (D' * D)^(-1) * D' * X
    %
    % Use a linear-system solution rather than an explicit inverse.
    if rcond(DtD) > eps
        A_hat = DtD \ DtX;
    else
        warning('RxDAWN:IllConditionedDesignMatrix', ...
            ['D''*D is singular or poorly conditioned. ', ...
             'Using the pseudoinverse to estimate A_hat.']);

        A_hat = pinv(D) * X;
    end

    %% Construct the matrices in the RxDAWN objective

    % Psi = A_hat' * D' * D * A_hat
    Psi = A_hat' * DtD * A_hat;

    Cx = X' * X;

    % Remove small numerical asymmetries.
    Psi = (Psi + Psi') / 2;
    Cx  = (Cx  + Cx')  / 2;
    Cv  = (Cv  + Cv')  / 2;

    %% Normalize Cx and Cv

    maxCxDiagonal = max(abs(diag(Cx)));
    maxCvDiagonal = max(abs(diag(Cv)));

    if maxCxDiagonal <= eps
        error('Cx cannot be normalized because its diagonal is zero.');
    end

    if maxCvDiagonal <= eps
        error('Cv cannot be normalized because its diagonal is zero.');
    end

    Cx = Cx / maxCxDiagonal;
    Cv = Cv / maxCvDiagonal;

    %% Solve the generalized eigenvalue problem

    B = (1 - alpha) * Cx + alpha * Cv;
    B = (B + B') / 2;

    % A*u = lambda*B*u
    [U, eigenvalueMatrix] = eig(Psi, B);

    eigenvalues = real(diag(eigenvalueMatrix));

    [eigenvalues, sortIndex] = sort(eigenvalues, 'descend');
    U = real(U(:, sortIndex));

end
