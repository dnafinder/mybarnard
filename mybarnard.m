function STATS = mybarnard(x, varargin)
%MYBARNARD  Barnard's Exact Test for 2×2 Contingency Tables (vectorized).
%
%   STATS = MYBARNARD(X)
%   STATS = MYBARNARD(X, PLTS)
%   STATS = MYBARNARD(X, PLTS, Tbx)
%
%   This function performs Barnard’s unconditional exact test for a 2×2
%   table using a fully vectorized implementation (no loops in the core
%   computational engine). This produces extremely fast evaluation even for
%   fine grids of nuisance parameters.
%
%   Barnard's test is strictly more powerful than Fisher’s exact test for
%   2×2 tables, but historically less used due to computational burden.
%   This implementation removes that limitation.
%
%   ----------------------------------------------------------------------
%   INPUTS
%   ----------------------------------------------------------------------
%   X      : 2×2 matrix of nonnegative integers.
%
%   PLTS   : (optional) plotting flag  
%            0 = no plots  (default)  
%            1 = show diagnostic plots
%
%   Tbx    : (optional) number of grid points for nuisance parameter np  
%            Default = 100  
%            Increasing Tbx increases precision at the cost of time.
%
%   ----------------------------------------------------------------------
%   OUTPUT
%   ----------------------------------------------------------------------
%   If no output is requested, prints a table with:
%       - Tables (Tbx)
%       - Size (column totals + 1)
%       - Wald statistic of observed table
%       - Nuisance parameter where P(np) is maximized
%       - One-tailed Barnard p-value = max(P(np))
%       - Two-tailed p-value = min(2*p, 1)
%
%   If STATS is requested:
%       STATS.TX0      → Wald statistic of observed table
%       STATS.p_value  → maximized p-value (one-tailed)
%       STATS.nuisance → nuisance parameter np*
%
%   ----------------------------------------------------------------------
%   EXAMPLE
%   ----------------------------------------------------------------------
%   X = [7 12; 8 3];
%   mybarnard(X)
%
%   ----------------------------------------------------------------------
%   REFERENCES
%   ----------------------------------------------------------------------
%   Barnard, G.A. (1945). "A New Test for 2×2 Tables". Nature.
%   Fisher, R.A. (1925). Statistical Methods for Research Workers.
%
%   ----------------------------------------------------------------------
%   CITATION
%   ----------------------------------------------------------------------
%   If you use this function in academic or technical work, please cite:
%
%     Cardillo G. (2009)
%     "MyBarnard: a very compact routine for Barnard's exact test on 2×2 matrix".
%     Available from GitHub:
%     https://github.com/dnafinder/mybarnard
%
%   ----------------------------------------------------------------------
%   METADATA
%   ----------------------------------------------------------------------
%   Author:  Giuseppe Cardillo
%   Email:   giuseppe.cardillo.75@gmail.com
%   GitHub:  https://github.com/dnafinder
%   Created: 2009-01-01
%   Updated: 2025-11-21
%   Version: 2.0.0
%
%   ----------------------------------------------------------------------

% -------------------------------------------------------------------------
% Input parsing and validation
% -------------------------------------------------------------------------
parser = inputParser;
parser.FunctionName = mfilename;

addRequired(parser, 'x', @(v) validateattributes(v, {'numeric'}, ...
    {'real','finite','integer','nonnegative','nonnan','size',[2 2]}, ...
    mfilename, 'x', 1));

addOptional(parser, 'plts', 0, @(v) isnumeric(v) && isscalar(v) && ...
    isfinite(v) && isreal(v) && any(v == [0 1]));

addOptional(parser, 'Tbx', 100, @(v) validateattributes(v, {'numeric'}, ...
    {'real','finite','integer','positive','nonnan','scalar'}, ...
    mfilename, 'Tbx', 3));

parse(parser, x, varargin{:});
plotFlag = parser.Results.plts;
Tbx      = parser.Results.Tbx;

% -------------------------------------------------------------------------
% Basic quantities and possible column reordering
% -------------------------------------------------------------------------
colSums = sum(x, 1);          % column totals
N       = sum(colSums);       % total sample size

% Reorder columns if necessary (ascending column totals) for numerical
% convenience; test is symmetric with respect to column swapping.
if ~issorted(colSums)
    x       = fliplr(x);
    colSums = sort(colSums);
end

% -------------------------------------------------------------------------
% Wald statistic on the full grid and observed table
% -------------------------------------------------------------------------
% Define all possible (I,J) pairs:
%   I = 0:colSums(1), J = 0:colSums(2)
[I, J] = ndgrid(0:colSums(1), 0:colSums(2));

% Wald statistic:
%   T(X) = (p(a) - p(b)) / sqrt(p*(1-p)*(1/C1 + 1/C2))
% where p(a) = I/C1, p(b) = J/C2, p = (I+J)/N.
pHat    = (I + J) ./ N;
denom   = sqrt( pHat .* (1 - pHat) .* sum(1 ./ colSums) );
TX      = (I./colSums(1) - J./colSums(2)) ./ denom;

% Fix singularities: I=0,J=0 and I=Cs1,J=Cs2 (0/0 => set to 0)
TX(1)   = 0;
TX(end) = 0;

% Wald statistic for the observed table.
% NOTE: follow the original implementation:
%   I0 = x(1,1)
%   J0 = x(1,2) == x(3) in linear indexing of the 2x2 table
TX0 = abs(TX(x(1) + 1, x(3) + 1));

% Indicator for all tables with T(X) >= T(X0)
idx = (TX >= TX0);

% -------------------------------------------------------------------------
% Log-probabilities under the nuisance parameter np
% -------------------------------------------------------------------------
% Under H0, with nuisance np:
%
%   P(X | np) = C(Cs1, I) * C(Cs2, J) * np^(I+J) * (1 - np)^(N - (I+J))
%
% We work in log-space and then exponentiate to avoid overflow.
%
% log P(X | np) = CF(I,J) + (I+J)*log(np) + [N - (I+J)]*log(1-np)
%
% where CF(I,J) collects all terms that do not depend on np.

A    = [1 1 Tbx];                 % replication pattern for 3D boxes
Bdim = colSums + 1;               % dimensionality of the I,J grid
npa  = linspace(0.0001, 0.9999, Tbx);

logNp    = log(npa);
log1mNp  = log(1 - npa);

E = repmat(I + J, A);             % (Cs1+1) x (Cs2+1) x Tbx
F = N - E;                        % same size as E

% Constant factorial/log-gamma term (independent of np)
CF2D = sum(gammaln(colSums + 1)) ...
       - (gammaln(I + 1) + gammaln(J + 1) ...
       + gammaln(colSums(1) - I + 1) + gammaln(colSums(2) - J + 1));
CF   = repmat(CF2D, A);           % 3D "box" replicated Tbx times

% Build the terms involving np as 3D arrays:
logNpBox   = repmat(reshape(logNp, A), Bdim);
log1mNpBox = repmat(reshape(log1mNp, A), Bdim);

% Final probabilities for each (I,J,np): S is 3D
S = exp(CF + E .* logNpBox + F .* log1mNpBox);

% -------------------------------------------------------------------------
% Compute P(np) = sum_{T(X) >= T(X0)} P(X | np)
% -------------------------------------------------------------------------
% We need, for each np, the sum of P(X|np) over all tables with TX>=TX0.
L = nnz(idx);  % number of tables where TX>=TX0

% S(repmat(idx,A)) is a vector of length L*Tbx; reshape to L-by-Tbx and sum
P = sum(reshape(S(repmat(idx, A)), L, Tbx), 1);

% Maximized p-value (Barnard's p) and corresponding nuisance parameter
[PV, maxIdx] = max(P);
npStar       = npa(maxIdx);

% -------------------------------------------------------------------------
% Display results
% -------------------------------------------------------------------------
resultsTable = table( ...
    Tbx, ...
    colSums + 1, ...
    TX0, ...
    npStar, ...
    PV, ...
    min(2*PV, 1), ...
    'VariableNames', { ...
    'Tables', ...
    'Size', ...
    'Wald_stat', ...
    'Nuisance', ...
    'one_tailed_p_value', ...
    'two_tailed_p_value'});

disp(resultsTable);

% Optional structured output
if nargout > 0
    STATS.TX0      = TX0;
    STATS.p_value  = PV;
    STATS.nuisance = npStar;
end

% -------------------------------------------------------------------------
% Optional plots
% -------------------------------------------------------------------------
if plotFlag
    scrsz = get(groot,'ScreenSize');
    figure('Color',[1 1 1],'OuterPosition',scrsz);

    % --- Plot 1: p-value as a function of nuisance parameter np ---------
    subplot(1,2,1);
    hold on;
    patch([0 npa 1], [0 P 0], 'b');
    plot([0 npStar], [PV PV], 'k--');
    plot([npStar npStar], [0 PV], 'w--');
    plot(npStar, PV, 'ro', 'MarkerFaceColor', 'r');
    hold off;
    title('Barnard''s exact p-value vs nuisance parameter np', ...
        'FontName','Arial', 'FontSize',14, 'FontWeight','Bold');
    xlabel('Nuisance parameter (np)', ...
        'FontName','Arial', 'FontSize',14, 'FontWeight','Bold');
    % <-- QUI il fix: solo testo semplice, niente \ge o \mid
    ylabel('p-value  [P(T(X) >= T(X_0) | np)]', ...
        'FontName','Arial', 'FontSize',14, 'FontWeight','Bold');
    axis square;

    % --- Plot 2: distribution of T(X) under np* --------------------------
    txtLabel = sprintf('P[T(X) | np = %.5g]', npStar);

    % Flatten TX and the slice of S corresponding to np*
    Aflat = TX(:);
    Sflat = reshape(S(:,:,maxIdx), numel(TX), 1);

    % Group by unique T(X) values and sum probabilities
    [Asorted, sortIdx] = sort(Aflat);
    Bsorted            = Sflat(sortIdx);

    uniqueTX = unique(Asorted);
    Yplot    = zeros(size(uniqueTX));

    for k = 1:numel(uniqueTX)
        Yplot(k) = sum(Bsorted(Asorted == uniqueTX(k)));
    end

    subplot(1,2,2);
    bar(Yplot);  % T(X) values are encoded on the x-axis indices
    axis square tight;

    xt = get(gca, 'XTick');
    xt = xt(xt >= 1 & xt <= numel(uniqueTX));  % versione corretta con & vettoriale
    xtl = arrayfun(@(v) sprintf('%.1f', uniqueTX(v)), xt, 'UniformOutput', false);
    set(gca, 'XTick', xt, 'XTickLabel', xtl);

    title('Distribution of Wald Statistic for Barnard''s exact test', ...
        'FontName','Arial', 'FontSize',12, 'FontWeight','Bold');
    xlabel('T(X)', ...
        'FontName','Arial', 'FontSize',12, 'FontWeight','Bold');
    ylabel(txtLabel, ...
        'FontSize',12, 'FontWeight','Bold');
end

end
