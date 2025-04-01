function [F, p, resid, df_model, df_error] = F_test_full_vs_red(y, X, Xred, df)
%function [F, p, resid, df_model, df_error] = F_test_full_vs_red(y, X, Xred, px, pxred)
% :Usage:
% ::
%
%     [F, p, resid, df_model, df_error] = F_test_full_vs_red(y, X, Xred, px, pxred)
%
% :Examples:
% ::
%
%    X = randn(100, 3); Xred = X(:,1); y = X(:,2) + randn(100, 1);
%    px = pinv(X); pxred = pinv(Xred);
%    [F, p, resid] = F_test_full_vs_red(y, X, Xred, px, pxred);
%
%
%    % Test full-model F-value against regress.m
%    Xred = X(:,end); % intercept only
%    px = pinv(X);
%    pxred = pinv(Xred);
%    [F, p, resid, dfm, dfe] = F_test_full_vs_red(y, X, Xred, px, pxred); % full model F-test
%    [b, bint, r, rint, stats] = regress(y, X);
%
% ..
%    Tested OK on 11/27/07, by tor
%   From Tor Wager
% .. Taken from CanlabCore on 2022/01/26

    % Precompute these since they are used regularly
    px = pinv(X);
    pxred = pinv(Xred);
    
    T = length(y);      % Length of time course

    k = size(px, 1);        % predictors: full model
    kred = size(pxred, 1);  % predictors: reduced model

    % Degrees of freedom: model: Full - reduced
    df_model = k - kred;                        % degrees of freedom for model (param - 1)
    df_error = T - k;                           % error DF, full model
    
    df_model = df; % No way to programmatically determine this, so passing as argument
    df_error = T-k; % error DF, full model


    % preallocate space
    [F, p] = deal(NaN(size(y, 1), size(y, 2)));
    h = waitbar(0,'Please wait ... computing stat mask');
    for row = 1:size(y, 1)
        waitbar(row / size(y, 1))
        for col = 1:size(y, 2)
            % Step 1: Find the OLS solution for the FULL model
            % ---------------------------------------------------
            beta = px * squeeze(y(row, col, :));                             % Beta values
            resid = squeeze(y(row, col, :)) - X * beta;                       % Residuals

            % Sums of squares
            %SST = y' * y;
            SSE = resid' * resid;
            %SSfull = SST - SSE;
            var_est = SSE / df_error;                   % Estimate of Sigma^2

            % Step 2: Find the OLS solution for the REDUCED model
            % F-test for full vs. reduced
            % ---------------------------------------------------
            betared = pxred * squeeze(y(row, col, :));                          % Beta values
            residred = squeeze(y(row, col, :)) - Xred * betared;
            SSEred = residred' * residred;
            %SSred = betared' * Xred' * y;                      % Full model sum of squares

            % F stat
            % (SSred - SSfull) ./ (var * df_model)
            F(row, col) = (SSEred - SSE) / (var_est * df_model);         % F-statistic - compare with F-distribution with (param-1, df) degrees of freedom

            p(row, col) = 1 - fcdf(F(row, col), df_model, df_error);
        end
    end
     close(h)
end
