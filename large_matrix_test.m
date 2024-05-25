function large_matrix_test(m, n)
    % Generate matrix A for different powers p
    powers = [1, 2, 3]; % Corresponds to k, k^2, k^3
    num_trials = 10; % Number of trials for averaging
    % Start the timer
    total_timer = tic;
    
    for p = powers
        fprintf('Considering power: %d:\n', p);
        
        avg_sv_errors = zeros(num_trials, 1);
        avg_U_errors = zeros(num_trials, 1);
        avg_V_errors = zeros(num_trials, 1);
        
        for trial = 1:num_trials
            [U, Sigma, V, A] = generate_complex_Matrix(m, n, p);
            [Uout, singvals, Vout] = qdwhsvd(A);
            
            % Compute singular value errors
            sv_errors = (diag(Sigma) - diag(singvals)).^2;
            avg_sv_errors(trial) = mean(sv_errors);
            
            % Compute errors in left singular vectors
            U_error = abs(U(:, 1:n) - Uout).^2;
            avg_U_errors(trial) = mean(U_error(:));
            
            % Compute errors in right singular vectors
            V_error = abs(V - Vout).^2;
            avg_V_errors(trial) = mean(V_error(:));
        end
        
        % Compute average errors over all trials
        final_avg_sv_error = mean(avg_sv_errors);
        final_avg_U_error = mean(avg_U_errors);
        final_avg_V_error = mean(avg_V_errors);
        
        % Display the results
        fprintf('Results for power %d:\n', p);
        fprintf('Average squared error in singular values: %f\n', final_avg_sv_error);
        fprintf('Average squared error in left singular vectors: %f\n', final_avg_U_error);
        fprintf('Average squared error in right singular vectors: %f\n', final_avg_V_error);
        fprintf('Elapsed time: %f seconds.\n', toc(total_timer));
    end
end

function [U, Sigma, V, A] = generateMatrix(m, n, p)
    % Generate random orthogonal matrices U (mxm) and V (nxn)
    [U, ~] = qr(randn(m));
    [V, ~] = qr(randn(n));
    
    % Create the singular values decreasing according to 1/k^p
    k = 1:n; % Smaller dimension for the diagonal matrix
    sigma = 1 ./ (k.^p);
    Sigma = diag(sigma);
    
    % Construct the matrix A with dimensions mxn
    A = U(:, 1:n) * Sigma * V'; % U reduced to mxn by taking only the first n columns
end

function [U, Sigma, V, A] = generate_complex_Matrix(m, n, p)
    % Generate random complex orthogonal matrices U (mxm) and V (nxn)
    [U, ~] = qr(randn(m) + 1i*randn(m));
    [V, ~] = qr(randn(n) + 1i*randn(n));
    
    % Create the singular values decreasing according to 1/k^p
    k = 1:n; % Smaller dimension for the diagonal matrix
    sigma = 1 ./ (k.^p);
    Sigma = diag(sigma);
    
    % Construct the matrix A with dimensions mxn
    A = U(:, 1:n) * Sigma * V'; % U reduced to mxn by taking only the first n columns
end

function dims(U)
    [m, n] = size(U);
    
    % Print the dimensions
    fprintf('Matrix dimensions: %d x %d\n', m, n);
end
