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
        avg_error_ratios = zeros(num_trials, 1); % For storing error ratios
        iterations = zeros(num_trials, 1); % For storing the number of iterations
        max_quantities = zeros(num_trials, 1); % For storing max quantities
        backward_errors = zeros(num_trials, 1); % For storing backward errors
        
        for trial = 1:num_trials
            [U, Sigma, V, A] = generate_complex_Matrix(m, n, p);
            [Uini, H, it, Uout, singvals, Vout] = svd_og(A);
            [V1, V2] = first_call(H);
            
            % Store the number of iterations
            iterations(trial) = it;
            
            % Compute singular value errors
            sv_errors = (diag(Sigma) - diag(singvals)).^2;
            avg_sv_errors(trial) = mean(sv_errors);
            
            % Compute errors in left singular vectors
            U_error = abs(U(:, 1:n) - Uout).^2;
            avg_U_errors(trial) = mean(U_error(:));
            
            % Compute errors in right singular vectors
            V_error = abs(V - Vout).^2;
            avg_V_errors(trial) = mean(V_error(:));
            
            % Compute the error ratio
            avg_error_ratios(trial) = computeErrorRatio(V1, V2, H);
            
            % Compute max orthogonality error
            max_quantities(trial) = computeMaxOrthogonality(Uout, Vout);
            
            % Compute backward error
            backward_errors(trial) = computeBackwardError(A, Uout, Sigma, Vout);
        end
        
        % Compute average errors over all trials
        final_avg_sv_error = mean(avg_sv_errors);
        final_avg_U_error = mean(avg_U_errors);
        final_avg_V_error = mean(avg_V_errors);
        final_avg_error_ratio = mean(avg_error_ratios); % Average error ratio
        final_avg_iterations = mean(iterations); % Average number of iterations
        final_avg_max_quantity = mean(max_quantities); % Average max quantity
        final_avg_backward_error = mean(backward_errors); % Average backward error
        
        % Display the results
        fprintf('Results for power %d:\n', p);
        fprintf('Average squared error in singular values: %f\n', final_avg_sv_error);
        fprintf('Minimum squared error in singular values: %f\n', min(avg_sv_errors));
        fprintf('Maximum squared error in singular values: %f\n', max(avg_sv_errors));
        
        fprintf('Average squared error in left singular vectors: %f\n', final_avg_U_error);
        fprintf('Minimum squared error in left singular vectors: %f\n', min(avg_U_errors));
        fprintf('Maximum squared error in left singular vectors: %f\n', max(avg_U_errors));
        
        fprintf('Average squared error in right singular vectors: %f\n', final_avg_V_error);
        fprintf('Minimum squared error in right singular vectors: %f\n', min(avg_V_errors));
        fprintf('Maximum squared error in right singular vectors: %f\n', max(avg_V_errors));
        
        fprintf('Average berr: %e\n', final_avg_error_ratio);
        fprintf('Minimum berr: %e\n', min(avg_error_ratios));
        fprintf('Maximum berr: %e\n', max(avg_error_ratios));
        
        fprintf('Average max orthogonality error: %e\n', final_avg_max_quantity);
        fprintf('Minimum max orthogonality error: %e\n', min(max_quantities));
        fprintf('Maximum max orthogonality error: %e\n', max(max_quantities));
        
        fprintf('Average backward error: %e\n', final_avg_backward_error);
        fprintf('Minimum backward error: %e\n', min(backward_errors));
        fprintf('Maximum backward error: %e\n', max(backward_errors));
        
        fprintf('Average number of iterations: %f\n', final_avg_iterations);
        fprintf('Minimum number of iterations: %f\n', min(iterations));
        fprintf('Maximum number of iterations: %f\n', max(iterations));
        
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

function ratio = computeErrorRatio(V1, V2, A)
    % Compute the block matrix
    V = [V1, V2];
    
    % Compute the matrix product
    AV = V' * A * V;
    
    % Extract blocks from AV
    n1 = size(V1, 2);
    n2 = size(V2, 2);
    
    A1 = AV(1:n1, 1:n1);
    E = AV(1:n1, n1+1:n1+n2);
    E_star = AV(n1+1:n1+n2, 1:n1);
    A2 = AV(n1+1:n1+n2, n1+1:n1+n2);
    % Check if E has any nonzero elements
    
    % Compute the Frobenius norm of E
    normE = norm(E, 'fro');

    
    % Compute the Frobenius norm of A
    normA = norm(A, 'fro');

    
    % Compute the ratio
    ratio = normE / normA;

end

function backward_error = computeBackwardError(A, U, Sigma, V)
    % Compute the reconstructed matrix from U, Sigma, V
    A_reconstructed = U * Sigma * V';
    
    % Compute the Frobenius norm of the error matrix
    error_matrix = A - A_reconstructed;
    norm_error = norm(error_matrix, 'fro');
    
    % Compute the Frobenius norm of the original matrix A
    norm_A = norm(A, 'fro');
    
    % Compute the backward error
    backward_error = norm_error / norm_A;
end

function max_quantity = computeMaxOrthogonality(U, V)
    % Get the number of columns in U and V
    n = 10000;
    m = size(U, 2);
    
    % Compute the orthogonality errors
    U_error = U' * U - eye(m);
    V_error = V' * V - eye(m);
    
    % Compute the Frobenius norm of the errors
    norm_U_error = norm(U_error, 'fro');
    norm_V_error = norm(V_error, 'fro');
    
    % Compute the scaled errors
    scaled_U_error = norm_U_error / sqrt(n);
    scaled_V_error = norm_V_error / sqrt(n);
    
    % Compute the maximum of the scaled errors
    max_quantity = max(scaled_U_error, scaled_V_error);
    
    % Print the result
end

