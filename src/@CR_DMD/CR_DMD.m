classdef CR_DMD < handle
    % CR_DMD : Class for Comprehensive Robust Dynamic Mode Decomposition

    properties
        params                % Struct : epsilon, eta, w, mu, max_iter and criteria for pre, me, and dr, respective algorithms
        use_gpu = false       % Boolean : Flag to indicate whether to use GPU acceleration
        verbose = true        % Boolean : Flag to indicate whether to print verbose output
        plot_progress = false % Boolean : Flag to indicate whether to plot progress during iterations
        display_freq = 10     % Integer : Frequency of displaying progress
        prev_msg_len = 0;     % Integer : Previous message length for clearing
    end

    methods
        function obj = CR_DMD(params)
            % Check if GPU is available
            if canUseGPU()
                obj.use_gpu = true;
                if isfield(params, 'verbose'), obj.verbose = params.verbose; end
                if obj.verbose, fprintf('GPU detected and will be used for computations.\n'); end
            else
                obj.use_gpu = false;
                if isfield(params, 'verbose'), obj.verbose = params.verbose; end
                if obj.verbose, fprintf('GPU not detected. Using CPU.\n'); end
            end

            % Initialize parameters
            obj.params = obj.init_params(params);
        end

        function p = init_params(obj, cp)
            % --- Preprocessin Default Parameters ---
            p.pre.w = 0.5;
            p.pre.epsilon = 0.1; 
            p.pre.eta = 0.1;
            p.pre.max_iter = 5000;
            p.pre.criteria = 1e-4;
            p.pre.auto_step = true;

            % --- Mode Extraction Defaults ---
            p.me.r = 21; % Target rank
            
            % --- Dimensional Reduction Defaults ---
            p.dr.w = 0.5;
            p.dr.epsilon = 0.1;
            p.dr.eta = 0.1;
            p.dr.mu = 1;
            p.dr.max_iter = 5000;
            p.dr.criteria = 1e-4;
            p.dr.auto_step = true;

            % Override defaults with user-provided parameters
            if nargin > 1 && isstruct(cp)
                if isfield(cp, 'verbose'), obj.verbose = cp.verbose; end
                if isfield(cp, 'pre')
                    p.pre = obj.merge_structs(p.pre, cp.pre);
                    if isfield(cp.pre, 'gamma_x') || isfield(cp.pre, 'gamma_s') || isfield(cp.pre, 'gamma_z1') || isfield(cp.pre, 'gamma_z2')
                        warning('CR_DMD:init_params', 'Setting gamma parameters directly is deprecated. COmform if your step sizes meet convergence criteria.');
                        p.pre.auto_step = false;
                    end
                end
                if isfield(cp, 'me'), p.me = obj.merge_structs(p.me, cp.me); end
                if isfield(cp, 'dr')
                    p.dr = obj.merge_structs(p.dr, cp.dr); 
                    if isfield(cp.dr, 'gamma_x') || isfield(cp.dr, 'gamma_s') || isfield(cp.dr, 'gamma_z1') || isfield(cp.dr, 'gamma_z2')
                        warning('CR_DMD:init_params', 'Setting gamma parameters directly is deprecated. COmform if your step sizes meet convergence criteria.');
                        p.dr.auto_step = false;
                    end
                end
            end
        end

        function [x_clean, s_outliers, Phi, Lambda, b] = preprocess(obj, x_noisy, h, v)
            % Algorithm 1 : Preprocessing
            % Inputs: x_noisy - Noisy observations (m x n)
            %         h - Number of rows in spatial dimension
            %         v - Number of columns in spatial dimension
            % Outputs: x_clean - Cleaned data (m x n)
            %          s_outliers - Detected outliers (m x n)
            %          Phi - DMD modes
            %          Lambda - DMD eigenvalues
            if obj.verbose, fprintf('Starting Algorithm 1: Preprocessing...\n'); end

            % Transport the input data if necessary
            if obj.use_gpu, x_noisy = gpuArray(x_noisy); end
            [N, M] = size(x_noisy);

            % Initialize variables
            x = x_noisy;
            s = zeros(N, M, 'like', x_noisy);
            z_1 = zeros(3*N, M, 'like', x_noisy);
            z_2 = x + s;

            % Calculate step sizes
            if obj.params.pre.auto_step
                steps = obj.calculate_auto_steps('pre');
                gamma_x = steps.gamma_x;
                gamma_s = steps.gamma_s;
                gamma_z1 = steps.gamma_z1;
                gamma_z2 = steps.gamma_z2;
            else
                gamma_x = obj.params.pre.gamma_x;
                gamma_s = obj.params.pre.gamma_s;
                gamma_z1 = obj.params.pre.gamma_z1;
                gamma_z2 = obj.params.pre.gamma_z2;
            end

            % Hyperparameters
            w = obj.params.pre.w;
            epsilon = obj.params.pre.epsilon;
            eta = obj.params.pre.eta;

            % Main iteration loop
            if obj.verbose, fprintf('\n'); end
            obj.prev_msg_len = 0; % Reset message length
            max_res = 0;
            for t = 1:obj.params.pre.max_iter
                x_old = x;
                s_old = s;

                % Update x
                grad_x = adj_diff_st(z_1, w, h, v, M) + z_2;
                x = x - gamma_x * grad_x;

                % Update s
                s = proj_l1ball(s - gamma_s * z_2, eta);

                % Extrapolation
                x_hat = 2 * x - x_old;
                s_hat = 2 * s - s_old;

                % Update z_1
                z1_tmp = z_1 + gamma_z1 * diff_st(x_hat, w, h, v);
                z_1 = z1_tmp - gamma_z1 * prox_l12(z1_tmp / gamma_z1, 1 / gamma_z1);

                % Update z_2
                z2_tmp = z_2 + gamma_z2 * (x_hat + s_hat);
                z_2 = z2_tmp - gamma_z2 * proj_l2ball(z2_tmp / gamma_z2, x_noisy, epsilon);

                % Compute residual
                res_x = norm(x - x_old, 'fro') / (norm(x_old, 'fro'));
                res_s = norm(s - s_old, 'fro') / (norm(s_old, 'fro'));
                res = max(res_x, res_s);

                if t == 1
                    max_res = res;
                else
                    max_res = max(max_res, res);
                end

                % --- Display the progress ---
                if mod(t, obj.display_freq) == 0 || t == obj.params.pre.max_iter
                    obj.update_status('Preproc.', t, obj.params.pre.max_iter, res, obj.params.pre.criteria, max_res);
                end

                % --- Stopping criteria ---
                if res < obj.params.pre.criteria
                    break;
                end
            end

            x_clean = gather(x);
            s_outliers = gather(s);

            % --- DMD Mode Extraction ---
            [Phi, Lambda, b] = obj.dmd(x_clean);
        end
        
        function [xi, s_outliers] = dimensional_reduction(obj, x_noisy, Phi, Lambda, b, h, v)
            % Algorithm 2 : Dimensionality Reduction
            % Inputs: x_noisy - Noisy observations
            %         Phi - DMD modes
            %         Lambda - DMD eigenvalues
            %         b - Mode amplitudes
            %         h - Number of rows in spatial dimension
            %         v - Number of columns in spatial dimension
            % Outputs: xi - Reduced coefficients
            %          s_outliers - Detected outliers

            if obj.verbose, fprintf('Starting Algorithm 2: Dimensionality Reduction...\n'); end
            % Transport the input data if necessary
            if obj.use_gpu
                x_noisy = gpuArray(x_noisy);
                Phi = gpuArray(Phi);
                Lambda = gpuArray(Lambda);
                b = gpuArray(b);
            end
            [N, M] = size(x_noisy);
            r = size(Phi, 2);
            % Rank modes
            ranking = obj.rank_modes(Phi, Lambda, b, M);
            
            % Initialize variables (Matching Rec_DMD.m)
            xi = zeros(r, 1, 'like', x_noisy); % Corrected initialization (Vector)
            s = zeros(N, M, 'like', x_noisy);
            z_1 = zeros(3*N, M, 'like', x_noisy);
            z_2 = zeros(N, M, 'like', x_noisy);
            C = obj.vand(Lambda, M);

            % Calculate step sizes
            if obj.params.dr.auto_step
                steps = obj.calculate_auto_steps('red', Phi, C);
                gamma_xi = steps.gamma_xi;
                gamma_s = steps.gamma_s;
                gamma_z1 = steps.gamma_z1;
                gamma_z2 = steps.gamma_z2;
            else
                gamma_xi = obj.params.dr.gamma_xi;
                gamma_s = obj.params.dr.gamma_s;
                gamma_z1 = obj.params.dr.gamma_z1;
                gamma_z2 = obj.params.dr.gamma_z2;
            end

            % Hyperparameters
            w = obj.params.dr.w;
            epsilon = obj.params.dr.epsilon;
            eta = obj.params.dr.eta;
            mu = obj.params.dr.mu;

            % Main iteration loop
            if obj.verbose, fprintf('\n'); end
            obj.prev_msg_len = 0; % Reset message length
            max_res = 0;
            for t = 1:obj.params.dr.max_iter
                xi_old = xi;
                s_old = s;

                % Update xi
                grad_xi = sum(conj(C) .* (Phi' * (adj_diff_st(z_1, w, h, v, M) + z_2)), 2);
                xi = prox_l1(xi - gamma_xi * grad_xi, mu .* ranking * gamma_xi);

                % Update s
                s = proj_l1ball(s - gamma_s * z_2, eta);

                % Extrapolation
                xi_hat = 2 * xi - xi_old;
                s_hat = 2 * s - s_old;

                % Update z_1
                z1_tmp = z_1 + gamma_z1 * diff_st(real(Phi * diag(xi_hat) * C), w, h, v);
                z_1 = z1_tmp - gamma_z1 * prox_l12(z1_tmp ./ gamma_z1, 1 ./ gamma_z1);

                % Update z_2
                z2_tmp = z_2 + gamma_z2 * (real(Phi * diag(xi_hat) * C) + s_hat);
                z_2 = z2_tmp - gamma_z2 * proj_l2ball(z2_tmp ./ gamma_z2, x_noisy, epsilon);

                % Compute residual
                res_xi = norm(xi - xi_old) / norm(xi_old);
                res_s = norm(s - s_old) / norm(s_old);
                res = max(res_xi, res_s);

                if t == 1
                    max_res = res;
                else
                    max_res = max(max_res, res);
                end

                % --- Display the progress ---
                if mod(t, obj.display_freq) == 0 || t == obj.params.dr.max_iter
                    obj.update_status('Dim. Red.', t, obj.params.dr.max_iter, gather(res), obj.params.dr.criteria, gather(max_res), 1500);
                end

                % --- Stopping criteria ---
                if t > 1500 && res < obj.params.dr.criteria
                    break;
                end
            end

            xi = gather(xi);
            s_outliers = gather(s);

        end

        function [Phi, Lambda, b] = dmd(obj, x_clean)
            % Standard Mode Extraction (DMD)
            % Inputs: x_clean - Cleaned data
            % Outputs: Phi - DMD modes
            %          Lambda - DMD eigenvalues
            %          b - Mode amplitudes

            if obj.verbose, fprintf('Extracting Modes...\n'); end
            [U, S, V] = svd(x_clean(:, 1:end-1), 'econ');
            if obj.params.me.r > size(S, 1)
                r = size(S, 1);
                warning('CR_DMD:RankReduced', 'Target rank %d is larger than available rank %d. Using %d.', obj.params.me.r, size(S, 1), r);
            else
                r = obj.params.me.r;
            end
            
            U_r = U(:, 1:r);
            S_r = S(1:r, 1:r);
            V_r = V(:, 1:r);
            A_tilde = U_r' * x_clean(:, 2:end) * V_r / S_r;
            [W, D] = eig(A_tilde);
            Phi = x_clean(:, 2:end) * ( V_r * (S_r \ W));
            Lambda = diag(D);
            Vand = obj.vand(Lambda, size(x_clean, 2)-1);
            P = (W' * W) .* (Vand * Vand');
            q = diag(conj(Vand * V_r * S_r' * W));
            b = P \ q;
        end

        function ranking = rank_modes(~, Phi, Lambda, b, M)
            % Rank DMD modes based on energy contribution
            % Inputs: Phi - DMD modes
            %         Lambda - DMD eigenvalues
            %         b - Mode amplitudes
            % Outputs: ranking - Indices of modes sorted by energy contribution

            Lambda_abs = abs(Lambda);
            b_abs = abs(b);
            Phi_abs = sum(abs(Phi.^2), 1)';
            importance = Phi_abs .* b_abs .* ((1 - Lambda_abs .^ M) ./ (1 - Lambda_abs));
            importance = 1 ./ importance;
            ranking = importance / sum(importance);
        end
    end

    methods (Access = private)
        function s1 = merge_structs(obj, s1, s2)
            % Recursively merge s2 into s1. Fields in s2 overwrite s1.
            f = fieldnames(s2);
            for i = 1:length(f)
                if isstruct(s2.(f{i})) && isfield(s1, f{i}) && isstruct(s1.(f{i}))
                    s1.(f{i}) = obj.merge_structs(s1.(f{i}), s2.(f{i}));
                else
                    s1.(f{i}) = s2.(f{i});
                end
            end
        end

        function steps = calculate_auto_steps(obj, stage, varargin)
            % Calculate automatic step sizes based on convergence criteria
            % Inputs: stage - 'pre' or 'dr' indicating the algorithm stage
            %         varargin - Additional parameters as needed
            % Outputs: steps - Struct containing calculated step sizes

            if strcmp(stage, 'pre')
                w = obj.params.pre.w;
                L_norm_sq = 1 + 8*w^2 + 4*(1-w)^2; 
                steps.gamma_x = 1 / L_norm_sq;
                steps.gamma_s = 1;
                steps.gamma_z1 = 0.5;
                steps.gamma_z2 = 0.5;
                
            elseif strcmp(stage, 'red')
                w = obj.params.dr.w;
                Phi = varargin{1};
                C = varargin{2};
                sigma_Phi = norm(Phi); 
                sigma_C = norm(C);
                L_norm_sq = (sigma_Phi^2 * sigma_C^2) * (1 + 8*w^2 + 4*(1-w)^2);
                
                steps.gamma_xi = 1 / L_norm_sq;
                steps.gamma_s = 1;
                steps.gamma_z1 = 0.5;
                steps.gamma_z2 = 0.5;
            end
        end

        function update_status(obj, stage, t, max_t, res, criteria, max_res, min_iter)
            if nargin < 8
                min_iter = 0;
            end
            if ~obj.verbose
                return;
            end
            % Display progress in the console
            len = 15;

            % Double Log-scale progress (log10(log10))
            if max_res > criteria
                % Use log10(abs(log10(x))) to visualize order of magnitude progress
                % Clamp x < 1 (e.g. 0.99) to ensure monotonicity in log-log domain
                get_val = @(x) log10(abs(log10(min(max(x, 1e-20), 0.99))));

                current_val = get_val(res);
                start_val = get_val(max_res);
                target_val = get_val(criteria);

                if abs(start_val - target_val) < 1e-10
                    ratio = 1.0;
                else
                    % Value increases as residual decreases (e.g. 1e-1 -> -1 -> 0, 1e-4 -> -4 -> 0.6)
                    ratio = (current_val - start_val) / (target_val - start_val);
                end
                
                % Clamp ratio between 0 and 1
                ratio = max(0, min(1, ratio));
                pos = floor(ratio * len);
            else
                % Fallback
                pos = 0;
            end

            bar = repmat('=', 1, pos);
            if pos < len
                bar = [bar, '>', repmat(' ', 1, len - pos - 1)];
            end

            % Construct the message
            msg = sprintf('%s [%s] Iter: %5d/%d, Res: %.2e (Goal: %.2e)', ...
                    pad(stage, 12), bar, t, max_t, res, criteria);

            % Clear previous line using backspaces and print new message
            % Use \b to delete previous characters
            reverseStr = repmat('\b', 1, obj.prev_msg_len);
            fprintf(reverseStr);
            fprintf(msg);
            
            % Update previous message length
            obj.prev_msg_len = length(msg);
            
            if t == max_t || (res < criteria && t > min_iter)
                fprintf('\n'); % Print newline at the end
            end
        end

        function C = vand(~, lambda, m)
            % Generate Vandermonde matrix for eigenvalues lambda and m time steps
            C = lambda .^ (0:m-1);
        end
    end
end