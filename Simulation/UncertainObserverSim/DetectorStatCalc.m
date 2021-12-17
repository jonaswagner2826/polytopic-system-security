function Z = DetectorStatCalc(results, sim_params)
    
    % DetectorStatCalc - Determines z_k prediction...
    %results is the output for a PolytopicSysSim...
    
%     % Parameters
%     z_k = @(X_0, M, idx_x_0, idx_k) X_0{idx_x_0}' * M{idx_k} * X_0{idx_x_0};

    M = DetectorStatCalcMatrix(results,sim_params);
    
    M = M(:);
    X_0 = {results.x_0}';
    
    
    Z = zeros(numel(X_0), sim_params.N);
    for idx_x_0 = 1:numel(X_0)
        for idx_k = 1:sim_params.N
            Z(idx_x_0, idx_k) = ...
                X_0{idx_x_0}' * M{idx_x_0}{idx_k} * X_0{idx_x_0};
        end
    end
    
    Z = reshape(num2cell(Z,2), size(results));
%     
%     Z = reshape(Z,[size(results))
    
    
%     Z = arrayfun(@(idx_x_0, idx_k) z_k(X_0, M, idx_x_0, idx_k),...
%         1:numel(X_0), 1:sim_params.N)
    
    
%     Z = cell(size(M));
%     Z = cell(sim_params.N,1);
%     for k = 1:sim_params.N
%         Z = arrayfun(@(k) ...
%             cellfun(@(x_0) ...
%             cellfun(@(M) z_k(M, x_0, k),...
%             M, 'UniformOutput', false),...
%             X_0, 'UniformOutput', false),...
%             1:sim_params.N, 'UniformOutput', false);
    % end
    
%     Z1 = cell2mat([Z{:,1}])
    
    
    
end