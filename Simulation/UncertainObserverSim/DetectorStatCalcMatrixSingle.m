function M = DetectorStatCalcMatrixSingle(sys_real,sys_hat, L, sim_params)
    
    % DetectorStatCalcMatrixSingle - Determines M_ for z_k prediction...
    %fed from the other function...
    
    % Parameters
    N = sim_params.N;
    SigmaInv = sim_params.SigmaInv;
    
    % System
    A = sys_real.A;
    A_hat = sys_hat.A;
    C = sys_real.C;
    C_hat = sys_hat.C;
    
    M = cell(N,1);
    for k = 1:N
        M{k} = zeros(size(A,1));
        for i = (1:k)-1
            for j = (1:k)-1
                M{k} = M{k} ...
                    + (A^(k-j))' * (A - A_hat)' * (A_hat - L*C_hat)^(j) ...
                        * (C' * SigmaInv * C)...
                        * ((A_hat - L*C_hat)^(i))' * (A-A_hat) * (A^(k-i));
            end
        end
    end
end