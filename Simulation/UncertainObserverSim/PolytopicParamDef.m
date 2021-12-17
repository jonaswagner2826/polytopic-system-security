function [ALPHA_real, ALPHA_hat] = PolytopicParamDef(alphaSelect, m)
    arguments
        alphaSelect
        m = 4% Assumes m = 4...
    end
    
    if alphaSelect == 0
        ALPHA_real = [1,0,0,0];%[0.25, 0.25, 0.25, 0.25];
        % Alpha_hat = [0,0,0,1];%[0.375, 0.125, 0.125, 0.125];
        ALPHA_hat = [0.999,0.001,0,0];
        
    elseif alphaSelect == 1
        %Random Alphas
        num_extra_alpha_real = 0;%500;
        num_extra_alpha_hat = 3;
        ALPHA_real = [eye(m), normalize(ones(m,1),1,'norm',1),...
            normalize(rand(m, num_extra_alpha_real),1,'norm',1)];
        ALPHA_hat = [(1/m)*ones(m,1),...[eye(m), [0.25;0.25;0.25;0.25],...
            normalize(rand(m, num_extra_alpha_hat),1,'norm',1)];
        
    elseif alphaSelect == 2
        % Random Edge Alphas
        num_extra_alpha_real = 5; % per subsys edge
        num_extra_alpha_hat = 0;
        ALPHA_real = [eye(m), (1/m) * ones(m,1)];
        ALPHA_hat = [(1/m) * ones(m,1), eye(m)]; %
        
        for idx = nchoosek(1:m,2)'
            for k = 1:num_extra_alpha_real
                temp = zeros(m,1);
                temp(idx(1)) = k;
                temp(idx(2)) = (num_extra_alpha_real - k);
                temp = normalize(temp,1,'norm',1);
                ALPHA_real = [ALPHA_real, temp];
            end
            for k = 1:num_extra_alpha_hat
                temp = zeros(m,1);
                temp(idx(1)) = k;
                temp(idx(2)) = num_extra_alpha_hat - k;
                temp = normalize(temp,1,'norm',1);
                ALPHA_hat = [ALPHA_hat, temp];
            end
        end
    end
end