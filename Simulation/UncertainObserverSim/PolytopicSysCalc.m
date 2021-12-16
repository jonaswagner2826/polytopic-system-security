function A = PolytopicSysCalc(A_all, alpha)
    %POLYTOPICSYSCALC Calculates A(\alpha) from A_all
    arguments
        A_all 
        alpha (:,1) double
    end
    
    if size(A_all, 3) == 1
        A = A_all;
        return
    elseif size(alpha,1) ~= size(A_all,3)
        error('A_all not compatable with alpha (m \neq m)')
    end
    
    A = zeros(size(A_all,1,2));
    for i = 1:size(alpha,1)
        A = A + alpha(i) * A_all(:,:,i);
    end
end

