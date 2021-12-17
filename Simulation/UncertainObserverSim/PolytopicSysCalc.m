function sys_all = PolytopicSysCalc(sys_poly, alpha_all)
    %POLYTOPICSYSCALC Calculates A(\alpha) from A_all
    arguments
        sys_poly
        alpha_all
    end
    
    m = size(alpha_all,1);
    
    M = struct2cell(sys_poly);
    for i = 1:4 %A,B,C,D
        if size(M{i},3) == m
            continue
        elseif size(M{i},3) == 1
            M{i} = M{i} .* ones(1,1,m);
        else
            error('issue with input... m \neq m');
        end        
    end
    
    sys_all = cell(size(alpha_all,1),1).';
    for i = 1:size(alpha_all,2)
        sys_all{i}.A = PolytopicSysMatrixCalc(M{1}, alpha_all(:,i));
        sys_all{i}.B = PolytopicSysMatrixCalc(M{2}, alpha_all(:,i));
        sys_all{i}.C = PolytopicSysMatrixCalc(M{3}, alpha_all(:,i));
        sys_all{i}.D = PolytopicSysMatrixCalc(M{4}, alpha_all(:,i));
    end
    
    % Account for single Alphas
    if size(alpha_all,2)==1, sys_all = sys_all{:,1}; end
    
%         
%     
%     
%     % do the thing with calling individual matrix calc on it...
%     
%     if size(A_all, 3) == 1
%         A = A_all;
%         return
%     elseif size(alpha,1) ~= size(A_all,3)
%         error('A_all not compatable with alpha (m \neq m)')
% %     elseif size(alpha,2) > 1
% %         error(['Calc designed to run only on a single set of alphas...',...
% %                 'Run as a cell function for alpha instead'])
%     end
%     
%     if size(alpha,2) == 1
%         A = zeros(size(A_all,1,2));
%         for i = 1:size(alpha,1)
%             A = A + alpha(i) * A_all(:,:,i);
%         end
%     else
%         A = cell(1,size(alpha,2));
%         for j = 1:size(alpha,2)
%             A{j} = zeros(size(A_all,1,2));
%             for i = 1:size(alpha,1)
%                 A{j} = A{j} + alpha(i,j) * A_all(:,:,i);
%             end
%         end
%     end
            
    
end

