
% M = cell(366, 6);   
% for i = 1:366      
%     M{i, 1} = A(i, :);  
%     M{i, 2} = B(i, :);  
%     M{i, 3} = C(i, :);  
%     M{i, 4} = D(i, :);  
%     M{i, 5} = E(i, :);  
%     M{i, 6} = F(i, :);  
% end
% % 初始化
% num_objects = size(M, 1);
% distance_matrix = zeros(num_objects, num_objects, 6);
% 
% for att_idx = 1:6
%     for i = 1:num_objects
%         for j = i+1:num_objects
%           
%             x = cell2mat(M(i, att_idx));
%             y = cell2mat(M(j, att_idx));
%             
%             M_norm = x / (sum(x) + eps);
%             N_norm = y / (sum(y) + eps);
%             
%             
%             max_len = max(length(M_norm), length(N_norm));
%             M_padded = [M_norm, zeros(1, max_len - length(M_norm))];
%             N_padded = [N_norm, zeros(1, max_len - length(N_norm))];
%             
%             % 计算 JS 散度
%             JS_divergence = 0;
%             for k = 1:max_len
%                 m_k = M_padded(k);
%                 n_k = N_padded(k);
%                 if m_k > 0 && n_k > 0
%                     term = m_k * log2(2 * m_k / (m_k + n_k + eps)) + n_k * log2(2 * n_k / (m_k + n_k + eps));
%                 elseif m_k > 0
%                     term = m_k * log2(2);
%                 elseif n_k > 0
%                     term = n_k * log2(2);
%                 else
%                     term = 0;
%                 end
%                 JS_divergence = JS_divergence + term;
%             end
%             
%             % 计算距离并存储
%             distance = sqrt(0.5 * max(JS_divergence, 0));
%             distance_matrix(i, j, att_idx) = real(distance);
%             distance_matrix(j, i, att_idx) = real(distance);
%         end
%     end
% end 
% for att_idx = 1:6
%     fprintf('Distance matrix for attribute %d (first 5 rows and columns):\n', att_idx);
%     disp(distance_matrix(1:5, 1:5, att_idx));
% end
% %得到二元模糊关系矩阵
% relation_matrix = mean(distance_matrix, 3);


fuzzy_classes = relation_matrix;  

Ci = mean(V_norm, 2); 
Ci_sum_row = reshape(Ci', 1, []);   % 1 x 366 行向量

%  每行求和
denominator = sum(fuzzy_classes, 2);    

% numerator: 每行与 C_sum_row 按元素乘积后求和
numerator = sum(fuzzy_classes .* repmat(Ci_sum_row, size(fuzzy_classes, 1), 1), 2);

% 计算概率 Pro
Pro = zeros(size(denominator));
nonzero_idx = denominator > 0;
Pro(nonzero_idx) = numerator(nonzero_idx) ./ denominator(nonzero_idx);

disp('计算得到的Pro向量为：');
disp(Pro);




