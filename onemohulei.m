% m = size(relation_matrix, 1);
% fuzzy_classes = cell(m, 1);
% 
% for i = 1:m
%     fuzzy_classes{i} = relation_matrix(i, :);
% end
% 
% C = mean(V_norm, 2); 
% % 将 366*1 的矩阵转换为一行矩阵
% C_sum_row = reshape(C', 1, []);
% disp('转换为一行后的矩阵为：');
% disp(C_sum_row);
% 
% % 假设 fuzzy_classes 是 366x366 的数组
% % C_sum_row 是 1x366 的行向量
% 
% % denominator：对 fuzzy_classes 每一行求和，得到 366x1 向量
% denominator = sum(fuzzy_classes, 2);  % 沿行求和
% 
% % numerator：先将 C_sum_row 转置为列向量方便矩阵乘法，或者用点乘实现逐元素乘后求和
% % 这里用矩阵乘法：(fuzzy_classes .* (repmat(C_sum_row, 366, 1))) 逐元素乘法
% % 再按行求和得到 366x1 向量
% 
% numerator = sum(fuzzy_classes .* repmat(C_sum_row, size(fuzzy_classes, 1), 1), 2);
% 
% % 计算概率 Pro，逐元素除
% Pro = zeros(size(denominator));
% nonzero_idx = denominator > 0;
% Pro(nonzero_idx) = numerator(nonzero_idx) ./ denominator(nonzero_idx);
% Pro = zeros(m, 1);  % 结果向量，长度为模糊类数量m
% 
% for i = 1:m
%     numerator = sum(C_sum_row .* fuzzy_classes{i});    % 按元素乘积再求和，分子
%     denominator = sum(fuzzy_classes{i});           % 分母是模糊类向量求和
% 
%     if denominator > 0
%         Pro(i) = numerator / denominator;          % 比值作为概率
%     else
%         Pro(i) = 0;                                % 避免除零
%     end
% end
% 
% disp('计算得到的Pro向量为：');
% disp(Pro);
% 

% 直接用 relation_matrix 就代表 fuzzy_classes，关系矩阵是 m×n 的普通数组
fuzzy_classes = relation_matrix;   % 366 x 366 数组，不用 cell

Ci = mean(V_norm, 2); 
Ci_sum_row = reshape(Ci', 1, []);   % 1 x 366 行向量

% denominator: 每行求和
denominator = sum(fuzzy_classes, 2);    

% numerator: 每行与 C_sum_row 按元素乘积后求和
numerator = sum(fuzzy_classes .* repmat(Ci_sum_row, size(fuzzy_classes, 1), 1), 2);

% 计算概率 Pro
Pro = zeros(size(denominator));
nonzero_idx = denominator > 0;
Pro(nonzero_idx) = numerator(nonzero_idx) ./ denominator(nonzero_idx);

disp('计算得到的Pro向量为：');
disp(Pro);

