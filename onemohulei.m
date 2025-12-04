
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


