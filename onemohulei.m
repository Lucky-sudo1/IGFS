x3=[0.1,0.4,0.5];
y3=[0.2,0.2,0.6];
M3=(x3./(sum(x3)));
N3=(y3./(sum(y3)));
d3=(0.5.*(M3(1,1).*log2(((2.*M3(1,1))./(M3(1,1)+N3(1,1))))+M3(1,2).*log2(((2.*M3(1,2))./(M3(1,2)+N3(1,2))))+M3(1,3).*log2(((2.*M3(1,3))./(M3(1,3)+N3(1,3))))+N3(1,1).*log2(((2.*N3(1,1))./(M3(1,1)+N3(1,1))))+N3(1,2).*log2(((2.*N3(1,2))./(M3(1,2)+N3(1,2))))+N3(1,3).*log2(((2.*N3(1,3))./(M3(1,3)+N3(1,3))))))^0.5
   


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



