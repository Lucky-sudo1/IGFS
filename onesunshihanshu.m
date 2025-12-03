% 假设已经存在名为M的366行6列的单元数组，每个单元格直接存储3个数值的数组  
% 这里不进行M的初始化，实际使用时请确保M已正确赋值  
% 将六个矩阵的每一行合并到 M 中  
for i = 1:366  
    % 将每个矩阵的第 i 行赋值给 M  
    M{i, 1} = A(i, :);  
    M{i, 2} = B(i, :);  
    M{i, 3} = C(i, :);  
    M{i, 4} = D(i, :);  
    M{i, 5} = E(i, :);  
    M{i, 6} = F(i, :);  
end
% 获取M的行数和列数  
numRows = size(M, 1);  
numCols = size(M, 2);  

% 初始化结果单元数组MT  
MT = cell(numRows, numCols);  

% 遍历M的每一个元素  
for i = 1:numRows  
    for j = 1:numCols  
        % 获取当前属性下的数值数组  
        currentArray = M{i, j}; % 直接访问数值数组  
        
        % 确保 currentArray 是 1x3 数值数组  
        if isnumeric(currentArray) && isequal(size(currentArray), [1, 3])  
            % 对数值数组进行从大到小排序  
            sortedArray = sort(currentArray, 'descend');  
            
            % 将排序后的数值数组转换回单元数组形式  
            sortedCell = num2cell(sortedArray);  
            
            % 将排序后的单元数组保存到MT中  
            MT{i, j} = sortedCell;  
        else  
            error('M{%d, %d} 的内容不符合预期，期望为一个长度为3的数值数组，而实际内容为：%s', i, j, mat2str(currentArray));  
        end  
    end  
end

% 获取MT的行数和列数
numRows = size(MT, 1);
numCols = size(MT, 2);

% 初始化正理想点和负理想点的单元数组
h_plus = cell(1, numCols);
h_minus = cell(1, numCols);

% 遍历MT的每一列
for j = 1:numCols
    first_elements = zeros(numRows, 1);
    second_elements = zeros(numRows, 1);
    third_elements = zeros(numRows, 1);
    
    % 提取每一列中每个犹豫模糊集的三个元素
    for i = 1:numRows
        currentCell = MT{i, j};
        first_elements(i) = currentCell{1};
        second_elements(i) = currentCell{2};
        third_elements(i) = currentCell{3};
    end
    
    % 计算正理想点
    max_first = max(first_elements);
    max_second = max(second_elements);
    max_third = max(third_elements);
    h_plus{j} = num2cell([max_first, max_second, max_third]);
    
    % 计算负理想点
    min_first = min(first_elements);
    min_second = min(second_elements);
    min_third = min(third_elements);
    h_minus{j} = num2cell([min_first, min_second, min_third]);
end

% 初始化距离矩阵
distance_to_zheng = zeros(numRows, numCols);
distance_to_fu = zeros(numRows, numCols);

% 计算每个元素到正理想点和负理想点的距离
for att_idx = 1:numCols
    for i = 1:numRows
        % 获取当前元素
        x = cell2mat(MT{i, att_idx}); % 修改为花括号访问内部元胞数组
        
        % 获取正理想点
        plus_array = cell2mat(h_plus{att_idx});
        % 归一化（避免除零）
        M_norm = x / (sum(x) + eps);
        plus_norm = plus_array / (sum(plus_array) + eps);
        
        % 对齐长度
        max_len = max(length(M_norm), length(plus_norm));
        M_padded = [M_norm, zeros(1, max_len - length(M_norm))];
        plus_padded = [plus_norm, zeros(1, max_len - length(plus_norm))];
        
        % 计算 JS 散度
        JS_divergence_zheng = 0;
        for k = 1:max_len
            m_k = M_padded(k);
            p_k = plus_padded(k);
            if m_k > 0 && p_k > 0
                term = m_k * log2(2 * m_k / (m_k + p_k + eps)) + p_k * log2(2 * p_k / (m_k + p_k + eps));
            elseif m_k > 0
                term = m_k * log2(2);
            elseif p_k > 0
                term = p_k * log2(2);
            else
                term = 0;
            end
            JS_divergence_zheng= JS_divergence_zheng + term;
        end
        
        % 计算距离并存储
        distance_zheng = sqrt(0.5 * max(JS_divergence_zheng, 0));
        distance_to_zheng(i, att_idx) = real(distance_zheng);
        
        % 获取负理想点
        minus_array = cell2mat(h_minus{att_idx});
        % 归一化（避免除零）
        minus_norm = minus_array / (sum(minus_array) + eps);
        
        % 对齐长度
        max_len = max(length(M_norm), length(minus_norm));
        minus_padded = [minus_norm, zeros(1, max_len - length(minus_norm))];
        
        % 计算 JS 散度
        JS_divergence_fu = 0;
        for k = 1:max_len
            m_k = M_padded(k);
            n_k = minus_padded(k);
            if m_k > 0 && n_k > 0
                term = m_k * log2(2 * m_k / (m_k + n_k + eps)) + n_k * log2(2 * n_k / (m_k + n_k + eps));
            elseif m_k > 0
                term = m_k * log2(2);
            elseif n_k > 0
                term = n_k * log2(2);
            else
                term = 0;
            end
            JS_divergence_fu = JS_divergence_fu + term;
        end
        
        % 计算距离并存储
        distance_fu = sqrt(0.5 * max(JS_divergence_fu, 0));
        distance_to_fu(i, att_idx) = real(distance_fu);
    end
end

% Parameters  
m = 6; % Number of alternatives (c1 to c8)  
n = 366; % Number of attributes (a1 to a4)  
xi = 0.88; % Parameter for loss calculation  
v=2.25;
rho = 0.88; % Parameter for gain calculation   
%  omega=[0.05,0.3,0.2,0.05,0.1,0.3]%dfenleiquanzhong
 omega=[0.1666,0.1684, 0.1661,0.1675,0.1661,0.1653]
% 计算损失 I^- 和收益 I^+  
% 计算损失 I^- 和收益 I^+  
I_fu = -sum(v.* distance_to_zheng.^ xi .* repmat( omega, numRows, 1), 2);   
I_zheng = sum(distance_to_fu .^ rho .* repmat( omega, numRows, 1), 2);  

%计算总指标 I(z_t)  
Si =0.5.*(I_fu+I_zheng);

yi=Si;

 % 找出chi里面的最大和最小值
    yi_max = max(Si);
    yi_min = min(Si);


% 设置阈值 s
%    s = 0.43; %diyipian分类
  s=0.3%SCC


    num_alternatives = length(yi); % 假设备选方案数量等于chi的长度

    % 初始化相对效用值矩阵，num_alternatives行6列，分别对应6个相对效用函数
    relative_utility_values = zeros(num_alternatives, 6);

    for i = 1:num_alternatives
        % 对于hPP函数
        relative_utility_values(i, 1) = 0;
        % 对于hPN函数
        relative_utility_values(i, 2) = yi_max - yi(i);
        % 对于hBP函数
        relative_utility_values(i, 3) = s * (yi(i) - yi_min);
        % 对于hBN函数
        relative_utility_values(i, 4) = s * (yi_max - yi(i));
        % 对于hNP函数
        relative_utility_values(i, 5) = yi(i) - yi_min;
        % 对于hNN函数
        relative_utility_values(i, 6) = 0;
    end

    disp('每个备选方案的6个相对效用值:');
    disp(relative_utility_values);
    %将每个对象的六个效用值赋值给U
    U = relative_utility_values;
    % 计算每个对象的期望效用
    expected_utilities = zeros(num_alternatives, 3); % 初始化期望效用矩阵

    % 将good_state_probs赋值给probabilities
    probabilities = Pro; 

    % 遍历每个对象
    for i = 1:num_alternatives
        % 对于每个对象，根据对应的概率集计算期望效用
        % 采取 h_p 行动的期望效用
        expected_utilities(i, 1) = probabilities(i) * U(i, 1) + (1 - probabilities(i)) * U(i, 2);
        
        % 采取 h_b 行动的期望效用
        expected_utilities(i, 2) = probabilities(i) * U(i, 3) + (1 - probabilities(i)) * U(i, 4);
        
        % 采取 h_n 行动的期望效用
        expected_utilities(i, 3) = probabilities(i) * U(i, 5) + (1 - probabilities(i)) * U(i, 6);
    end

    % 显示每个对象的期望效用
    disp('每个对象的期望损失：');
    disp(expected_utilities);

    classification_result = cell(num_alternatives, 1);

    % 进行分类
    for i = 1:num_alternatives
        if expected_utilities(i, 1) < expected_utilities(i, 2) && expected_utilities(i, 1) <expected_utilities(i, 3)
            classification_result{i} = 'POS(X)';
        elseif expected_utilities(i, 2) < expected_utilities(i, 1) && expected_utilities(i, 2)< expected_utilities(i, 3)
            classification_result{i} = 'BND(X)';
        elseif expected_utilities(i, 3) < expected_utilities(i, 1) && expected_utilities(i, 3) < expected_utilities(i, 2)
            classification_result{i} = 'NEG(X)';
        end
    end

    % 显示分类结果
    disp('每个备选方案的分类结果：');
    for i = 1:num_alternatives
        fprintf('备选方案 %d: %s\n', i, classification_result{i});
    end

    % 对每个分类内的替代方案根据期望效用值进行排序
    pos_indices = find(strcmp(classification_result, 'POS(X)'));
    [~, sort_order_pos] = sort(expected_utilities(pos_indices, 1), 'ascend');
    bnd_indices = find(strcmp(classification_result, 'BND(X)'));
    [~, sort_order_bnd] = sort(expected_utilities(bnd_indices, 2), 'ascend');
    neg_indices = find(strcmp(classification_result, 'NEG(X)'));
    [~, sort_order_neg] = sort(expected_utilities(neg_indices, 3), 'ascend');
    
  % 输出分类集合结果
    pos_set = pos_indices(sort_order_pos);
    bnd_set = bnd_indices(sort_order_bnd);
    neg_set = neg_indices(sort_order_neg);
    
    % 计算各分类的基数
    pos_cardinality = length(pos_set);
    bnd_cardinality = length(bnd_set);
    neg_cardinality = length(neg_set);
    
    fprintf('POS(X)={');
    fprintf('%d,', pos_set);
    fprintf('}, 基数: %d\n', pos_cardinality);
    fprintf('BND(X)={');
    fprintf('%d,', bnd_set);
    fprintf('}, 基数: %d\n', bnd_cardinality);
    fprintf('NEG(X)={');
    fprintf('%d,', neg_set);
    fprintf('}, 基数: %d\n', neg_cardinality);

    
   
    % 合并排序结果
    sorted_indices = [pos_indices(sort_order_pos); bnd_indices(sort_order_bnd); neg_indices(sort_order_neg)];
    % 显示合并后的排序结果
    disp('合并后的排序结果：');
    disp(sorted_indices);
   % 创建一个数组来记录每个备选方案的位置
position_array = zeros(num_alternatives, 1);
for i = 1:length(sorted_indices)
    position_array(sorted_indices(i)) = i;
end

% 显示每个备选方案的位置
disp('每个备选方案的位置：');
disp(position_array);
