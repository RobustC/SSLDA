function SSCM = spatialSignCovariance(X, median)
    % 计算空间符号协方差矩阵（SSCM）
    % X: 数据矩阵，每行是一个数据点
    % median: 已知的空间中位数

    n = size(X, 1);   % 数据点数量
    d = size(X, 2);   % 数据维度

    % 初始化 SSCM 矩阵
    SSCM = zeros(d, d);
    
    % 计算每个数据点的符号向量，并累加到 SSCM
    for i = 1:n
        diff = X(i, :) - median; % 数据点相对于中位数的差值
        normDiff = norm(diff);   % 差值的范数
        if normDiff ~= 0
            u = diff / normDiff; % 符号向量
            SSCM = SSCM + (u' * u); % 更新 SSCM 矩阵
        end
    end

    % 计算最终 SSCM（取均值）
    SSCM = SSCM / n;
end
