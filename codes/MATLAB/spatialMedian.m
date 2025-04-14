
%v = 1:10; % 1到6的向量
%X = reshape(v, [5,2])
function median = spatialMedian(X, tol, maxIter)
    % 计算空间中位数（Spatial Median）
    % X: 数据矩阵，每行是一个数据点
    % tol: 迭代的容差（默认 1e-6）
    % maxIter: 最大迭代次数（默认 100）

    if nargin < 2
        tol = 1e-6;
    end
    if nargin < 3
        maxIter = 100;
    end

    % 初始中位数位置：使用数据中心的均值作为初始值
    median = mean(X, 1);
    n = size(X, 1);

    for iter = 1:maxIter
        % 计算到每个点的距离
        distances = sqrt(sum((X - median).^2, 2));

        % 检查是否有距离为零的点
        zeroDist = (distances == 0);
        if any(zeroDist)
            % 如果有距离为零的点，直接返回该点为中位数
            median = X(zeroDist, :);
            return;
        end

        % 更新中位数位置
        weights = 1 ./ distances;
        newMedian = sum(weights .* X) / sum(weights);

        % 检查收敛条件
        if norm(newMedian - median) < tol
            median = newMedian;
            return;
        end

        % 更新当前中位数位置
        median = newMedian;
    end

    warning('Spatial median did not converge within the maximum number of iterations.');
end
