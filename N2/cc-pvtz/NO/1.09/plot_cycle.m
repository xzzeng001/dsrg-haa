% 读取one-body数据
fid = fopen('one-body.txt', 'r');
data_one = textscan(fid, '%d%f%f', 'Delimiter', ',');
fclose(fid);
ids_one = data_one{1};
values_one = data_one{3};
N = length(ids_one);

% 生成均匀分布的角度
theta = linspace(0, 2*pi, N+1)';
theta(end) = [];
R = 1; % 大圆半径
x = R * cos(theta);
y = R * sin(theta);

% 创建ID到索引的映射
id_map = containers.Map(ids_one, 1:N);

% 归一化点大小（范围20-100）
abs_values_one = abs(values_one);
max_size = max(abs_values_one);
min_size = min(abs_values_one);
if max_size == min_size
    sizes = 50 * ones(N,1);
else
    sizes = ((abs_values_one - min_size) / (max_size - min_size)) * 600 + 20;
end

% 设置小圆的颜色（选择蓝色调或绿色调，醒目且统一）
circle_color = [0.2, 0.6, 0.2];  % 绿色调，显眼且统一 

% 绘制散点图
figure;
hold on;

% 绘制小圆
scatter(x, y, sizes, 'filled', 'MarkerFaceColor', circle_color);

% 标注编号（沿径向偏移）
text_offset = 0.07;
for i = 1:N
    angle = theta(i);
    dx = text_offset * cos(angle);
    dy = text_offset * sin(angle);
    text(x(i)+dx, y(i)+dy, num2str(ids_one(i)), ...
        'HorizontalAlignment', 'center', 'FontSize', 20, 'Color', 'k');
end

% 读取two-body数据
fid = fopen('two-body.txt', 'r');
data_two = textscan(fid, '%d%d%f%f', 'Delimiter', ',');
fclose(fid);
id1 = data_two{1};
id2 = data_two{2};
line_values = abs(data_two{4});

% 选择前K大的连接线
K = 40;  % 可以调整K的值
[~, sorted_indices] = sort(line_values, 'descend');
top_k_indices = sorted_indices(1:K);

% 归一化线宽（范围1-5）
max_lw = max(line_values(top_k_indices));
min_lw = min(line_values(top_k_indices));
if max_lw == min_lw
    line_widths = 2 * ones(K,1);
else
    line_widths = ((line_values(top_k_indices) - min_lw) / (max_lw - min_lw)) * 4 + 1;
end

% 线条颜色由浅灰色到深灰色渐变
line_colors = linspace(0.2, 0.9, K);  % 从浅灰色(0.6)到深灰色(0.2)

% 绘制连接线
for k = 1:K
    idx = top_k_indices(k);
    if isKey(id_map, id1(idx)) && isKey(id_map, id2(idx))
        idx1 = id_map(id1(idx));
        idx2 = id_map(id2(idx));
        x1 = x(idx1); y1 = y(idx1);
        x2 = x(idx2); y2 = y(idx2);
        
        % 计算渐变的灰色
        line_color = [line_colors(k) line_colors(k) line_colors(k)];
        
        % 绘制线条
        plot([x1, x2], [y1, y2], 'Color', line_color, 'LineWidth', line_widths(k));
    end
end

% 绘制大圆（外圆）
theta_large = linspace(0, 2*pi, 100);
x_large = 1.2 * R * cos(theta_large);  % 增大大圆半径（1.2倍小圆半径）
y_large = 1.2 * R * sin(theta_large);

plot(x_large, y_large, 'k-', 'LineWidth', 2);  % 大圆边框

% 图形修饰
axis equal;
axis off;
%title('Network Visualization with Large Circle');
hold off;
