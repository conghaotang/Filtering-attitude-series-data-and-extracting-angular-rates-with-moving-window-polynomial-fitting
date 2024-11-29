function q = g2q(g)  

% Gibbs 转换为 四元数

   gg = 1/sqrt(1 + g' * g);  % 这里的系数 暂时设定为 +1
   q = gg * [g; 1]; %相对的 四元数
   
end