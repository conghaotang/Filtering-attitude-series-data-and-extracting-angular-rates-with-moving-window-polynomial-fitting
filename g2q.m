function q = g2q(g)  

% Gibbs ת��Ϊ ��Ԫ��

   gg = 1/sqrt(1 + g' * g);  % �����ϵ�� ��ʱ�趨Ϊ +1
   q = gg * [g; 1]; %��Ե� ��Ԫ��
   
end