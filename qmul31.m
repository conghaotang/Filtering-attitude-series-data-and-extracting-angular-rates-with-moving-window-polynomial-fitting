function q = qmul31(q1, q2)
% Quaternion multiplication: q = q1*q2.
% 
% Prototype: q = qmul(q1, q2)
% Inputs: q1, q2 - input quaternion
% Output: q - output quaternion ,such that q = q1*q2

%四元数乘法 （四元数的表示：前三为矢量 第四为 标量）

  v_=q1(1:3);   s_=q1(4);
  v =q2(1:3);   s =q2(4);   
  
  q3 = s * v_ + s_ * v - askew(v_) * v;
  q4 = s * s_ - v_' * v;
  
  q=[q3; q4];