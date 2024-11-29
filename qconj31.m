function qo = qconj31(qi)
% Quaternion conjugation.
% 
% Prototype: qo = qconj(qi)
% Input: qi - input quaternion
% Output: qo - output quaternion ,if qi = [qi(1); qi(2:4)]
%              then qo = [qi(1); -qi(2:4)]

%四元数的共轭 （四元数的表示：前三为矢量 第四为 标量）

 qo = [-qi(1:3); qi(4)];