function qo = qconj31(qi)
% Quaternion conjugation.
% 
% Prototype: qo = qconj(qi)
% Input: qi - input quaternion
% Output: qo - output quaternion ,if qi = [qi(1); qi(2:4)]
%              then qo = [qi(1); -qi(2:4)]

%��Ԫ���Ĺ��� ����Ԫ���ı�ʾ��ǰ��Ϊʸ�� ����Ϊ ������

 qo = [-qi(1:3); qi(4)];