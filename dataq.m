function [qz,qz_,etaZ,alfaZ]=dataq(m,theta1,theta2)

%������ʵ���� �� �������� ����Ԫ����ʽ��ʾ

  qz = zeros(4,m);
  q0 = [0 0 0 1]';
  
  qz(:,1)=q0; %������ʼֵ
  
  % ���� ��ֵ ����
  etaZ=[];
  for i=2:m
      eta=theta1*randn(3,1); 
      etaZ=[etaZ eta]; %�������

      qz0=[eta;1];
      qzz=qmul31(qz0,q0);
      
      %����Ԫ�����й�һ������
      qzz1 = qnormlz(qzz);
      
      q0=qzz1;      %Ϊ��һ��Ԫ��׼��
      qz(:,i)=qzz1; %�洢��Ԫ�� 
  end
   
  % ���� ���� ����
  alfaZ=[];
  for i=1:m
      alfa=theta2*randn(3,1);
      alfaZ=[alfaZ  alfa]; %�������
      
      qz0_=[alfa;1];
      qzz_=qmul31(qz0_,qz(:,i));
      
      %����Ԫ�����й�һ������
      qzz1_ = qnormlz(qzz_);
      qz_(:,i)=qzz1_; %�洢��Ԫ��
  end
  
end