function [RZZ,RZZ_]=data(m,theta1,theta2)
% clc;
% clear;
% randn('state',1e+03);

%   m=100; n=3; theta1=0.01; theta2=0.01;
  RZ=zeros(3,3,m);   RZ(:,:,1)=eye(3);
  R_=zeros(3,3,m);
  
  for i =1:100
      rr=randn(3,1);
  end
  
  % ������ֵ
  betaZ=[];
  for i=2:m
      beta=theta1*randn(3,1); 
      betaZ=[betaZ beta];
      RZ(:,:,i)=(eye(3)+askew(beta))*RZ(:,:,i-1);
  end
  
  
  %�����ɵ� ���� ���� ������ ����  �� ���� ��Ԫ���� ����Ϊ1 ������ �Ծ������У��
  RZZ=zeros(3,3,m); % mm=[]; %attZ=zeros(3,m); attZZ_ = zeros(3,m);
   for k=1:m 
       [qx,max_num]=m2q(RZ(:,:,k));  %���� ת ��Ԫ��
       C = q2mat(qx);   % ��Ԫ�� ת ����
       RZZ(:,:,k) = C;   %У�� �� �� ��ֵ����
   end
    
%    for i=1:102
%        rr=randn(3,1);
%    end
%    
  % ����У����� ��ֵ���� ���ɹ�������
  for i=1:m
      alfa=theta2*randn(3,1);
      R_(:,:,i) = (eye(3)+askew(alfa))^-1 * RZZ(:,:,i);
  end
  
  
 %�Թ������� ����У�� 
  RZZ_=zeros(3,3,m); 
   for k=1:m 
       [q_x,max_num_]=m2q(R_(:,:,k)); 
       C_ = q2mat(q_x);
       RZZ_(:,:,k) = C_;   
   end
   
