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
  
  % 生成真值
  betaZ=[];
  for i=2:m
      beta=theta1*randn(3,1); 
      betaZ=[betaZ beta];
      RZ(:,:,i)=(eye(3)+askew(beta))*RZ(:,:,i-1);
  end
  
  
  %对生成的 矩阵 进行 正交化 处理  即 利用 四元数的 范数为1 的特性 对矩阵进行校正
  RZZ=zeros(3,3,m); % mm=[]; %attZ=zeros(3,m); attZZ_ = zeros(3,m);
   for k=1:m 
       [qx,max_num]=m2q(RZ(:,:,k));  %矩阵 转 四元数
       C = q2mat(qx);   % 四元数 转 矩阵
       RZZ(:,:,k) = C;   %校正 后 的 真值矩阵
   end
    
%    for i=1:102
%        rr=randn(3,1);
%    end
%    
  % 利用校正后的 真值矩阵 生成估计序列
  for i=1:m
      alfa=theta2*randn(3,1);
      R_(:,:,i) = (eye(3)+askew(alfa))^-1 * RZZ(:,:,i);
  end
  
  
 %对估计序列 进行校正 
  RZZ_=zeros(3,3,m); 
   for k=1:m 
       [q_x,max_num_]=m2q(R_(:,:,k)); 
       C_ = q2mat(q_x);
       RZZ_(:,:,k) = C_;   
   end
   
