% function [RZ,R_,alfaZ]=data(m,theta1,theta2,theta3,n)
clc;
clear;
randn('state',1e+03);

  m=100; n=3; theta1=0.01; theta2=0.01;
  RZ=zeros(3,3,m);   RZ(:,:,1)=eye(3);
  R_=zeros(3,3,m);
  
  % 生成真值
  betaZ=[];
  for i=2:m
      beta=theta1*randn(3,1);
      betaZ=[betaZ beta];
      RZ(:,:,i)=(eye(3)+askew(beta))*RZ(:,:,i-1);
  end
  
  % 生成估计序列
  alfaZ=zeros(3,100);
  for i=1:m
      alfa=theta2*randn(3,1);
      alfaZ(:,i)= alfa;
      R_(:,:,i) = (eye(3)+askew(alfa))^-1 * RZ(:,:,i);
  end
  
  for k=1:m
      Lk2 =  R_(:,:,k) * RZ(:,:,k)';
      Lkk2=  RZ(:,:,k) * R_(:,:,k)';
      alfa_g(:,k)=iaskew(1/2*(Lk2-Lkk2));   
  end
  
  RR_=zeros(3,3,m);
  
  for k=n+1:m
      
      if k==n+1
          
          for j=k:-1:k-n %
              
              S_(:,:,j) = R_(:,:,j)*R_(:,:,k)';
              S_j=S_(:,:,j);
              
              [q_x,maxm]=m2q(S_j);
              Cnb = q2mat(q_x);
              
              RR_(:,:,j) = ((eye(3)+alfa_g(:,j)) * Cnb)'* RZ(:,:,j);
              R_(:,:,j)= RR_(:,:,j);
          end
        
      else
          j=k; %:-1:k-n
          
          S_(:,:,j) = R_(:,:,j)*R_(:,:,k)';
          S_j=S_(:,:,j);
          
          [q_x,maxm]=m2q(S_j);
          Cnb = q2mat(q_x);
          
          RR_(:,:,j) = ((eye(3)+alfa_g(:,j)) * Cnb)'* RZ(:,:,j);
          R_(:,:,j)= RR_(:,:,j);
      end
  end
  
%   save('R_.mat','R_')
%   save('RZ.mat','RZ')
%   for k=1:m
%       Lk2 =  R_(:,:,k) * RZ(:,:,k)';
%       Lkk2=  RZ(:,:,k) * R_(:,:,k)';
%       alfa_g(:,k)=iaskew(1/2*(Lk2-Lkk2));   
%   end
  
  
