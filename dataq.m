function [qz,qz_,etaZ,alfaZ]=dataq(m,theta1,theta2)

%生成真实数据 和 量测数据 以四元数形式表示

  qz = zeros(4,m);
  q0 = [0 0 0 1]';
  
  qz(:,1)=q0; %给定初始值
  
  % 生成 真值 序列
  etaZ=[];
  for i=2:m
      eta=theta1*randn(3,1); 
      etaZ=[etaZ eta]; %乘性误差

      qz0=[eta;1];
      qzz=qmul31(qz0,q0);
      
      %对四元数进行归一化处理
      qzz1 = qnormlz(qzz);
      
      q0=qzz1;      %为下一历元做准备
      qz(:,i)=qzz1; %存储四元数 
  end
   
  % 生成 测量 序列
  alfaZ=[];
  for i=1:m
      alfa=theta2*randn(3,1);
      alfaZ=[alfaZ  alfa]; %乘性误差
      
      qz0_=[alfa;1];
      qzz_=qmul31(qz0_,qz(:,i));
      
      %对四元数进行归一化处理
      qzz1_ = qnormlz(qzz_);
      qz_(:,i)=qzz1_; %存储四元数
  end
  
end