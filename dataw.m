function [wqzI,qzt_,etaZ,alfaZ]=dataw(m,theta1w,theta2)
%输入
%m 序列长度
%theta1w 角速率变化的中误差
%theta2  四元数变化的中误差
%输出
%wqzI ： 真实值 时间 角速率 四元数
%qzt_ ： 时间 测量的四元数
%etaZ ：角速率的变化
%alfaZ ：量测的 乘性误差

%生成真实数据 和 量测数据 以四元数形式表示

  wqz = zeros(8,m);
  q0 = [0 0 0 1]';
  w0=[0 0 0]';
  
  wqz(:,1)=[0.001; w0; q0]; %给定初始值
  td = 0.001;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 生成 真值  角速率 序列
  etaZ=[];
  for i = 2 : m*100
      eta = theta1w * randn(3,1); 
      etaZ=[etaZ eta];
      
      wi = w0 + eta;
      err = [wi * td; 1];
      err1 = qnormlz(err);
      qzi = qmul31(err1,q0);
      
      %对四元数进行归一化处理
      qzz = qnormlz(qzi);
      
      w0=wi;       %为下一历元做准备
      q0=qzz;      %为下一历元做准备
      
      wqz(:,i)=[td*i; wi; qzz];  %存储 真值 角速率 【时间 角速率 四元数】   
  end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  %降采样 采样率为0.1s 采样数为200
  for k=1:m
      t = 0.1*k;
      %查找 时间t 在 wqz 中的位置
      wqzt=wqz(1,:)-t;
      [figI,imin] = min(abs(wqzt));  %figI 数值 ,imin 位置
      wqzI(:,k) = wqz(:,imin);%降采样后真实的 时间 角速率和四元数
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % 生成 测量 序列
  alfaZ=[];
  for i=1:m
     alfa=theta2*randn(3,1);
      alfaZ=[alfaZ  alfa]; %乘性误差
      
      qz0_=[alfa;1];
      qzi = wqzI(5:8,i);
      qzz_=qmul31(qz0_,qzi);
      
      %对四元数进行归一化处理
      qzz1_ = qnormlz(qzz_);
      qzt_(:,i)=[wqzI(1,i); qzz1_]; %存储四元数
  end
end