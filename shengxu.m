function [AICcz,max_num,max_value]=shengxu(AICcz,N)

 %函数功能：将 RM1Z 依据 RM1Z(2,1,j) 从小到大排列  %%% 注意要两次循环
 
  N=N-1;    %页数-1
  
  for i=0:N
      for j=1:N-i
          if AICcz(1,j) > AICcz(1,j+1)
              x = AICcz(:,j);
              AICcz(:,j) = AICcz(:,j+1);
              AICcz(:,j+1) = x;
          end
      end
  end
  
  max_num=AICcz(2,4);
  max_value = AICcz(1,4);
