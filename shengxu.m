function [AICcz,max_num,max_value]=shengxu(AICcz,N)

 %�������ܣ��� RM1Z ���� RM1Z(2,1,j) ��С��������  %%% ע��Ҫ����ѭ��
 
  N=N-1;    %ҳ��-1
  
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
