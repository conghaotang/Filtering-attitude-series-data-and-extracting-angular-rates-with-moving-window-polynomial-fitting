function [rmse, STD,  Meanz, maxz,minz,rmsez4] = accuracy (alfa) 

% input：
% output：rmse, STD, Mean, max,min 精度评定指标  3*1

%输出均方根、标准差 平均值 最大值 最小值
% length(alfa(1,:))
    rmsez1=sqrt(sum(alfa(1,:).^2)/length(alfa(1,:)));
    rmsez2=sqrt(sum(alfa(2,:).^2)/length(alfa(2,:)));
    rmsez3=sqrt(sum(alfa(3,:).^2)/length(alfa(3,:)));
    rmsez4=sqrt(mean(alfa(1,:).^2 + alfa(2,:).^2 + alfa(3,:).^2));
    rmse=[rmsez1; rmsez2; rmsez3];
    
 %%注意 默认命令std（A），输出的是每一列的标准差，此时除以的是N-1。
 %%如果是std（A,1），此时除以的是N。   
    STD1=std(alfa(1,:),1);  
    STD2=std(alfa(2,:),1);
    STD3=std(alfa(3,:),1);
    STD=[STD1;STD2;STD3];
    
    Mean1=mean(alfa(1,:));
    Mean2=mean(alfa(2,:));
    Mean3=mean(alfa(3,:));
    Meanz=[Mean1;Mean2;Mean3];
    
%     STDz1=sqrt(mean((alfa(1,:)-Mean1).^2));  
%     STDz2=sqrt(mean((alfa(2,:)-Mean2).^2));
%     STDz3=sqrt(mean((alfa(3,:)-Mean3).^2));
%     STDz=[STDz1;STDz2;STDz3];
    
    max1 = max(alfa(1,:));
    max2 = max(alfa(2,:));
    max3 = max(alfa(3,:));
    maxz=[max1;max2;max3];
    
    min1= min(alfa(1,:));
    min2 = min(alfa(2,:));
    min3 = min(alfa(3,:));
    minz=[min1;min2;min3];