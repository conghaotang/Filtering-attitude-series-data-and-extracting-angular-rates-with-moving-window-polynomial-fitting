clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num = xlsread('E:\FIR\DATA20210110C1b.xls');  %注意 去掉 xls中的加载项datam.xls
times=num(:,2)';
qz_1 = num(:,3:6)'; %测量的四元数
mz=length(times);
qz_11 = [times; qz_1]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j=1;
for i=1:10:mz
    qz_t(:,j) = qz_11(:,i);
    j=j+1;
end
timeq= qz_t(1,:);
td = timeq(2) - timeq(1);
qz_ = qz_t(2:5,:);
m = length(timeq);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=15;
P = eye(3);  %Pk=P1; P1=theta2^2*eye(3); [1000 0 0;0 1000 0; 0 0 1000]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=zeros(9,m);   %未知数  
GG=zeros(3,9,m);   %未知数的系数
g_est=zeros(3,m);  %Gibbs 的估计 
qt_est=zeros(5,m);  %真实四元数的估计值
for k=n+1:m   %k需要从 n+1开始 因为需要 历元1作为参考
    
    i=1; % yk 中 g_ 所在的位置
    
    %窗口长度内 观测数据   
    NNz=zeros(9);  % 法矩阵 BTPB
    Wz=zeros(9,1); %  BTPL
    G=zeros(3,9,n+1); %mm=[];
    
    q_kn = qz_(:,k-n);
    qz_KNconj = qconj31(q_kn);
    
    for j= k:-1: (k-n+1)
        p_j = qmul31(qz_(:,j), qz_KNconj);%不需要进行归一化 单位四元数乘单位四元数 仍为单位四元数  
        p_jj = qnormlz(p_j);
        g_j = q2g(p_jj);
        
        Fj=eye(3) + askew(g_j) + g_j * g_j';
        Qj=Fj * P * Fj';
        
        G1=eye(3);  G2=(j-k) * td * eye(3);  G3=((j-k) * td)^2 * eye(3); %令时间间隔为 td s
        Gj=[G1 G2 G3];
                
        NNz = NNz + Gj' * Qj^-1 * Gj;
        Wz = Wz + Gj' * Qj^-1 * g_j;
        i=i+1;
    end

    beta =NNz^-1 * Wz;
    g_est = beta(1:3);  %估计的 Gibbs 矢量
    
    p_est_k = g2q(g_est);
%     gg = 1/sqrt(1 + g_estk' * g_estk);  % 这里的系数 暂时设定为 +1
%     p_est_k = gg * [g_estk; 1]; %相对的 四元数
    q_est_k = qmul31(p_est_k,q_kn);
    qt_est(:,k) = [timeq(1,k);q_est_k];
end

attz = [];
for i=n+1:m
    att1 = q2att31(qt_est(2:end,i));
    att2 = q2att31(qz_(:,i));
    att12 = [qt_est(1,i); att1; att2];
    attz = [attz att12];
end

figureqdata(attz(2:4,:),attz(5:7,:),length(attz),1) %滤波的 原始的
attd = attz(2:4,:)-attz(5:7,:);

Pxx_mea = periodogram(att_mea(1,n+1:end));   
Pxx_est = periodogram(att_est(1,n+1:end));
zzz=[Pxx_mea Pxx_est];


% [Pxx,f]=psd(attz(2,:),Nfft,fs,window,noverlap,dflag);
% [Pxx,f]=psd(x,Nfft,fs,window,noverlap,dflag);
%  matlab中的psd用法:PSD 是做功率谱bai密度的函数du。
% x是信号；zhiNfft快速傅里叶变换点数；fs是采样频dao率；
%window是加的zhuan窗函数；noverlap是指没有重叠率shu（根据采样定理可以算出最小采样频率）；dflag好像是判断前边这个noverlap是否有重叠，有的话就假（不继续），如果没有重叠为真（继续做）。

%   epochs=n+1:m;
%   figure 
%   set(gcf,'color',[1,1,1])
%   subplot(4,1,1);
%   plot(epochs, qt_est(2,n+1:end),'-b','Linewidth',2.5);   %ylabel('Yaw [arcmin]'); %xlabel('epoch')  %,'Linewidth',1
%   hold on
%   plot(epochs,qz_(1,n+1:end),'-r','Linewidth',2.5);
%   set(gca,'Fontsize',10,'Fontname','Times New Roman');
% %   set(gca,'XLim',[0 105]);
% %   set(gca,'XTick',[0 20 40 60 80 100]);
%    
%   subplot(4,1,2);
%   plot(epochs, qt_est(3,n+1:end),'b','Linewidth',2.5); ylabel('Error [rad]'); 
%   hold on
%   plot(epochs,qz_(2,n+1:end),'r','Linewidth',2.5);
%   set(gca,'Fontsize',10,'Fontname','Times New Roman');
% %   set(gca,'XLim',[0 105]);
% %   set(gca,'XTick',[0 20 40 60 80 100]);
%   
%   subplot(4,1,3);
%   plot(epochs, qt_est(4,n+1:end),'b','Linewidth',2.5);  xlabel('Time [s]')  %,'Linewidth',1 ,'Linewidth',1
%   hold on
%   plot(epochs,qz_(3,n+1:end),'r','Linewidth',2.5);
%   
%   subplot(4,1,4);
%   plot(epochs, qt_est(5,n+1:end),'b','Linewidth',2.5);  xlabel('Time [s]')  %,'Linewidth',1 ,'Linewidth',1
%   hold on
%   plot(epochs,qz_(4,n+1:end),'r','Linewidth',2.5);
% %   set(gca,'XLim',[0 105]);
% %   set(gca,'XTick',[0 20 40 60 80 100]);
%   legend('Original','Reconstructed') %intial
%   set(gca,'Fontsize',10,'Fontname','Times New Roman');
  
  
%不存在每一个历元 真实的姿态数据 故无法计算乘性误差四元数
%该方法降低了 高频噪声


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%四元数 到 姿态角
%  figurea(qt_est(2:5,n+1:end),qz_(:,n+1:end),m,n)


