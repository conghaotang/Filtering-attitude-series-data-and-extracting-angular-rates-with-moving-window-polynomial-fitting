clear all;
close all;
clc;
clo = 996;  % 1994 2000 1e+03;  1e+03; %1000;sum(100*clock) sum(10*clock)
randn('state',clo); %0  1e+03     2.1041705 2.10414e+05
%%%%%qwcase1 clo = 1000;  n=16;  theta1=0.001;   theta2=0.05;
%%%%%qwcase2 clo = 996;  n=14;  theta1=0.0001;   theta2=0.01;

m= 200;  %序列总长度
n=14;  %窗口长度 
theta1=0.0001;   theta2=0.01; %theta3=0.05;
[wqzI,qzt_,etawZ,alfaZ]=dataw(m,theta1,theta2);
td = qzt_(1,2) - qzt_(1,1);
qz  = wqzI(5:8,:);
qz_ = qzt_(2:5,:);
P=theta2^2*eye(3);  %Pk=P1; P1=theta2^2*eye(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta=zeros(9,m);   %未知数  
GG=zeros(3,9,m);   %未知数的系数
g_est=zeros(3,m);  %Gibbs 的估计 
q_est=zeros(4,m);  %真实四元数的估计值
wqkz = zeros(8,m); 

for k=n+1:m   %k需要从 n+1开始 因为需要 历元1作为参考
    
    i=1; % yk 中 g_ 所在的位置
    
    %窗口长度内 观测数据   
    NNz=zeros(9);  % 法矩阵 BTPB
    Wz=zeros(9,1); %  BTPL
    G=zeros(3,9,n+1); %mm=[];
    
    q_kn = qz_(:,k-n);
    h_kn = q2g(q_kn);
    qz_KNconj = qconj31(q_kn);
   
    for j= k:-1: (k-n+1) 
        p_j  = qmul31(qz_(:,j), qz_KNconj);   
        p_jj = qnormlz(p_j);  % 归一化
        g_j = q2g(p_jj);
        
        Fj=eye(3) + askew(g_j) + g_j * g_j';
        Qj=Fj * P * Fj';
        
        G1=eye(3);  G2=(j-k)*td*eye(3);  G3=((j-k)*td)^2 * eye(3); %令时间间隔为 1 s
        Gj=[G1 G2 G3];
       
        NNz = NNz + Gj' * Qj^-1 * Gj;
        Wz = Wz + Gj' * Qj^-1 * g_j;
        i=i+1;
    end

    beta=NNz^-1 * Wz;
    g_est =  beta(1:3);  %估计的 Gibbs 矢量  (:,k):
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %提取 角速率
    hk_est = gmul(g_est,h_kn);
    g_est_f = beta(4:6);   % Gibbs的微分
%     hk_est_f = gmul(g_est_f,h_kn); % h估计的微分
    hk_est_f = gmulf(g_est,h_kn,g_est_f);
    wkfz = 2/(1 + hk_est' * hk_est);
    wkfm = (hk_est_f - askew(hk_est) * hk_est_f);
    wk_est = wkfz * wkfm ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %计算 四元数
    p_est_k = g2q(g_est);
    q_est_k = qmul31(p_est_k,q_kn);
    
    wqkz(:,k)= [qzt_(1,k); wk_est; q_est_k];  %估计值 时间 角速率 四元数 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%对比 方法 计算 角速率
qz_est = wqkz([1,5:8],:);
wdz=zeros(4,m); qdjz=[];
for k=n+2:m 
    qdq = (qz_est(2:5,k) - qz_est(2:5,k-1))/td;
    qdh = qconj31(qz_est(2:5,k));     %;qz_est(2:5,k)
    qd = 2 * qmul31(qdq,qdh);  %无需进行归一化处理 前三位 为w 后一位接近于0
%     qdj = qnormlz(qd);
%     qdjz=[qdjz qdj];
    wd = qd(1:3);
    wdz(:,k)=[qz_est(1,k); wd];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%两种方法的 精度
wz_est  = wqkz(1:4,:); %时间 角速率
dwz_est = wqzI(1:4,:) - wz_est;
[wrmse_est, wSTD_est, wMeanz_est, wmaxz_est,wminz_est] = accuracy (dwz_est(2:4,n+2:end));
accw_est=[wrmse_est, wminz_est, wmaxz_est, wMeanz_est, wSTD_est];

dwdz = wqzI(1:4,:) - wdz;
[rmse_wd, STD_wd, Meanz_wd, maxz_wd,minz_wd] = accuracy (dwdz(2:4,n+2:end));
acc_wd=[rmse_wd,minz_wd, maxz_wd, Meanz_wd, STD_wd];

% w = wqzI(2:4,:);   
% figurew3(w(:,n+2:end),wdz(2:4,n+2:end),wz_est(2:4,n+2:end),m,n+2)
figurewerr(dwdz(2:4,n+2:end),dwz_est(2:4,n+2:end),m,n+2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%计算四元数误差
alfa_est = zeros(4,m);
alfa_mea = zeros(4,m);
for k=n+1:m
    %计算乘性四元数
    qzkconj = qconj31(qz(:,k));
    alfa4 = qmul31(qz_est(2:5,k), qzkconj); 
    alfa4_ = qmul31(qz_(:,k), qzkconj);  

    %归一化 保存 乘性四元数
    alfa41 = qnormlz(alfa4);
    alfa_est(:,k)= alfa41;
     
    alfa41_ = qnormlz(alfa4_);
    alfa_mea(:,k)= alfa41_;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rmse_est, STD_est, Meanz_est, maxz_est,minz_est] = accuracy (alfa_est(:,n+1:end));
accq_est=[rmse_est, minz_est, maxz_est, Meanz_est, STD_est];

[rmse_mea, STD_mea, Meanz_mea, maxz_mea,minz_mea] = accuracy (alfa_mea(:,n+1:end));
accq_mea=[rmse_mea, minz_mea, maxz_mea, Meanz_mea, STD_mea];

[rmseZ, STDZ, MeanzZ, maxzZ,minzZ] = accuracy (alfaZ(:,n+1:end));
accZ=[rmseZ,minzZ, maxzZ, MeanzZ, STDZ];

figureqerr(alfa_mea(:,n+1:end),alfa_est(:,n+1:end),m,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pxx_true = periodogram(wqzI(2,:));
% Pxx_mea = periodogram(alfa_mea(1,n+1:end));   
% Pxx_est = periodogram(alfa_est(1,n+1:end));
% zzz=[Pxx_mea Pxx_est];

% h=0.1;
% attz = [];
% for i=n+1:m
%     att1 = q2att31(qz_est(2:5,i));
%     att2 = q2att31(qz_(:,i));
%     att12 = [qz_est(1,i); att1; att2];
%     attz = [attz att12];  % 时间 估计的姿态角 原始的姿态
% end
% att_est  = attz(2:4,:)./pi*180; %arcsec  
% att_orig = attz(5:7,:)./pi*180; %arcsec ./pi*180
% [psd_Aest,psd_Aest_f,psd_Aorig,psd_Aorig_f]= psdcom (att_est,att_orig,h);
% figurepsd(psd_Aest,psd_Aorig,psd_Aest_f)

% Pxx_est = periodogram(wz_est(2,:));


% 
% save Pxx_true Pxx_true
% save Pxx_est Pxx_est