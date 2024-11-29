clear all;
close all;
clc;
%qcase1中  clo = 998;  n=8;  theta1=0.01;   theta2=0.01;
%qcase2中  clo = 998;  n=18;  theta1=0.01;   theta2=0.02;
clo = 998; %sum(100*clock); %1000;
randn('state',clo); %0  1e+03   2.11e+05  2.1041705 2.10414e+05

m=200;  %序列总长度
n=18;  %窗口长度 
theta1=0.01;   theta2=0.02; %theta3=0.05;
[qz,qz_,etaZ,alfaZ]=dataq(m,theta1,theta2);
P=theta2^2*eye(3);  %Pk=P1; P1=theta2^2*eye(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
beta=zeros(9,m);   %未知数  
GG=zeros(3,9,m);   %未知数的系数
g_est=zeros(3,m);  %Gibbs 的估计 
q_est=zeros(4,m);  %真实四元数的估计值
for k=n+1:m   %k需要从 n+1开始 因为需要 历元1作为参考
    
    i=1; % yk 中 g_ 所在的位置
    
    %窗口长度内 观测数据   
    NNz=zeros(9);  % 法矩阵 BTPB
    Wz=zeros(9,1); %  BTPL
    G=zeros(3,9,n+1); %mm=[];
    
    q_kn = qz_(:,k-n);
    qz_KNconj = qconj31(q_kn);
    for j= k:-1: (k-n+1)
        
        p_j = qmul31(qz_(:,j), qz_KNconj); %是否需要进行 归一化  
        p_jj = qnormlz(p_j);
        g_j = q2g(p_jj);
        
        Fj=eye(3) + askew(g_j) + g_j * g_j';
        Qj=Fj * P * Fj';
        
        G1=eye(3);  G2=(j-k)*eye(3);  G3=(j-k)^2*eye(3); %令时间间隔为 1 s
        G(:,:,i)=[G1 G2 G3];
        
        Gj=G(:,:,i);        
        NNz = NNz + Gj' * Qj^-1 * Gj;
        Wz = Wz + Gj' * Qj^-1 * g_j;
        i=i+1;
    end

    GG(:,:,k)= G(:,:,1);
    beta(:,k)=NNz^-1 * Wz;
    g_est(:,k) = GG(:,:,k) * beta(:,k);  %估计的 Gibbs 矢量
    
    g_estk = g_est(:,k);
    gg = 1/sqrt(1 + g_estk' * g_estk);  % 这里的系数 暂时设定为 +1
    
    p_est_k = gg * [g_estk; 1]; %相对的 四元数
    q_est_k = qmul31(p_est_k,q_kn);
    q_est(:,k) = q_est_k;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alfa_est = zeros(4,m);
alfa_mea = zeros(4,m);
for k=n+1:m
    %计算乘性误差四元数
    qzkconj = qconj31(qz(:,k));
    alfa4 = qmul31(q_est(:,k), qzkconj); 
    alfa4_ = qmul31(qz_(:,k), qzkconj);  %应该与alfaz相同

    %归一化 保存 乘性四元数
    alfa41 = qnormlz(alfa4);
    alfa_est(:,k)= alfa41;
     
    alfa41_ = qnormlz(alfa4_);
    alfa_mea(:,k)= alfa41_;
end

[rmse_est, STD_est, Meanz_est, maxz_est,minz_est] = accuracy (alfa_est(:,n+1:end));
acc_est=[rmse_est, STD_est, Meanz_est, maxz_est,minz_est];

[rmse_mea, STD_mea, Meanz_mea, maxz_mea,minz_mea] = accuracy (alfa_mea(:,n+1:end));
acc_mea=[rmse_mea, STD_mea, Meanz_mea, maxz_mea,minz_mea];

[rmseZ, STDZ, MeanzZ, maxzZ,minzZ] = accuracy (alfaZ(:,n+1:end));
accZ=[rmseZ, STDZ, MeanzZ, maxzZ,minzZ];

figurea(alfa_mea(:,n+1:end),alfa_est(:,n+1:end),m,n)

