clear all;
close all;
clc;
% clo = 996 ;  % 1994 2000 1e+03;  1e+03; %1000;sum(100*clock) sum(10*clock)
% randn('state',clo); %0  1e+03     2.1041705 2.10414e+05
%%%%%qwcase1 clo = 1000;  n=18;  theta1=0.001;   theta2=0.05;
%%%%%qwcase1 clo = 996;  n=14;  theta1=0.0001;   theta2=0.01;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��ȡ���ݲ����ò�����10s 5s
%{
num = xlsread('E:\FIR\D3a.xls');  %ע�� ȥ�� xls�еļ�����DATA20210110C1b.xls
%��ȡ���� �ٽ�������ɸѡ
%�ѽ������е���Ԫ�� ʸ�� ������ʽ
qSCA123 = num;
% save qSCA123 qSCA123
% load('qSCA123.mat');
a = find(qSCA123(:,3)==3);  %ѡȡSCA3
qSCA1 = qSCA123(a,:);

j=1; D3a=[];
for i=1:length(qSCA1)   
  if qSCA1(j,1)==qSCA1(j+1,1)
      D3a=[D3a; qSCA1(j,:)];
      j=j+2;
  else
      D3a=[D3a; qSCA1(j+1,:)];
      j=j+1;
  end
end
% save zqSCA1 zqSCA1 %663511946
% aa = find(zqSCA1(:,1) == 663511946);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save D3a D3a
load('D3a.mat');
num = D3a(:,[1:2,5:7,4]); %����Ԫ�� ʸ�� ����
num(:,2)=num(:,1)-663500000;

qz_1 = num(:,2:6)'; %ʱ�� ��������Ԫ��
mz=length(qz_1);

% j=1; 
% for i=1:h:200
%     qz_t(:,j) = qz_1(:,i);
%     j=j+1;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qz_t = qz_1;
times= qz_t(1,:);
td = times(2) - times(1); %ʱ����
qz_ = qz_t(2:5,:);        %ԭʼ��Ԫ��
m = 3600;
h=td;
n=15; %���ڳ��� 
P = [1 0 0;0 1 0;0 0 1]*eye(3);  %(10^6)* (10^-3)*Pk=P1; P1=theta2^2*eye(3); [1000 0 0;0 1000 0; 0 0 1000]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beta=zeros(9,m);   %δ֪��  
% g_est=zeros(3,m);  %Gibbs �Ĺ��� 
% q_est=zeros(4,m);  %��ʵ��Ԫ���Ĺ���ֵ
% wqkz = zeros(8,m); 

for k=n+1:m   %k��Ҫ�� n+1��ʼ ��Ϊ��Ҫ ��Ԫ1��Ϊ�ο�
    
    i=1; % yk �� g_ ���ڵ�λ��
    %���ڳ����� �۲�����   
    NNz=zeros(9);  % ������ BTPB
    Wz=zeros(9,1); %  BTPL
    
    q_kn = qz_(:,k-n);
    h_kn = q2g(q_kn);
    qz_KNconj = qconj31(q_kn);
   
    for j= k:-1: (k-n+1) 
        p_j  = qmul31(qz_(:,j), qz_KNconj);   
        p_jj = qnormlz(p_j);  % ��һ��
        g_j = q2g(p_jj);
        
        Fj=eye(3) + askew(g_j) + g_j * g_j';
        Qj=Fj * P * Fj';
        
        G1=eye(3);  G2=(j-k)*td*eye(3);  G3=((j-k)*td)^2 * eye(3); %��ʱ����Ϊ 1 s
        Gj=[G1 G2 G3];
       
        NNz = NNz + Gj' * Qj^-1 * Gj;
        Wz = Wz + Gj' * Qj^-1 * g_j;
        i=i+1;
    end

    beta = NNz^-1 * Wz;
    g_est = beta(1:3);  %���Ƶ� Gibbs ʸ��  (:,k):
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  %��ȡ ������
    hk_est = gmul(g_est,h_kn);
    g_est_f = beta(4:6);   % Gibbs��΢��
%     hk_est_f = gmul(g_est_f,h_kn); % h���Ƶ�΢��
    hk_est_f = gmulf(g_est,h_kn,g_est_f);% h���Ƶ�΢��  gmul(g_est_f,h_kn);
  
    wkfz = 2/(1 + hk_est' * hk_est);
    wkfm = ( hk_est_f - askew(hk_est) * hk_est_f);
    wk_est = wkfz * wkfm ;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %���� ��Ԫ��
    p_est_k = g2q(g_est);
    q_est_k = qmul31(p_est_k,q_kn);
    
    wqkz(:,k)= [qz_t(1,k); wk_est; q_est_k];  %����ֵ ʱ�� ������ ��Ԫ�� 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�Ա� ���� ���� ������
wz_est  = wqkz(1:4,:); %ʱ�� ������
qz_est = wqkz([1,5:8],:); %ʱ�� ��Ԫ��
for k=n+2:m 
    qdq = (qz_est(2:5,k) - qz_est(2:5,k-1))/(qz_est(1,k)-qz_est(1,k-1));
    qdh = qconj31(qz_est(2:5,k));     
    qd = 2 * qmul31(qdq,qdh);  %������й�һ������ ǰ��λ Ϊw ��һλ�ӽ���0

    wd = qd(1:3);
    wdz(:,k)=[qz_est(1,k); wd];  %��ַ�������Ľ��ٶ���Ϣ  ʱ�� ���ٶ�
end

wz_estq = wz_est(2:4,n+2:m)./pi*180;   %  ��/s        qz_est(1,n+2)
wdzq = wdz(2:4,n+2:m)./pi*180;         %  ��/s 
figureWdata1(wz_estq,wdzq,n+2,length(qz_est),h) %�˲��� ԭʼ��
[psd_West,psd_West_f,psd_Worig,psd_Worig_f]= psdcom (wz_estq,wdzq,h);
figureWpsd(psd_West,psd_Worig,psd_West_f) 

[Wrmse_est, WSTD_est, WMeanz_est, Wmaxz_est,Wminz_est] = accuracy (wz_estq);
[Wrmse_d, WSTD_d, WMeanz_d, Wmaxz_d,Wminz_d] = accuracy (wdzq);
acc_W_est = [Wrmse_est, WSTD_est, WMeanz_est, Wmaxz_est,Wminz_est]; 
acc_W_d = [Wrmse_d, WSTD_d, WMeanz_d, Wmaxz_d,Wminz_d];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�����̬�仯
qqr_ = []; qqr_est=[];
for i=n+1:m-1   % n+1 Ϊ ��̬ �� ��һ����
    qz_KNconj = qconj31(qz_est(2:5,i));
    qq =  qmul31(qz_est(2:5,i+1), qz_KNconj);
    qqz = [qz_est(1,i+1);qq];
    qqr_est = [qqr_est qqz];
    
    qz_iconj = qconj31(qz_(:,i));
    qq_ =  qmul31(qz_(:,i+1), qz_iconj);
    qqz_=[qz_t(1,i+1); qq_];
    qqr_ = [qqr_ qqz_];
end

attzd=[];
for i=1:length( qqr_ )
    att1 = q2att31(qqr_est(2:5,i));
    att2 = q2att31(qqr_(2:5,i));
    att12 = [qqr_(1,i); att1; att2];
    attzd = [attzd att12];  % ʱ�� ���Ƶ���̬�� ԭʼ����̬
end
att_estd  = attzd(2:4,:)./pi*180*3600; %�����̬�仯 arcsec 
att_origd = attzd(5:7,:)./pi*180*3600;
[Armse_estd, ASTD_estd, AMeanz_estd, Amaxz_estd,Aminz_estd] = accuracy (att_estd);
[Armse_origd, ASTD_origd, AMeanz_origd, wmaxz_origd,Aminz_origd] = accuracy (att_origd);

acc_att_estd = [Armse_estd, ASTD_estd, AMeanz_estd, Amaxz_estd,Aminz_estd]; 
acc_att_origd = [Armse_origd, ASTD_origd, AMeanz_origd, wmaxz_origd,Aminz_origd];

figureqdata(att_estd,att_origd,n+2,length(qz_est),h) %�˲��� ԭʼ��
%�� semilogy ��ʽ���ƹ�����
[psd_Aest,psd_Aest_f,psd_Aorig,psd_Aorig_f]= psdcom (att_estd,att_origd,h);
figureApsd(psd_Aest,psd_Aorig,psd_Aest_f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�������ܶ�ֱ�Ӽ��㷽�� �Ϲ�����ʦ �� + ����  �����������������ݷ���-�Ϲ���
% h = 1;
% N=1024;
% Xk = fft(att_est(1,:),N);
% fs=1/h;
% S1 = abs(Xk).^2/N/fs;  S1(2:end) = S1(2:end)*2; % ���㵥�߹�����
% N2 = floor(N/2);
% semilogy([0:N2-1]*fs/N, sqrt(S1(1:N2))); grid, xlabel('f / Hz'); ylabel('PSD');
% plot([0:N2-1]*fs/N, log10(sqrt(S1(1:N2))));

% w_est = wqkz(2:4,:);
% w_d   = wdz(2:4,:);

%��̬�����ٶȱ�ʾ ���� ������̬�����ٶȵĹ������ܶ� PSD
%ŷ����

attz = [];
for i=n+1:m
    att1 = q2att31(qz_est(2:5,i));
    att2 = q2att31(qz_(:,i));
    att12 = [qz_est(1,i); att1; att2];
    attz = [attz att12];  % ʱ�� ���Ƶ���̬�� ԭʼ����̬
end
att_est  = attz(2:4,:)./pi*180; %arcsec  
att_orig = attz(5:7,:)./pi*180; %arcsec ./pi*180
% attzheng = round(att_orig);
% att_est_zheng = att_est - attzheng;
% att_orig_zheng = att_orig - attzheng;
h=td;
figureqdata(att_est,att_orig,qz_est(1,n+1),qz_est(1,end),h) %�˲��� ԭʼ��

epochs=1:length(theta_att);
figure
plot(epochs,qz_(1,n+1:end),'.','Linewidth',2.5);
hold on
plot(epochs,qz_est(2,n+1:end),'.','Linewidth',1);

figure
plot(epochs, qz_(2,n+1:end),'.','Linewidth',2.5);
hold on
plot(epochs, qz_est(3,n+1:end),'.','Linewidth',1);

figure
plot(epochs, qz_(3,n+1:end),'.','Linewidth',2.5);
hold on
plot(epochs, qz_est(4,n+1:end),'.','Linewidth',1);

figure
plot(epochs, qz_(4,n+1:end),'.','Linewidth',2.5);
hold on
plot(epochs, qz_est(5,n+1:end),'.','Linewidth',1);


% figureqdata(qz_est(3:5,n+1:end),qz_(2:4,n+1:end),qz_est(1,n+1),qz_est(1,end),h) %�˲��� ԭʼ��
    
theta_att = att_est-att_orig;
epochs=1:length(theta_att);
figure
plot(epochs, theta_att(1,:),'.','Linewidth',2.5);
figure
plot(epochs, theta_att(2,:),'.','Linewidth',2.5);
figure
plot(epochs, theta_att(3,:),'.','Linewidth',2.5);

%}
% figureqdata(att_est_zheng,att_orig_zheng,qz_est(1,n+1),qz_est(1,end),h) %�˲��� ԭʼ��
% att_orig = dwdz(2:4,n+2:end);% alfa_mea(:,n+1:end)./pi*180*3600;
% att_est = dwz_est(2:4,n+2:end); %alfa_est(:,n+1:end)./pi*180*3600;

% figureqerr(alfa_mea(:,n+1:end),alfa_est(:,n+1:end),m,n)
% figurewerr(dwdz(2:4,n+2:end),dwz_est(2:4,n+2:end),m,n+2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%���ַ����� ����
% wz_est  = wqkz(1:4,:); %ʱ�� ������
% dwz_est = wqzI(1:4,:) - wz_est;
% [wrmse_est, wSTD_est, wMeanz_est, wmaxz_est,wminz_est] = accuracy (dwz_est(2:4,n+2:end));
% accw_est=[wrmse_est, wminz_est, wmaxz_est, wMeanz_est, wSTD_est];
% 
% dwdz = wqzI(1:4,:) - wdz;
% [rmse_wd, STD_wd, Meanz_wd, maxz_wd,minz_wd] = accuracy (dwdz(2:4,n+2:end));
% acc_wd=[rmse_wd,minz_wd, maxz_wd, Meanz_wd, STD_wd];

% w = wqzI(2:4,:);   
% figurew3(w(:,n+2:end),wdz(2:4,n+2:end),wz_est(2:4,n+2:end),m,n+2)
% figurewerr(dwdz(2:4,n+2:end),dwz_est(2:4,n+2:end),m,n+2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%������Ԫ�����
% alfa_est = zeros(4,m);
% alfa_mea = zeros(4,m);
% for k=n+1:m
%     %���������Ԫ��
%     qzkconj = qconj31(qz(:,k));
%     alfa4 = qmul31(qz_est(2:5,k), qzkconj); 
%     alfa4_ = qmul31(qz_(:,k), qzkconj);  
% 
%     %��һ�� ���� ������Ԫ��
%     alfa41 = qnormlz(alfa4);
%     alfa_est(:,k)= alfa41;
%      
%     alfa41_ = qnormlz(alfa4_);
%     alfa_mea(:,k)= alfa41_;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [rmse_est, STD_est, Meanz_est, maxz_est,minz_est] = accuracy (alfa_est(:,n+1:end));
% accq_est=[rmse_est, minz_est, maxz_est, Meanz_est, STD_est];
% 
% [rmse_mea, STD_mea, Meanz_mea, maxz_mea,minz_mea] = accuracy (alfa_mea(:,n+1:end));
% accq_mea=[rmse_mea, minz_mea, maxz_mea, Meanz_mea, STD_mea];
% 
% [rmseZ, STDZ, MeanzZ, maxzZ,minzZ] = accuracy (alfaZ(:,n+1:end));
% accZ=[rmseZ,minzZ, maxzZ, MeanzZ, STDZ];

% figureqerr(alfa_mea(:,n+1:end),alfa_est(:,n+1:end),m,n)



% 
% save Pxx_true Pxx_true
% save Pxx_est Pxx_est