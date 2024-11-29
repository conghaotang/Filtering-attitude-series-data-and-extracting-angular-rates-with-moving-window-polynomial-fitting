clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num = xlsread('E:\FIR\DATA20210110C1b.xls');  %ע�� ȥ�� xls�еļ�����datam.xls
times=num(:,2)';
qz_1 = num(:,3:6)'; %��������Ԫ��
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
beta=zeros(9,m);   %δ֪��  
GG=zeros(3,9,m);   %δ֪����ϵ��
g_est=zeros(3,m);  %Gibbs �Ĺ��� 
qt_est=zeros(5,m);  %��ʵ��Ԫ���Ĺ���ֵ
for k=n+1:m   %k��Ҫ�� n+1��ʼ ��Ϊ��Ҫ ��Ԫ1��Ϊ�ο�
    
    i=1; % yk �� g_ ���ڵ�λ��
    
    %���ڳ����� �۲�����   
    NNz=zeros(9);  % ������ BTPB
    Wz=zeros(9,1); %  BTPL
    G=zeros(3,9,n+1); %mm=[];
    
    q_kn = qz_(:,k-n);
    qz_KNconj = qconj31(q_kn);
    
    for j= k:-1: (k-n+1)
        p_j = qmul31(qz_(:,j), qz_KNconj);%����Ҫ���й�һ�� ��λ��Ԫ���˵�λ��Ԫ�� ��Ϊ��λ��Ԫ��  
        p_jj = qnormlz(p_j);
        g_j = q2g(p_jj);
        
        Fj=eye(3) + askew(g_j) + g_j * g_j';
        Qj=Fj * P * Fj';
        
        G1=eye(3);  G2=(j-k) * td * eye(3);  G3=((j-k) * td)^2 * eye(3); %��ʱ����Ϊ td s
        Gj=[G1 G2 G3];
                
        NNz = NNz + Gj' * Qj^-1 * Gj;
        Wz = Wz + Gj' * Qj^-1 * g_j;
        i=i+1;
    end

    beta =NNz^-1 * Wz;
    g_est = beta(1:3);  %���Ƶ� Gibbs ʸ��
    
    p_est_k = g2q(g_est);
%     gg = 1/sqrt(1 + g_estk' * g_estk);  % �����ϵ�� ��ʱ�趨Ϊ +1
%     p_est_k = gg * [g_estk; 1]; %��Ե� ��Ԫ��
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

figureqdata(attz(2:4,:),attz(5:7,:),length(attz),1) %�˲��� ԭʼ��
attd = attz(2:4,:)-attz(5:7,:);

Pxx_mea = periodogram(att_mea(1,n+1:end));   
Pxx_est = periodogram(att_est(1,n+1:end));
zzz=[Pxx_mea Pxx_est];


% [Pxx,f]=psd(attz(2,:),Nfft,fs,window,noverlap,dflag);
% [Pxx,f]=psd(x,Nfft,fs,window,noverlap,dflag);
%  matlab�е�psd�÷�:PSD ����������bai�ܶȵĺ���du��
% x���źţ�zhiNfft���ٸ���Ҷ�任������fs�ǲ���Ƶdao�ʣ�
%window�Ǽӵ�zhuan��������noverlap��ָû���ص���shu�����ݲ���������������С����Ƶ�ʣ���dflag�������ж�ǰ�����noverlap�Ƿ����ص����еĻ��ͼ٣��������������û���ص�Ϊ�棨����������

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
  
  
%������ÿһ����Ԫ ��ʵ����̬���� ���޷�������������Ԫ��
%�÷��������� ��Ƶ����


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%��Ԫ�� �� ��̬��
%  figurea(qt_est(2:5,n+1:end),qz_(:,n+1:end),m,n)


