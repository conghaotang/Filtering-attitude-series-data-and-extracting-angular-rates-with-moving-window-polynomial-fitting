clear all;
close all;
clc;
% ss=1e+03; %sum(100*clock);
randn('state',1e+03); %0  1e+03   2.11e+05  2.1041705 2.10414e+05

m=100;  %序列总长度
theta1=0.01;   theta2=0.01; theta3=0.05;
[RZ,R_,alfaZ]=data(m,theta1,theta2,theta3);  %真值 估计值 估计值的乘性误差

P=zeros(3,3,m);
for i=1:m
    P(:,:,i)=theta2^2*eye(3);  %Pk=P1; P1=theta2^2*eye(3);
end
for i=[51 55 59]
    P(:,:,i)=theta3^2*eye(3);
end

n=3;  %窗口长度 为 n+1
S_=zeros(3,3,n+1); g=zeros(3,n+1);
beta=zeros(9,m);  GG=zeros(3,9,m);
g_est=zeros(3,m);  S_est=zeros(3,3,m);  R_est=zeros(3,3,m); 
g_mea=zeros(3,m);  mmz=[];

for k=n+1:m
    
    i=1; %论文中的j
 
    %窗口长度内 观测数据   
    NNz=zeros(9);  Wz=zeros(9,1); G=zeros(3,9,n+1); mm=[];
    for j=k:-1:k-n  
        S_(:,:,j) = R_(:,:,j)*R_(:,:,k)';
        S_j=S_(:,:,j);
%         q_ = m2qua(S_j);
%         g_(:,i)=qua2g(q_);  %观测量顺序为 gk  gk-n
        
        [q_x,maxm]=m2q(S_j);
        g_(:,i)=q2g(q_x);
%         mm=[mm maxm];
        g_i=g_(:,i);
        Fi=-1/2*(1+norm(g_i)^2)*(eye(3)-askew(g_i))^-1;
        Q=Fi*P(:,:,j)*Fi';
        
        G1=eye(3);  G2=(j-k)*eye(3);  G3=(j-k)^2*eye(3); %令时间间隔为 1 s
        G(:,:,i)=[G1 G2 G3];
        
        Gi=G(:,:,i);        
        NNz = NNz + Gi' * Q^-1 * Gi;
        Wz = Wz + Gi' * Q^-1 * g_i;
        i=i+1;
    end
%     mmz=[mmz mm];
    g_mea(:,k)=g_(:,1);
    GG(:,:,k)= G(:,:,1);
    beta(:,k)=NNz^-1 * Wz;
    g_est(:,k) = GG(:,:,k) * beta(:,k);
end

alfa_y=zeros(3,m);   alfa_y2=zeros(3,m);
for k=n+1:m
    eyea= eye(3) + askew(g_est(:,k));
    eyed= eye(3) - askew(g_est(:,k));
    S_est(:,:,k) = eyea^-1 * eyed;
    R_est(:,:,k) = S_est(:,:,k) * R_(:,:,k);
    
    Lk =  R_est(:,:,k) * RZ(:,:,k)'; 
    Lkk=  RZ(:,:,k) * R_est(:,:,k)'; 
    alfa_y(:,k)=iaskew(1/2*(Lk-Lkk));
    
%     Lk2 =  R_(:,:,k) * RZ(:,:,k)';
%     Lkk2=  RZ(:,:,k) * R_(:,:,k)';
%     alfa_y2(:,k)=iaskew(1/2*(Lk2-Lkk2));

end

[rmse, STD, Meanz, maxz,minz] = accuracy (alfaZ(:,5:end));
 acc=[rmse, STD, Meanz, maxz,minz];
 
 [rmse_y, STD_y, Meanz_y, maxz_y,minz_y] = accuracy (alfa_y(:,5:end));
 acc_y=[rmse_y, STD_y, Meanz_y, maxz_y,minz_y];

%  figurea(alfaZ,alfa_y,m)

%  qnb=[1 0 0 0];
%  Cnb = q2mat(qnb)
 
%  R = RZ(:,:,10)
%  q_ = m2qua(R);
%  Cnb = q2mat(q_)
%  
%  arc=R-Cnb 
% q_ = m2qua(arc)
% qv=q_(2:4,1);
%   g=qv/q_(1)
%   
%   g=[0.1 0.3 0.2]';
%   daoxu = (eye(3)+askew(g))^-1 * (eye(3)-askew(g))
%   true =  (eye(3)-askew(g)) * (eye(3)+askew(g))^-1

