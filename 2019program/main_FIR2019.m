clear all;
close all;
clc;
% ss=sum(100*clock); %1001
randn('state',2.11e+05); %0  1e+03   2.11e+05  2.1041705 2.10414e+05

% for i =1:10
%       rr(:,i)=normrnd(0,1,3,1);
% end

m=100;  %序列总长度
n=6;  %窗口长度 为 n+1
theta1=0.001;   theta2=0.1; %theta3=0.05;
[RZ,R_]=data(m,theta1,theta2);  %真值 估计值

P=zeros(3,3,m);
for i=1:m
    P(:,:,i)=theta2^2*eye(3);  %Pk=P1; P1=theta2^2*eye(3);
end

S_=zeros(3,3,n+1); g=zeros(3,n+1); %S_est=zeros(3,3,m);
beta=zeros(9,m);  GG=zeros(3,9,m);
g_est=zeros(3,m);    R_est=zeros(3,3,m); 
g_mea=zeros(3,m);  mmz=[];

for k=n+1:m

    %窗口长度内 观测数据   
    GG=[]; QZ=[]; gz=[];
    for j= k:-1: (k-n)
        S_j = R_(:,:,j)*R_(:,:,k)';
        
        [q_x,maxm]=m2q(S_j); %第四项 为 标量
        g_i=q2g(q_x);

        Fi=-1/2 * (1+g_i(1)^2+g_i(2)^2+g_i(3)^2) * (eye(3)-askew(g_i))^-1;
        Q=Fi * P(:,:,j) * Fi';
        
        G1=eye(3);  G2=(j-k)*eye(3);  G3=(j-k)^2*eye(3); %令时间间隔为 1 s
        Gi=[G1 G2 G3];
        GG=[GG; Gi];
        QZ=[QZ;  Q];
        gz=[gz; g_i];
%         NNz = NNz + Gi' * Q^-1 * Gi;
%         Wz = Wz + Gi' * Q^-1 * g_i;
    end
    QZdiag = diagq(QZ);
    NNz= GG'* QZdiag^-1 * GG;
    Wz = GG'* QZdiag^-1 * gz;
    beta=NNz^-1 * Wz;  
    VV=GG * beta - gz;
    g_est(:,k) = beta(1:3);
end

alfa_y=zeros(3,m);   alfa_y2=zeros(3,m);
% att_est=zeros(3,m); 
for k=n+1:m
    eyea= eye(3) + askew(g_est(:,k));
    eyed= eye(3) - askew(g_est(:,k));
    S_est = eyea^-1 * eyed;
    R_est(:,:,k) = S_est * (R_(:,:,k)')^-1;
    
%     att_est(:,k) = m2att(R_est(:,:,k));
%     att_Z(:,k) = m2att(RZ(:,:,k));
%     att_(:,k)= m2att(R_(:,:,k));
    Lk =  R_est(:,:,k) * RZ(:,:,k)'; 
    Lkk=  RZ(:,:,k) * R_est(:,:,k)'; 
    alfa_y(:,k)=iaskew(1/2*( Lk - Lkk ));
    
    Lk2 =  R_(:,:,k) * RZ(:,:,k)';
    Lkk2=  RZ(:,:,k) * R_(:,:,k)';
    alfa_y2(:,k)=iaskew(1/2*( Lk2 - Lkk2 ));

end

% att_EZ = att_est - att_Z;
% att_YZ = att_ - att_Z;
% 
% rrr=R_est-R_;

[rmse, STD, Meanz, maxz,minz,rmsezintial] = accuracy (alfa_y2(:,n+1:end));
 acc=[rmse, STD, Meanz, maxz,minz];  %原始
 
 [rmse_y, STD_y, Meanz_y, maxz_y,minz_y,rmsezyou] = accuracy (alfa_y(:,n+1:end));
 acc_y=[rmse_y, STD_y, Meanz_y, maxz_y,minz_y]; %重构

%  figurea(alfa_y2(:,n+1:end),alfa_y(:,n+1:end),m,n)
 
       
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

