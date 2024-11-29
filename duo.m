
clc;
clear;
randn('state',1000); %0  1

betaZ=[]; mea=zeros(1,100); 
 truef=zeros(1,100);  truef(1,1)=1;
 theta1=0.001;   theta2=1; %theta3=0.05;
 for i=1:100
%      beta=theta1*randn(1,1);
%      truef(1,i)=truef(1,i-1)+ beta;
   truef(1,i)=sin(0.5*i+0.03);
 end
 
 for i=1:100
     alfa=theta2*randn(1,1);
     mea(1,i)= alfa+ truef(1,i);
 end
 
 n=3; P=theta2^2*eye(4) ; g_est=zeros(1,100); %即true
 for k=n+1:100
     GG=[]; meaz=[];
     for j=k:-1:k-n
         G1=1;  G2=(j-k);  G3=(j-k)^2; %令时间间隔为 1 s
         Gi=[G1 G2 G3];
         GG=[GG; Gi];
         meaz=[meaz;mea(:,j)];
     end
     NN= GG'* P^-1 * GG;
     Wz = GG'* P^-1 * meaz;
     beta=NN^-1 * Wz;
     VV=GG * beta - meaz;
     g_est(:,k) = beta(1);
 end
 meaerr = mea(4:100)-truef(4:100);
 esterr = g_est(4:100)-truef(4:100);
 
 rmsee=rms(meaerr);
 rmse1=rms(esterr);
 figure
%  plot(4:100,meaerr,'b')
  plot(1:100,truef,'b')
 hold on
%  plot(4:100,esterr,'r')
 plot(1:100,g_est,'r')
  hold on
  plot(1:100,mea,'g')
  
%  legend('meaerr','esterr')

%  [rmse, STD,  Meanz, maxz,minz] = accuracy (meaerr);
%  acc=[rmse, STD,  Meanz, maxz,minz];
%  [rmse1, STD1,  Meanz1, maxz1,minz1] = accuracy (esterr);
%  acc1=[rmse1, STD1,  Meanz1, maxz1,minz1];
 