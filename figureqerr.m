function  figureqerr(q_mea,q_est,m,n)
%四元数及角速率的误差图

  epochs=(n+1)*0.1:0.1:m*0.1;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  figure
  subplot(3,1,1);
  plot(epochs, q_mea(1,:),'-b','Linewidth',3.5);   %ylabel('Yaw [arcmin]'); %xlabel('epoch')  %,'Linewidth',1
  hold on
  plot(epochs,q_est(1,:),'-r','Linewidth',3.5);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
  set(gca,'XLim',[0 21]);
  set(gca,'YLim',[-0.06 0.06]);
  set(gca,'YTick',[-0.05 0 0.05]);
  legend('Original','Filtered') %intial
   
  subplot(3,1,2);
  plot(epochs, q_mea(2,:),'b','Linewidth',3.5); ylabel('Attitude Error [rad]'); 
  hold on
  plot(epochs,q_est(2,:),'r','Linewidth',3.5);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
  set(gca,'XLim',[0 21]);
  set(gca,'YLim',[-0.06 0.06]);
  set(gca,'YTick',[-0.05 0 0.05]);
  
  subplot(3,1,3);
  plot(epochs, q_mea(3,:),'b','Linewidth',3.5);  xlabel('Time [0.1s]')  %,'Linewidth',1 ,'Linewidth',1
  hold on
  plot(epochs,q_est(3,:),'r','Linewidth',3.5);
  set(gca,'XLim',[0 21]);
  set(gca,'YLim',[-0.06 0.06]);
  set(gca,'YTick',[-0.05 0 0.05]);
  
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
