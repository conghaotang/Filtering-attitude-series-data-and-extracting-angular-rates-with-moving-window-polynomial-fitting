function  figurewerr(w_mea,w_est,m,n)

%½ÇËÙÂÊµÄÎó²îÍ¼
  epochs=n*0.1:0.1:m*0.1;
  figure
  set(gcf,'color',[1,1,1])
  
  subplot(3,1,1);
  plot(epochs, w_mea(1,:),'Linewidth',3.5);   
  hold on
  plot(epochs,w_est(1,:),'Linewidth',3.5);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
  set(gca,'XLim',[0 21]);
  set(gca,'XTick',[0:2:20]);
  set(gca,'YLim',[-0.5 0.5]); %qcase2
  set(gca,'YTick',[-0.5 0 0.5]); %qcase2
  
  subplot(3,1,2);
  plot(epochs, w_mea(2,:),'Linewidth',3.5);  ylabel('Angular Rate Estimation Error [rad/s]'); 
  hold on
  plot(epochs,w_est(2,:),'Linewidth',3.5);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
  set(gca,'XLim',[0 21]);
  set(gca,'XTick',[0:2:20]);
  set(gca,'YLim',[-0.5 0.5]); %qcase2
  set(gca,'YTick',[-0.5 0 0.5]); %qcase2
    
  subplot(3,1,3);
  plot(epochs, w_mea(3,:),'Linewidth',3.5); xlabel('Time [0.1s]');  %ylabel('Roll'); %,'Linewidth',1 ,'Linewidth',1
  hold on
  plot(epochs,w_est(3,:),'Linewidth',3.5);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
  set(gca,'XLim',[0 21]);
  set(gca,'XTick',[0:2:20]);
  set(gca,'YLim',[-0.5 0.5]); %qcase2
  set(gca,'YTick',[-0.5 0 0.5]); %qcase2
%   ylabel('Roll');
  legend('Difference','Proposed')

%%%%%%%%%%%%%%%%%
 