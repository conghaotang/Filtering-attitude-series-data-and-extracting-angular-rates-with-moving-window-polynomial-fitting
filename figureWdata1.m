function  figureWdata1(w_est,w_mea,m,n,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%½ÇËÙÂÊµÄÎó²îÍ¼
  epochs=m:h:n;
%  epochs=1:length(w_mea);
  figure
  set(gcf,'color',[1,1,1])
  
  subplot(3,1,1);
  plot(epochs, w_mea(1,:),'Linewidth',1.5);   %,'.','Linewidth',2.5
  hold on
  plot(epochs,w_est(1,:),'Linewidth',1.5);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
  set(gca,'YLim',[0.016 0.022]);
  set(gca,'YTick',[0.016 0.018 0.02]);
  set(gca,'XLim',[0 3600]);
  set(gca,'XTick',[0:500:3500]); %'Angular Rate Estimation Error [rad/s]'
  
  subplot(3,1,2);
  plot(epochs, w_mea(2,:),'Linewidth',1.5);  ylabel('Angular Rate Estimation [¡ã/s]'); 
  hold on
  plot(epochs,w_est(2,:),'Linewidth',1.5);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
%   set(gca,'XLim',[0 5400]);
set(gca,'YLim',[0.034 0.04]);
  set(gca,'YTick',[0.035 0.04]);
  set(gca,'XLim',[0 3600]);
  set(gca,'XTick',[0:500:3500]);
  legend('Difference','Proposed')
    
  subplot(3,1,3);
  plot(epochs, w_mea(3,:),'Linewidth',1.5); xlabel('Time [1s]');  %ylabel('Roll'); %,'Linewidth',1 ,'Linewidth',1
  hold on
  plot(epochs,w_est(3,:),'Linewidth',1.5);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
%   set(gca,'XLim',[0 5400]);
set(gca,'YLim',[0.02 0.08]);
  set(gca,'YTick',[0.02 0.04 0.06 0.08]);
  set(gca,'XLim',[0 3600]);
  set(gca,'XTick',[0:500:3500]);


  
 