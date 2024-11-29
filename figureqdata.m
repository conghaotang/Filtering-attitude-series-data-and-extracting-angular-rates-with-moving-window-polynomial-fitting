function  figureqdata(w_est,w_mea,m,n,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%½ÇËÙÂÊµÄÎó²îÍ¼
  epochs=m:h:n;
%  epochs=1:length(w_mea);
  figure
  set(gcf,'color',[1,1,1])
  
  subplot(3,1,1);
  plot(epochs, w_mea(1,:),'-b','Linewidth',1.5);   %,'.','Linewidth',2.5
  hold on
  plot(epochs,w_est(1,:),'-r','Linewidth',1.5);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
%   set(gca,'XLim',[0 5400]);
     set(gca,'YLim',[-100 -40]);
     set(gca,'YTick',[-100 -50]);
  set(gca,'XLim',[0 3600]);
%   set(gca,'XTick',[0:2:20]);
  
  subplot(3,1,2);
  plot(epochs, w_mea(2,:),'-b','Linewidth',1.5);  ylabel('Relative Attitude [arcsec]'); 
  hold on
  plot(epochs,w_est(2,:),'-r','Linewidth',1.5);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
%   set(gca,'XLim',[0 5400]);
  set(gca,'XLim',[0 3600]);
  set(gca,'YLim',[-160 -100]);
     set(gca,'YTick',[-150 -100]);
     legend('Origial','Filtered')
%   set(gca,'XTick',[0:2:20]);
    
  subplot(3,1,3);
  plot(epochs, w_mea(3,:),'-b','Linewidth',1.5); xlabel('Time [1s]');  %ylabel('Roll'); %,'Linewidth',1 ,'Linewidth',1
  hold on
  plot(epochs,w_est(3,:),'-r','Linewidth',1.5);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
%   set(gca,'XLim',[0 5400]);
  set(gca,'XLim',[0 3600]);
  set(gca,'YLim',[-400 200]);
     set(gca,'YTick',[-400 -200 0 200]);
%   set(gca,'XTick',[0:2:20]);

  
 