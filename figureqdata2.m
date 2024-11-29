function  figureqdata2(w_est,w_mea,m,n,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�����ʵ����ͼ
  epochs=m:h:n;
%  epochs=1:length(w_mea);
  figure
  set(gcf,'color',[1,1,1])
  
  subplot(3,1,1);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
  plot(epochs, w_mea(1,:),'-b','Linewidth',1);   %,'.','Linewidth',2.5
  hold on
  plot(epochs,w_est(1,:),'-r','Linewidth',0.5);  
  set(gca,'XLim',[0 5400]);
    set(gca,'YLim',[-4 4]);
%   set(gca,'XLim',[0 3600]);
%   set(gca,'XTick',[0:2:20]);
  
  subplot(3,1,2);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
  plot(epochs, w_mea(2,:),'-b','Linewidth',1);  ylabel('Relative Attitude [arcsec]'); 
  hold on
  plot(epochs,w_est(2,:),'-r','Linewidth',0.5); 
  set(gca,'XLim',[0 5400]);
%   set(gca,'XLim',[0 3600]);
%   set(gca,'XTick',[0:2:20]);
    
  subplot(3,1,3);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
  plot(epochs, w_mea(3,:),'-b','Linewidth',1); xlabel('Time [1s]');  %ylabel('Roll'); %,'Linewidth',1 ,'Linewidth',1
  hold on
  plot(epochs,w_est(3,:),'-r','Linewidth',0.5);  
  set(gca,'XLim',[0 5400]);
%   set(gca,'XLim',[0 3600]);
%   set(gca,'XTick',[0:2:20]);

  legend('Origial','Filtered')
 