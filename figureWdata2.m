function  figureWdata2(w_est,w_mea,m,n,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%½ÇËÙÂÊµÄÎó²îÍ¼
  epochs=m:h:n;
%  epochs=1:length(w_mea);
  figure
  set(gcf,'color',[1,1,1])
  
  subplot(3,1,1);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
  plot(epochs, w_mea(1,:),'Linewidth',1);   %,'.','Linewidth',2.5
  hold on
  plot(epochs,w_est(1,:),'Linewidth',0.5);
  set(gca,'XLim',[0 5400]);
  set(gca,'YLim',[-4 4]);
%   set(gca,'XLim',[0 3600]);
%   set(gca,'XTick',[0:500:3500]); %'Angular Rate Estimation Error [rad/s]'
  
  subplot(3,1,2);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
  plot(epochs, w_mea(2,:),'Linewidth',1);  ylabel('Angular Rate Estimation [¡ã/h]'); 
  hold on
  plot(epochs,w_est(2,:),'Linewidth',0.5);
  set(gca,'XLim',[0 5400]);
%   set(gca,'XLim',[0 3600]);
%   set(gca,'XTick',[0:500:3500]);
    
  subplot(3,1,3);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
  plot(epochs, w_mea(3,:),'Linewidth',1); xlabel('Time [1s]');  %ylabel('Roll'); %,'Linewidth',1 ,'Linewidth',1
  hold on
  plot(epochs,w_est(3,:),'Linewidth',0.5);
  set(gca,'XLim',[0 5400]);
%   set(gca,'XLim',[0 3600]);
%   set(gca,'XTick',[0:500:3500]);

  legend('Difference','Proposed')
 