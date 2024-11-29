function  figurea(alfaZ,alfa_y1,m,n)
  %alfaZ,alfa_y1  ‘≠ º  π¿À„
  epochs=n+1:m;
  figure 

  set(gcf,'color',[1,1,1])
  subplot(3,1,1);
  plot(epochs, alfaZ(1,:),'-b','Linewidth',3.5);   %ylabel('Yaw [arcmin]'); %xlabel('epoch')  %,'Linewidth',1
  hold on
  plot(epochs,alfa_y1(1,:),'-r','Linewidth',3.5);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
  set(gca,'XLim',[0 210]);
  set(gca,'YLim',[-0.1 0.1]); %qcase2
  set(gca,'YTick',[-0.1 0 0.1]); %qcase2
  %set(gca,'YTick',[-0.1 0 0.1]); %qcase1
   
  subplot(3,1,2);
  plot(epochs, alfaZ(2,:),'b','Linewidth',3.5); ylabel('Attitude Error [rad]'); 
  hold on
  plot(epochs,alfa_y1(2,:),'r','Linewidth',3.5);
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
  set(gca,'XLim',[0 210]);
  set(gca,'YLim',[-0.1 0.1]); %qcase2
  set(gca,'YTick',[-0.1 0 0.1]); %qcase2
  %set(gca,'YTick',[-0.1 0 0.1]); %qcase1
  
  subplot(3,1,3);
  plot(epochs, alfaZ(3,:),'b','Linewidth',3.5);  xlabel('Time [1s]')  %,'Linewidth',1 ,'Linewidth',1
  hold on
  plot(epochs,alfa_y1(3,:),'r','Linewidth',3.5);
  set(gca,'XLim',[0 210]);
  set(gca,'YLim',[-0.1 0.1]); %qcase2
  set(gca,'YTick',[-0.1 0 0.1]); %qcase2
  %set(gca,'YTick',[-0.1 0 0.1]); %qcase1
  
  legend('Original','Filtered') %intial
  
  set(gca,'Fontsize',10,'Fontname','Times New Roman');
%   set(gca,'YLim',[-0.1 0.1]);
%   set(gca,'YTick',[-0.1 0 0.1]);