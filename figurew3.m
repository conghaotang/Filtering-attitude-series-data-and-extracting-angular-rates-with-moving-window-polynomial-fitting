function  figurew3(wtrue,w_mea,w_est,m,n)

%½ÇËÙÂÊµÄÎó²îÍ¼
  epochs=n:m;
  figure
  set(gcf,'color',[1,1,1])
  
  subplot(3,1,1);
  plot(epochs, wtrue(1,:),'Linewidth',2.5);
  hold on
  plot(epochs, w_mea(1,:),'Linewidth',2.5);   
  hold on
  plot(epochs,w_est(1,:),'Linewidth',2.5);
 
  subplot(3,1,2);
  plot(epochs, wtrue(2,:),'Linewidth',2.5);
  hold on
  plot(epochs, w_mea(2,:),'Linewidth',2.5); %ylabel('Pitch'); 
  hold on
  plot(epochs,w_est(2,:),'Linewidth',2.5);
  
  subplot(3,1,3);
  plot(epochs, wtrue(3,:),'Linewidth',2.5);
  hold on
  plot(epochs, w_mea(3,:),'Linewidth',2.5); xlabel('Epoch [0.1s]');  %ylabel('Roll'); %,'Linewidth',1 ,'Linewidth',1
  hold on
  plot(epochs,w_est(3,:),'Linewidth',2.5);

  legend('True','Difference','Proposed')