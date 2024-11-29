function figureWpsd(psd_Aest,psd_Aorig,psd_Aest_f)
  
  f = psd_Aest_f;
  
  P11=(psd_Aest(:,1)); %10*log10 sqrt
  P21=(psd_Aorig(:,1));
 
  P12=(psd_Aest(:,2)); %10*log10
  P22=(psd_Aorig(:,2));

  P13= (psd_Aest(:,3)); %10*log10 sqrt  sqrt
  P23= (psd_Aorig(:,3));
  
  figure
  subplot(3,1,1)
     set(gca,'Fontsize',10,'Fontname','Times New Roman');
     semilogy(f, (P21),'Linewidth',1.5); %semilogy(f, P1); % %semilogy(f, P1);  %(200:end)
     hold on;
     semilogy(f,(P11),'Linewidth',1.5); %semilogy(f, P2); % %semilogy(f, P2);
      set(gca,'YLim',[10^-20 10^10]);
subplot(3,1,2)
      set(gca,'Fontsize',10,'Fontname','Times New Roman');
      semilogy(f, (P22),'Linewidth',1.5); ylabel('PSD [бу^2/s]');
      %      ylabel('$$ \sqrt{PSD} $$','Interpreter','latex');
      
      hold on;
      semilogy(f,(P12),'Linewidth',1.5); % arcsec/sqrt(Hz)
      %      legend('Filtered','Original') (rad/s)^2/Hz   rad/s/sqrt(Hz)
       set(gca,'YLim',[10^-20 10^10]);
 subplot(3,1,3)
    set(gca,'Fontsize',10,'Fontname','Times New Roman');
    semilogy(f, (P23),'Linewidth',1.5); xlabel('Frequency [1Hz]');
    hold on;
    semilogy(f,(P13),'Linewidth',1.5);
    set(gca,'YLim',[10^-10 10^10]);
    set(gca,'YTick',[10^-10 10^0 10^10]);
    legend( 'Difference','Proposed')