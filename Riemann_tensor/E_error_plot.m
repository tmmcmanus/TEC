function[] = E_error_plot(h11_central,E,i,j,E_absolute_error)
h1  = subplot(2,1,1)
x=1:length(h11_central);
plot(x,E,'r',x,h11_central,'b-')
xlabel('index','FontSize',11,'Fontweight','bold')
xlim([1 length(h11_central)])
legend('E','E_c_a_l_c')
grid on
title(sprintf('E and E_c_a_l_c (%dx%d Atlas)',i-3,j-3),'FontSize',13,'Fontweight','bold')
h2 = subplot(2,1,2)
plot(E_absolute_error)
xlabel('index','FontSize',11,'Fontweight','bold')
xlim([1 length(h11_central)])
title(sprintf('E Absolute Error (%dx%d Atlas)',i-3,j-3),'FontSize',13,'Fontweight','bold')
grid on
end
