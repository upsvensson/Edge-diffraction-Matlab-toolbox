tftot64 = EDres64.tftot; 

tftot96 = EDres96.tftot; 

figure(2)
clf(2)
h = semilogx(controlparameters.frequencies,20*log10(abs([tftot64 tftot96 tftot96-tftot64])),'-o'); 
g = get(h(1),'Parent');
set(g,'FontSize',14);
set(h(1),'LineWidth',2);
g = xlabel('Frequency [Hz]');
set(g,'FontSize',14)
g = ylabel('TF magnitude re. 1m free-field [dB]');
set(g,'FontSize',14)
g = title('Freq. resp. of a noise barrier on a rigid ground'); 
set(g,'FontSize',14)
xlim([50 5000])
grid