clear all;
close all;
clc;
fname='al_jazeera_r1.wav';
[Speech, Fs]=audioread(fname);
[Breath_seg] = breath_detection(Speech,Fs);
plot(0:1/Fs:(length(Speech)-1)/Fs,Speech,'b',...
     0:1/Fs:(length(Breath_seg)-1)/Fs,0.5*Breath_seg,'r')
 title( 'Breath Detection in speech');
 xlabel('Time');
 ylabel('Amplitude');
 legend('Speech','Breath locations');
