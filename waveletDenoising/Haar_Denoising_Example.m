% Uses the Haar wavelet transform with local false discovery rate control
% to denoise signals contaminated with white Gaussian noise

y = Make_Signal('Bumps',10^4);
y(5000:50:6000) = y(5000:50:6000)+30;
x = linspace(0,1,length(y));
yn = y+4*randn(size(y));

yd1 = Haar_Denoise_LFDR(yn,0.5);
yd2 = Haar_Denoise_LFDR(yn,0.01);

subplot(1,2,1)
plot(x,yn,'Color',[0.5,0.8,1]);
hold on
plot(x,yd1,'Color',[1,0,0]);
hold off
subplot(1,2,2)
plot(x,yn,'Color',[0.5,0.8,1]);
hold on
plot(x,yd2,'Color',[1,0,0]);
hold off