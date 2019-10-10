clear all
close all
clc

t = 0:0.025:2.6;

x = cos(2*t)+sin(4*t);

P = 5;

dt = t(2)-t(1);
N=length(x); 
Df=1/(N*dt);
h=hanning(N);

X=fft(h'.*x)/N;
X=X(1:N/2);

Df=1/(N*dt);
Freq=((1:N/2)-1)*Df;

figure
stem(Freq,abs(X))

xIndex = find(X == max(X), 1, 'first');
Freqmax = Freq(xIndex);

a = find(X>0, 1, 'first');
if a == 1
    a = 1+1;
else a == a
end

Freq1 = Freq(a);



T = 8*pi/Freqmax;
Tm=T;

e = normrnd(0,0.05,[1,length(t)]);
x = x+e;

figure
plot(t,x)
hold on


ee = 100000;
ee2 = 100000;
e = 100000;

num = 0;

while ee2>2
    
T = normrnd(Tm,40,[1,length(t)]);


for l=1:1:length(t)
    
W = zeros(length(t),P*2+1);
    
for n=1:1:length(t)
   for k=1:1:P
       W(n,k)=exp(1i*2*pi*k*t(n)/T(l));
   end
end

for n=1:1:length(t)
   for k=P+1:1:P+1
       W(n,k)=exp(-1i*2*pi*0*t(n)/T(l));
   end
end

for n=1:1:length(t)
   for k=P+2:1:(2*P+1)
       W(n,k)=exp(-1i*2*pi*k*t(n)/T(l));
   end
end

Xhat = inv(conj(W')*W)*conj(W')*x';

xhat = W*Xhat;

for i=1:1:length(x)
    e=sum((xhat(i)-x(i))^2)/length(x);
end

   if abs(e)<abs(ee)
    ee=abs(e);
    Tm = T(l);
    ee2 = sum((abs(x')-abs(xhat)))/length(x);
   else
    ee=ee; 
   end

   num = num+1;
   
end

T = Tm;

end

hold on
plot(t,(xhat))
