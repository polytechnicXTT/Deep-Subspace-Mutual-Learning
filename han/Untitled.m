a=rand(1e4);
b=rand(1e4);
tic,c = a*b;
toc
tic,d=sin(a(:));
toc
tic,d = fft(d);
toc