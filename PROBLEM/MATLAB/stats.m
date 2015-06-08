function[] = stats(n1)

IN= '../CHGSTest1/OUT/nfncstats.1.dat'

figure(1);
clf;
figure(2);
clf;
theend = logical(0);
f = fopen(IN,'r');
hold on;

[A]=fscanf(f,'%f', [4,n1]);

fclose(f);

for i=1:n1
  ta(i)   = A(1,i);
  ma(i)   = A(2,i);
  ea(i)   = A(3,i);
  mesh(i) = A(4,i);
end;

format long

maxmass = max(ma)
minmass = min(ma)
trumass = ma(1)

absmassdiff = maxmass-minmass

figure(1);
plot(ta,ma-ma(1))

hold on;

figure(2)
plot(ta,ea)

hold on;

figure(3)
plot(ta,mesh)

avemeshsize = sum(mesh(1:n1))/n1

hold on;


