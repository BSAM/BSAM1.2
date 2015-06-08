function[] = cntrplot2d(s0,nn,var,toggle0,toggle1,toggle2)
%
% s0 is the string indicating the appropriate PROBLEM directory, eg, 'CH'
% nn is the frame number to be printed
% var is the variable number
% toggle0 = 0 for adaptive mesh
% toggle0 = 1 for uniform mesh
% toggle1 = 0 for no output
% toggle1 = 1 for jpg output
% toggle1 = 2 for eps output
% toggle2 = 0 for no bounding boxes
% toggle2 = 1 for bounding boxes

s1 = ['0000000' num2str(nn)];
s2 = s1((length(s1)-4):length(s1));

dir =['../' s0 '/OUT/']

if toggle0 == 0
  s3 = 'm'
else
  s3 = 'u'
end;

IN  = [dir s3 s2 '.dat']
if toggle1 == 1
  OUT = [dir s3 s2 '.jpg']
elseif toggle1 == 2
  OUT = [dir s3 s2 '.eps']
end;

clf reset;

set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.5 0.5 4 4]);

theend = logical(0);
f = fopen(IN,'r');
hold on;
ipatch=0;

[time,count] = fscanf(f, '%f', 1);
[maxlevel,count] = fscanf(f, '%d', 1);

while(~theend)

  [level,count] =  fscanf(f, '%d', 1);
  [ndim,count]  =  fscanf(f, '%d', 1);
  [r,count] =  fscanf(f, '%d', 1);
  [nrvars,count]  =  fscanf(f, '%d', 1);

  if count ~= 0

    [dx(1),count] =  fscanf(f, '%f', 1);
    [dx(2),count] =  fscanf(f, '%f', 1);

    [xl(1),count] =  fscanf(f, '%f', 1);
    [xl(2),count] =  fscanf(f, '%f', 1);

    [xu(1),count] =  fscanf(f, '%f', 1);
    [xu(2),count] =  fscanf(f, '%f', 1);

    [n(1),count] =  fscanf(f, '%d', 1);
    [n(2),count] =  fscanf(f, '%d', 1);
    
    [mg(1,1),count] =  fscanf(f, '%d', 1);
    [mg(1,2),count] =  fscanf(f, '%d', 1);
    [mg(2,1),count] =  fscanf(f, '%d', 1);
    [mg(2,2),count] =  fscanf(f, '%d', 1);

    ipatch = ipatch + 1;
    disp(sprintf('Processing patch number %d', ipatch));
    disp(sprintf('The patch size is %d,  %d', n(1), n(2)));

    xu = xl+dx.*n;
    if toggle2 == 1
      drawbox(xl,xu,level);
    end;
    
    if level == 0
      xlg(1) = xl(1);
      xug(1) = xu(1);
      xlg(2) = xl(2);
      xug(2) = xu(2);
    end;

    A = zeros(nrvars,(n(1)+2)*(n(2)+2),'double');
   
    [A]=fscanf(f,'%f', [nrvars,(n(1)+2)*(n(2)+2)]); % ghost layer included.
   
    if level == maxlevel
      drawcont(A,var,n,xl,xu,dx,xlg,xug);
    end;
  else
    theend = 1;
  end
end;

t = (xug-xlg)./4;
text(xlg(1)+t(1)/2,xug(2)+t(2)/4,['time = ' num2str(time)], 'FontSize', 14);

%grid on
set(gca,'XTick',xlg(1):t(1):xug(1))
set(gca,'YTick',xlg(2):t(2):xug(2))

if toggle1 == 1
  print('-djpeg','-r300',OUT)
elseif toggle1 == 2
  print('-deps','-r1200',OUT)
end;

ipatch
fclose(f);

function drawcont(A,var,n,xl,xu,dx,xlg,xug)

c = zeros(n(1)+2,n(2)+2,'double');
x = zeros(n(1)+2,n(2)+2,'double');
y = zeros(n(1)+2,n(2)+2,'double');

xll = xl(1)-dx(1)/2.0;
yll = xl(2)-dx(2)/2.0;

for j = 1:n(2)+2
  for i = 1:n(1)+2
    x(i,j) = xll+(i-1)*dx(1);
    y(i,j) = yll+(j-1)*dx(2);
    c(i,j) = A(var, (j-1)*(n(1)+2)+i);
  end;
end;

%colormap([0,0,0])
colormap('default')
%contour(x,y,c,0.5,'Linewidth',1)
%v = [ .1 .2 .3 .4 .5 .6 .7 .8 .9];
v = [-.1 0 .1];
%v = [ 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95];
%v = [0];
%v = [.1 .3 .5 .7 .9];
%contour(x,y,c,v)
%contourf(x,y,c)
%colorbar;
%contour(x,y,c,v)
%contour(x,y,c,[0.0,0.0])
%colormap([0,0,0])
%colormap('default')
surf(x,y,c,'LineStyle','none')
%contour(x,y,c,v,'Linewidth',1.25)
shading interp;
axis([xlg(1),xug(1),xlg(2),xug(2)])
axis equal
axis([xlg(1),xug(1),xlg(2),xug(2)])
%colorbar

function[] = drawbox(xl,xu,lvl)

switch lvl
  case 0
    linecolor = 'black';
    linethick = 4;
  case 1
    linecolor = 'blue';
    linethick = 2;
  case 2
    linecolor = 'red';
    linethick = 1;
  case 3
    linecolor = 'green';
    linethick = 0.5;
  case 4
    linecolor = 'cyan';
    linethick = 0.25;
  otherwise
    disp(sprintf('No color for this level: %d', level));
    stop;
end;


line([xl(1),xu(1)],[xl(2),xl(2)],'Color',linecolor,'LineWidth',linethick)
line([xl(1),xu(1)],[xu(2),xu(2)],'Color',linecolor,'LineWidth',linethick)
line([xl(1),xl(1)],[xl(2),xu(2)],'Color',linecolor,'LineWidth',linethick)
line([xu(1),xu(1)],[xl(2),xu(2)],'Color',linecolor,'LineWidth',linethick)




