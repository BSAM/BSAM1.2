function[] = surfplot3d(s0,nn,var,toggle1,toggle2)
%
% s0 is the string indicating the appropriate PROBLEM directory, eg, 'CH'
% nn is the iteration number to print
% var is the variable number
% Set toggle1 = 0 for no output
% Set toggle1 = 1 for jpg output
% Set toggle1 = 2 for eps output
% Set toggle2 = 0 to print only the base-level bounding box
% Set toggle2 = 1 to print the bounding boxes
% Set toggle2 to any other value to turn off all bounding boxes

s1 = ['0000000' num2str(nn)];
s2 = s1((length(s1)-4):length(s1));

dir =['../' s0 '/OUT/']

IN  = [dir 'm' s2 '.dat']

if toggle1 == 1
  OUT = [dir 'm' num2str(toggle2) s2 '.jpg']
elseif toggle1 == 2
  OUT = [dir 'm' num2str(toggle2) s2 '.eps']
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
    [dx(3),count] =  fscanf(f, '%f', 1);

    [xl(1),count] =  fscanf(f, '%f', 1);
    [xl(2),count] =  fscanf(f, '%f', 1);
    [xl(3),count] =  fscanf(f, '%f', 1);

    [xu(1),count] =  fscanf(f, '%f', 1);
    [xu(2),count] =  fscanf(f, '%f', 1);
    [xu(3),count] =  fscanf(f, '%f', 1);

    [n(1),count] =  fscanf(f, '%d', 1);
    [n(2),count] =  fscanf(f, '%d', 1);
    [n(3),count] =  fscanf(f, '%d', 1);
    
    [mg(1,1),count] =  fscanf(f, '%d', 1);
    [mg(1,2),count] =  fscanf(f, '%d', 1);
    [mg(2,1),count] =  fscanf(f, '%d', 1);
    [mg(2,2),count] =  fscanf(f, '%d', 1);
    [mg(3,1),count] =  fscanf(f, '%d', 1);
    [mg(3,2),count] =  fscanf(f, '%d', 1);

    ipatch = ipatch + 1;
    disp(sprintf('Processing patch number %d', ipatch));
    disp(sprintf('The patch size is %d,  %d, %d', n(1), n(2), n(3)));

    xu = xl+dx.*n;
    if ((toggle2 == 1) && (level>0))
      drawbox(xl,xu,level);
    end;

    if level == 0
      xlg(1) = xl(1);
      xug(1) = xu(1);
      xlg(2) = xl(2);
      xug(2) = xu(2);
      xlg(3) = xl(3);
      xug(3) = xu(3);
    end;

    A = zeros(nrvars,(n(1)+2)*(n(2)+2)*(n(3)+2),'double');

    [A]=fscanf(f,'%f', [nrvars,(n(1)+2)*(n(2)+2)*(n(3)+2)]); % ghost layer included.
   
    if level == maxlevel
      drawsurf(A,var,n,xl,xu,dx,xlg,xug);
    end;
  else
    theend = 1;
  end
end;

%view(115,20)
%view(3)
view(160,20)
%camlight headlight; material dull; lighting phong
%camlight(115,30); camlight(-65,-30)
camlight headlight;
material dull; lighting gouraud

t = (xug-xlg)./4;
text(xlg(1),xlg(2),xug(3)+t(3)/4,['time = ' num2str(time)], 'FontSize', 14);

grid on
set(gca,'XTick',xlg(1):t(1):xug(1))
set(gca,'YTick',xlg(2):t(2):xug(2))
set(gca,'ZTick',xlg(3):t(3):xug(3))

text(xlg(1)+2*t(1),xug(2)+1.5*t(2),xlg(3),'X', 'FontSize', 14);
text(xug(1)+t(1),xlg(2)+2*t(2),xlg(3),'Y', 'FontSize', 14);
text(xug(1)+.75*t(1),xlg(2),xlg(3)+2*t(3),'Z', 'FontSize', 14);

if toggle1 == 1
  print('-djpeg','-r400',OUT)
elseif toggle1 == 2
  print('-deps','-r400',OUT)
end;

ipatch
fclose(f);

function drawsurf(A,var,n,xl,xu,dx,xlg,xug)

%c = zeros(n(1)+2,n(2)+2,n(3)+2,'double');
c = zeros(n(2)+2,n(1)+2,n(3)+2,'double');

for k=1:n(3)+2
  for j=1:n(2)+2
    for i=1:n(1)+2

%      c(i,j,k) = A(1,((k-1)*(n(2)+2)+j-1)*(n(1)+2)+i);
      c(j,i,k) = A(var,((k-1)*(n(2)+2)+j-1)*(n(1)+2)+i);

    end;
  end;
end;

[xx,yy,zz]=meshgrid(linspace(xl(1)-dx(1)/2,xu(1)+dx(1)/2,n(1)+2), ...
                    linspace(xl(2)-dx(2)/2,xu(2)+dx(2)/2,n(2)+2), ...
                    linspace(xl(3)-dx(3)/2,xu(3)+dx(3)/2,n(3)+2));

p = patch(isosurface(xx, yy, zz, c, 0.5));
%isonormals(xx,yy,zz,c,p) % This fails to calculate the true normals.
set(p, 'FaceColor', [0.5,0.5,0.5], 'EdgeColor', 'none');
daspect([1 1 1])

axis([xlg(1),xug(1),xlg(2),xug(2),xlg(3),xug(3)])
axis equal
axis([xlg(1),xug(1),xlg(2),xug(2),xlg(3),xug(3)])

function[] = drawbox(xl,xu,lvl)

switch lvl
  case 0
    linecolor = 'black';
    linethick = 2;
  case 1
    linecolor = 'blue';
    linethick = 1;
  case 2
    linecolor = 'red';
    linethick = 0.5;
  case 3
    linecolor = 'green';
    linethick = 0.25;
  case 4
    linecolor = 'cyan';
    linethick = 0.25;
  otherwise
    disp(sprintf('No color for this level: %d', level));
    stop;
end;

line([xl(1),xu(1)],[xl(2),xl(2)],[xu(3),xu(3)],'Color',linecolor,'LineWidth',linethick)
line([xl(1),xu(1)],[xu(2),xu(2)],[xu(3),xu(3)],'Color',linecolor,'LineWidth',linethick)
line([xl(1),xl(1)],[xl(2),xu(2)],[xu(3),xu(3)],'Color',linecolor,'LineWidth',linethick)
line([xu(1),xu(1)],[xl(2),xu(2)],[xu(3),xu(3)],'Color',linecolor,'LineWidth',linethick)

line([xl(1),xu(1)],[xl(2),xl(2)],[xl(3),xl(3)],'Color',linecolor,'LineWidth',linethick)
line([xl(1),xu(1)],[xu(2),xu(2)],[xl(3),xl(3)],'Color',linecolor,'LineWidth',linethick)
line([xl(1),xl(1)],[xl(2),xu(2)],[xl(3),xl(3)],'Color',linecolor,'LineWidth',linethick)
line([xu(1),xu(1)],[xl(2),xu(2)],[xl(3),xl(3)],'Color',linecolor,'LineWidth',linethick)

line([xl(1),xl(1)],[xl(2),xl(2)],[xl(3),xu(3)],'Color',linecolor,'LineWidth',linethick)
line([xl(1),xl(1)],[xu(2),xu(2)],[xl(3),xu(3)],'Color',linecolor,'LineWidth',linethick)
line([xu(1),xu(1)],[xl(2),xl(2)],[xl(3),xu(3)],'Color',linecolor,'LineWidth',linethick)
line([xu(1),xu(1)],[xu(2),xu(2)],[xl(3),xu(3)],'Color',linecolor,'LineWidth',linethick)
