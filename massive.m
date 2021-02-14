X0=0;
Y0=-2;
XEnd=8;
YEnd=4;
SizeX=8;
SizeY=6;
hx=0.1;
hy=0.1;

% xv=[0 0 4 4 8 8 3 3 0];
% yv=[-1 3.5 3.5 2 2 -1.5 -1.5 -1 -1];

xv=[0 0 8 8 0 NaN 0 0 8 8 0];
yv=[-1 1 1 -1 -1 NaN 2 3 3 2 2];

[Xx,Yy] = ndgrid(X0:hx:XEnd, Y0:hy:YEnd);
[xc,yc]=ndgrid(X0+hx/2:hx:XEnd-hx/2,Y0+hy/2:hy:YEnd-hy/2);
[n,m]=size(Xx);
%plot(xc,yc);

in = inpolygon(xc,yc,xv,yv);
plot(xv,yv,'LineWidth',2);
axis equal

hold on
plot(xc(in),yc(in),'bo');
plot(xc(~in),yc(~in),'r+');

A=rot90(in);
%for i=1:m-1
%    for j=1:n-1
%        if A(i,j)==0
%          A(i,j)=1; 
%       else
%            A(i,j)=0;
%        end
%    end
%end

dlmwrite('tubes.txt',A,'delimiter',' ');