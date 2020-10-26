X = [0, 0; 0,4; 4,4; 4,2; 6,2; 6,0; 8,0; 8,-2; 5,-2; 5,1; 2,1; 2,0];
pg = polyshape(X(:,1),X(:,2));
plot(pg);
xv=X(:,1);
yv=X(:,2);
X0=0;
Y0=-2;
SizeX=8;
SizeY=6;
[Xx,Yy] = ndgrid(0:.2:8,-2:.2:4);
[xc,yc]=ndgrid(0.1:0.2:7.9,-1.9:.2:3.9);
[n,m]=size(Xx);

in = inpolygon(xc,yc,xv,yv);
plot(pg);
axis equal

hold on
plot(xc(in),yc(in),'r+');
plot(xc(~in),yc(~in),'bo');

file = fopen('file.txt', 'a+');
if file == -1
    error('File is not open:(');
end




