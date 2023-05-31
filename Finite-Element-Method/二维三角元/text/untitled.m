clear
clc

model=createpde;
R1=[3,4,0,1,1,0,0,0,1,1]';%[ a b x[4] y[4] ]
gn=[R1];
g=decsg(gn);
pdegplot(g,"EdgeLabels","on","FaceLabels","on")
xlim([-5,5]);
ylim([-5,5]);
axis equal
geometryFromEdges(model,g);

mesh=generateMesh(model,'GeometricOrder','linear','Hmax',0.2);
[p, e, t ]= meshToPet(mesh);
figure
pdemesh(p, e, t );


Np=length(p);
A=zeros(Np,Np);
b=zeros(Np,1);

x=p(1,:) ;
y=p(2,:) ;
numbers=[1:1:length(p)];
scatter(x,y);
text(x, y,string(numbers),'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom'); % 添加标签
xlim([-2,2]);
ylim([-2,2]);
axis equal


