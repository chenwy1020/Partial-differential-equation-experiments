%% 此程序用以展示如何在matlab上生成网格，并得到 [p e t]
% 通过此方式可以给定边界标号
% 具体关于 [p, e, t] 参见 "pdetriples.pdf"
% 网格画图可以用 pdeplot() 或 pdemesh(), 具体参见 "pdeplot.pdf", "pdemesh.pdf"

clear
%% 例：只生成正方形 [-0.5, 0.5]*[1, 2]
% 关于如何生成基本形状函数参看"2DGeometryCreation.pdf"
R1=[3, 4, -0.5, 0.5, 0.5, -0.5, 1 , 1 , 2, 2]';
gm=[R1];
g=decsg(gm);

%% 方式一： 由 initmesh 生成 [p, e, t] 并由 pdemesh 画图
% 具体关于 initmesh 参看 "initmesh.pdf"
[p, e, t]=initmesh(g);
pdemesh(p, e, t);

%% 方式1.1：由 refinemesh 加密网格
figure
% 具体关于 refinemesh() 参看 "refinemesh.pdf"
[p, e, t]=refinemesh(g, p, e, t);
pdeplot(p, e, t);


%% 方式二：由 geometryFromEdges（）生成网格和 generateMesh, meshtoPet
figure
model=createpde;

% 具体关于 geometryFromEdges() 参看 "GeometryFromEdge.pdf"
% 将几何信息绑定到model
geometryFromEdges(model,g);  %
pdegplot(model,'EdgeLabels','on')

% 具体关于 generateMesh() 参看 "generatemesh.pdf"
% 由genereateMesh 生成网格, 其有更高的自由度
% 注：GeometricOrder 有两个选项：'quadratic' 和 'linear'
%        'linear' 仅包含三角形顶点
%        'quadratic' 包含三角形顶点+中点
% 注：'Hmax' 为最大边长
% mesh=generateMesh(model, 'GeometricOrder', 'linear', 'Hmax', 0.1);
mesh=generateMesh(model, 'GeometricOrder', 'quadratic', 'Hmax', 0.2);
figure
% 网格转化为 p e t 
[p, e, t]=meshToPet(mesh);
pdemesh(p, e, t, 'NodeLabels', 'On');

