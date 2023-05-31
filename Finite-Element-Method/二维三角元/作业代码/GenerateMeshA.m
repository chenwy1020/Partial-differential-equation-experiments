%% �˳�������չʾ�����matlab���������񣬲��õ� [p e t]
% ͨ���˷�ʽ���Ը����߽���
% ������� [p, e, t] �μ� "pdetriples.pdf"
% ����ͼ������ pdeplot() �� pdemesh(), ����μ� "pdeplot.pdf", "pdemesh.pdf"

clear
%% ����ֻ���������� [-0.5, 0.5]*[1, 2]
% ����������ɻ�����״�����ο�"2DGeometryCreation.pdf"
R1=[3, 4, -0.5, 0.5, 0.5, -0.5, 1 , 1 , 2, 2]';
gm=[R1];
g=decsg(gm);

%% ��ʽһ�� �� initmesh ���� [p, e, t] ���� pdemesh ��ͼ
% ������� initmesh �ο� "initmesh.pdf"
[p, e, t]=initmesh(g);
pdemesh(p, e, t);

%% ��ʽ1.1���� refinemesh ��������
figure
% ������� refinemesh() �ο� "refinemesh.pdf"
[p, e, t]=refinemesh(g, p, e, t);
pdeplot(p, e, t);


%% ��ʽ������ geometryFromEdges������������� generateMesh, meshtoPet
figure
model=createpde;

% ������� geometryFromEdges() �ο� "GeometryFromEdge.pdf"
% ��������Ϣ�󶨵�model
geometryFromEdges(model,g);  %
pdegplot(model,'EdgeLabels','on')

% ������� generateMesh() �ο� "generatemesh.pdf"
% ��genereateMesh ��������, ���и��ߵ����ɶ�
% ע��GeometricOrder ������ѡ�'quadratic' �� 'linear'
%        'linear' �����������ζ���
%        'quadratic' ���������ζ���+�е�
% ע��'Hmax' Ϊ���߳�
% mesh=generateMesh(model, 'GeometricOrder', 'linear', 'Hmax', 0.1);
mesh=generateMesh(model, 'GeometricOrder', 'quadratic', 'Hmax', 0.2);
figure
% ����ת��Ϊ p e t 
[p, e, t]=meshToPet(mesh);
pdemesh(p, e, t, 'NodeLabels', 'On');

