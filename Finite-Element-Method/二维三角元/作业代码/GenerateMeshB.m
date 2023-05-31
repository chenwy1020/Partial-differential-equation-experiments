%% �˳�������չʾ�����matlab���������񣬲��õ� [p e t]
% ͨ���˷�ʽ���Ը����߽���
% ������� [p, e, t] �μ� "pdetriples.pdf"
% ����ͼ������ pdeplot() �� pdemesh(), ����μ� "pdeplot.pdf", "pdemesh.pdf"

clear
%% ����ֻ���������� [-0.5, 0.5]*[1, 2] �۳�Բ B_{0.3}(0, 1.5)
R1=[3, 4, -0.5, 0.5, 0.5, -0.5, 1 , 1 , 2, 2]';
C1=[1, 0, 1.5, 0.3]';
C1=[C1; zeros(length(R1)-length(C1), 1)];
gm=[R1, C1];
sf='R1-C1';

ns=char('R1', 'C1');
ns=ns';
g=decsg(gm, sf, ns);

%% ��ʽһ�� �� initmesh ���� [p, e, t] ���� pdemesh ��ͼ
% ������� initmesh �ο� "initmesh.pdf"
[p, e, t]=initmesh(g);
pdemesh(p, e, t);

%% ��ʽ1.1���� refinemesh ��������
figure
% ������� refinemesh() �ο� "refinemesh.pdf"
[p, e, t]=refinemesh(g, p, e, t);
pdemesh(p, e, t);


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
mesh=generateMesh(model, 'GeometricOrder', 'linear', 'Hmax', 0.1);

figure
% ����ת��Ϊ p e t 
[p, e, t]=meshToPet(mesh);
pdemesh(p, e, t, 'NodeLabels', 'On');

