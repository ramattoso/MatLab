clear; clc; close all;
%%% Get all statistics and plots from matriz adjacency.
cd ('Dados')
x = importdata('matrix.mat'); % matrix adjacency. 
nnodes = size(x,1);
G = graph(x); % build graph from x
%% Statistics
% 1 - Degree - total number of edges connected to a node
degree = zeros(1,nnodes);
for i=1:nnodes
    degree(i) = nnz(x(i,:));
end
% 1.1 mean degree 
mdegree = sum(degree)/nnodes;
% 1.2 relative degree frequency
% f(k) = n(k)/nnodes;
aux = unique(degree); % acha os valores únicos de grau
for i=1:size(aux,2)
    n = find(degree ==aux(i));
    f(i,1) = aux(i);
    f(i,2) = size(n,2)/nnodes;
end
% 1.3 Degree distribution
figure(1); scatter(f(:,1), f(:,2), 20,'filled', 'r'); title('Degree Distribution');
xlabel('$k$','interpreter','latex'); ylabel('$P_k$','interpreter','latex');

% 2 - Path
path = distances(G);
diameter = max(max(path));
nzele = nnz(path);
pathmedium = sum(sum(path))/nzele;

% 3 - Clusterization
neigh = {};
for i=1:nnodes
    neigh{i} = neighbors(G,i)';
end
for i=1:nnodes
    if size(neigh{i},2) ~= 0
        for j=1:size(neigh{i},2)
            inter(i,j) = size(intersect(neigh{i},neigh{neigh{i}(1,j)}),2);
        end
    else
        inter(i,1) = 0;
    end
end

%C_i = 2*E_i/(degree_i(degree_i-1))
for i=1:nnodes
    E(i,1) = max(inter(i,:));
    C(i,1) = 2*E(i,1)/(degree(1,i)*(degree(1,i)-1));
end

% 4 - Centrality
Centerg = centrality(G,'degree');
Centerc = centrality(G,'closeness');
Centerb = centrality(G,'betweenness');

%%% Ordening nodes considering the degree for change position
newposition = []; 
[auxnew, order] = sort(degree(1,1:nnodes));
alpha = ceil(size(auxnew,2)/2);
M = [auxnew; order];
% auxp1 = M(:,alpha:size(auxnew,2)); auxp2 = flip(M(:,1:alpha-1),2);
% newp = [auxp1 auxp2];
% newposition = [newposition newp];
k = size(M,2); newposition(:,alpha) = M(:,k); j=1;
for i=1:1:k/2-1
    newposition(:,alpha+i) = M(:,k-j);
    j = j+1;
    newposition(:,alpha-i) = M(:,k-j);
    j = j+1;
end
newposition(:,alpha+48) = M(:,k-j);
%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Node Colors and Positions
labels = importdata('AdjM.mat'); 
color = zeros(nnodes,3); linescolor = zeros(nnodes,3); 
A = 1000; xp = zeros(nnodes,1); yp = zeros(nnodes,1);
c = 0; d =A; a =0; b= A;
[auxxp, auxyp] = xydata(nnodes,c,d,a,b);

for i=1:4
    color(i,1) = 178/255;
    color(i,2) = 186/255;
    color(i,3) = 187/255;
end
for i=5:8
    color(i,1) = 240/255;
    color(i,2) = 178/255;
    color(i,3) = 122/255;
end
for i=9:57
    color(i,1) = 130/255;
    color(i,2) = 224/255;
    color(i,3) = 170/255;
end
for i=58:80
    color(i,1) = 187/255;
    color(i,2) = 143/255;
    color(i,3) = 206/255;
end
for i=81:83
    color(i,1) = 236/255;
    color(i,2) = 112/255;
    color(i,3) = 99/255;
end
for i=84:87
    color(i,1) = 154/255;
    color(i,2) = 125/255;
    color(i,3) = 10/255;
end
for i=88:nnodes
    color(i,1) = 137/255;
    color(i,2) = 232/255;
    color(i,3) = 242/255;
end
for j=1:nnodes
    xp(newposition(2,j),1) = auxxp(j,1);
    yp(newposition(2,j),1) = auxyp(j,1);
end
%% Edge Colors
j = 1; i=1; aux = table2array(G.Edges(j,1));
while j < size(G.Edges,1)
    while aux(1) == i
        coloredge(j,:) = color(i,:);
        j=j+1;
        aux = table2array(G.Edges(j,1));
    end
    i = i+1;
end
coloredge(j,:) = color(nnodes,:);
%% Nodes Names
for i=1:nnodes
    names{i} = char(splitlines(compose(labels{i,1})));
end
%% final plot
degree = degree +1; % desenhar todos os nós
figure(2); h = plot(G,'XData',xp,'YData',yp, 'NodeLabel',names,'MarkerSize', 2*degree/5,'NodeFontSize', degree/12, 'LineWidth', G.Edges.Weight, 'NodeColor', color, 'EdgeColor', coloredge, 'EdgeAlpha',0.1);
set(gca,'xtick',[]); set(gca,'ytick',[]); set(gcf,'PaperUnits','inches','PaperPosition',[0 0 20 20]); axis off; set(gca,'LooseInset',get(gca,'TightInset')); 
cd ..
