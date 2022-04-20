clear; clc; close all;
%%% Get all statistics and plots from matriz adjacency.
cd ('Dados')
x = importdata('matrix.mat'); % matrix adjacency. 
nnodes = size(x,1);
G = graph(x); % build graph from x

%% Statistics
% 1 - Degree - total number of edges connected to a node
for i=1:nnodes
    degree(i) = nnz(x(i,:)); %nnz(x(i,:)) -> in case of not pondered edges
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
k = [0:1:96];
pk = zeros(97,2);
for i=1:97
    pk(i,1) = i-1;
    if ismember(k(i),aux)
        m = find(aux == k(i));
        pk(i,2) = f(m,2);
    else
        pk(i,2) = 0;
    end
end
       
% 1.3 Degree distribution
figure(1);
histogram(degree,30,'Normalization','probability');
title('Degree Distribution');
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

%5 - Componente Gigante e Conexividade
bins = conncomp(G);j=1;
for i=1:size(bins,2)
    um = find(bins==i);
    if size(um,2) > 1
        componentes (j,1) = size(um,2);
        componentes (j,2) = i;
        j=j+1;
    end
end
% 6 - Modularidade
% Feito no R, por motivos. 
%% END OF STATISTICS

%%%%%%%%%%%%%%%%% Plotting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
labels = importdata('AdjM.mat'); 
color = zeros(nnodes,3); linescolor = zeros(nnodes,3); 
A = 10; xp = zeros(nnodes,1); yp = zeros(nnodes,1);
B = 18;
%% Ordening nodes considering the degree for change position
% 1 - Ordered by Detective and Degree
flag = 2;
if flag == 1
    auxnew = []; newposition = []; j=1;
    detectiveposition = [4;8;57;80;83;87;nnodes];
    for i=1:size(detectiveposition,1)
        [auxnew, order] = sort(degree(1,j:detectiveposition(i,1)));
        alpha = ceil(size(auxnew,2)/2);
        M = [auxnew; order];
        auxp1 = M(:,alpha:size(auxnew,2)); auxp2 = flip(M(:,1:alpha-1),2);
        newp = [auxp1 auxp2];
        newposition = [newposition newp];
        j = detectiveposition(i,1)+1;
    end
    % Node Colors and Positions
    a = 3*A/4; b = A; c = 3*B/5; d = B;
    [auxxp, auxyp] = xydata(4,c,d,a,b);
    for i=1:4
        color(i,1) = 178/255;
        color(i,2) = 186/255;
        color(i,3) = 187/255;
        xp(newposition(2,i),1) = auxxp(i,1);
        yp(newposition(2,i),1) = auxyp(i,1);
    end
    a = 3*A/4; b = A; c = 0; d = B/5;
    [auxxp, auxyp] = xydata(4,c,d,a,b); j=1;
    for i=5:8
        color(i,1) = 240/255;
        color(i,2) = 178/255;
        color(i,3) = 122/255;
        xp(newposition(2,i)+4,1) = auxxp(j,1);
        yp(newposition(2,i)+4,1) = auxyp(j,1);
        j=j+1;
    end
    a = 0; b =3*A/4; c = 0; d = 3*B/5;
    [auxxp, auxyp] = xydata(49,c,d,a,b); j=1;
    for i=9:57
        color(i,1) = 130/255;
        color(i,2) = 224/255;
        color(i,3) = 170/255;
        xp(newposition(2,i)+8,1) = auxxp(j,1);
        yp(newposition(2,i)+8,1) = auxyp(j,1);
        j=j+1;
    end
    a = A/4; b =3*A/4; c = 3*B/5; d = B; j=1;
    [auxxp, auxyp] = xydata(23,c,d,a,b);
    for i=58:80
        color(i,1) = 187/255;
        color(i,2) = 143/255;
        color(i,3) = 206/255;
        xp(newposition(2,i)+57,1) = auxxp(j,1);
        yp(newposition(2,i)+57,1) = auxyp(j,1);
        j=j+1;
    end
    a = 0; b =A/5; c = 3*B/5; d = 4*B/5;
    [auxxp, auxyp] = xydata(3,c,d,a,b); j=1;
    for i=81:83
        color(i,1) = 236/255;
        color(i,2) = 112/255;
        color(i,3) = 99/255;
        xp(newposition(2,i)+80,1) = auxxp(j,1);
        yp(newposition(2,i)+80,1) = auxyp(j,1);
        j=j+1;
    end
    a = 3*A/4; b =A; c = 9*B/20; d = 3*B/5;
    [auxxp, auxyp] = xydata(4,c,d,a,b); j=1;
    for i=84:87
        color(i,1) = 154/255;
        color(i,2) = 125/255;
        color(i,3) = 10/255;
        xp(newposition(2,i)+83,1) = auxxp(j,1);
        yp(newposition(2,i)+83,1) = auxyp(j,1);
        j=j+1;
    end
    a = 3*A/4; b =A; c = 2*B/10; d = 9*B/20;
    [auxxp, auxyp] = xydata(9,c,d,a,b); j=1;
    for i=88:nnodes
        color(i,1) = 137/255;
        color(i,2) = 232/255;
        color(i,3) = 242/255;
        xp(newposition(2,i)+87,1) = auxxp(j,1);
        yp(newposition(2,i)+87,1) = auxyp(j,1);
        j=j+1;
    end
end
% 2 - Ordered by degree
if flag == 2
    newposition = [];
    [auxnew, order] = sort(degree(1,1:nnodes));
    alpha = ceil(size(auxnew,2)/2);
    M = [auxnew; order];
    k = size(M,2); newposition(:,alpha) = M(:,k); j=1;
    for i=1:1:k/2-1
        newposition(:,alpha+i) = M(:,k-j);
        j = j+1;
        newposition(:,alpha-i) = M(:,k-j);
        j = j+1;
    end
    newposition(:,alpha+48) = M(:,k-j);
    c = 0; d =A; a =0; b= A;
    [auxxp, auxyp] = xydata(nnodes,c,d,a,b);
    for j=1:nnodes
        xp(newposition(2,j),1) = auxxp(j,1);
        yp(newposition(2,j),1) = auxyp(j,1);
    end
end
if flag == 2 || flag == 3
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
end

% 3 - Ordered by Comunitty and degree
if flag == 3
    commu1 = importdata('commu1.csv');
    commu2 = importdata('commu2.csv');
    commu3 = importdata('commu3.csv');
    t1 = size(commu1,2); t2 = size(commu2,2); t3 = size(commu3,2);
    for i=1:t1
        auxd(i,1) = commu1(1,i);
        auxd(i,2) = degree(1,auxd(i,1));
    end
    [~,order] = sort(auxd(:,2)); % sort just the first column
    M = auxd(order,:)'; 
    newposition = [];
    alpha = ceil(size(M,2)/2);
    k = size(M,2); newposition(:,alpha) = M(:,k); j=1;
    for i=1:1:k/2-1
        newposition(:,alpha+i) = M(:,k-j);
        j = j+1;
        newposition(:,alpha-i) = M(:,k-j);
        j = j+1;
    end
    newposition(:,t1) = M(:,1);
    a = 2/10*A; b =A; c = 0; d = 2/5*B;
    [auxxp, auxyp] = xydata(t1,c,d,a,b);
    for j=1:40
        xp(newposition(1,j),1) = auxxp(j,1);
        yp(newposition(1,j),1) = auxyp(j,1);
    end
    auxd = [];
    for i=1:t2
        auxd(i,1) = commu2(1,i);
        auxd(i,2) = degree(1,auxd(i,1));
    end
    [~,order] = sort(auxd(:,2)); % sort just the first column
    M = auxd(order,:)'; 
    newposition = [];
    alpha = ceil(size(M,2)/2);
    k = size(M,2); newposition(:,alpha) = M(:,k); j=1;
    for i=1:1:k/2-1
        newposition(:,alpha+i) = M(:,k-j);
        j = j+1;
        newposition(:,alpha-i) = M(:,k-j);
        j = j+1;
    end
    newposition(:,t2) = M(:,1);
    a = 0; b =1/5*A; c = 30/100*B; d = 45/100*B;
    [auxxp, auxyp] = xydata(t2,c,d,a,b);
    for j=1:t2
        xp(newposition(1,j),1) = auxxp(j,1);
        yp(newposition(1,j),1) = auxyp(j,1);
    end
    
    auxd = [];
    for i=1:t3
        auxd(i,1) = commu3(1,i);
        auxd(i,2) = degree(1,auxd(i,1));
    end
    [~,order] = sort(auxd(:,2)); % sort just the first column
    M = auxd(order,:)'; 
    newposition = [];
    alpha = ceil(size(M,2)/2);
    k = size(M,2); newposition(:,alpha) = M(:,k); j=1;
    for i=1:1:k/2-1
        newposition(:,alpha+i) = M(:,k-j);
        j = j+1;
        newposition(:,alpha-i) = M(:,k-j);
        j = j+1;
    end
    newposition(:,t3) = M(:,1);
    a = 0; b =A; c = 9/20*B; d = B;
    [auxxp, auxyp] = xydata(t3,c,d,a,b);
    for j=1:t3
        xp(newposition(1,j),1) = auxxp(j,1);
        yp(newposition(1,j),1) = auxyp(j,1);
    end
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
comcolor1(1,1) = 195/255;
comcolor1(1,2) = 195/255;
comcolor1(1,3) = 162/255;
comcolor2(1,1) = 235/255;
comcolor2(1,2) = 235/255;
comcolor2(1,3) = 224/255;
comcolor3(1,1) = 109/255;
comcolor3(1,2) = 109/255;
comcolor3(1,3) = 70/255;

degree = degree +1; % desenhar todos os nós
figure(2); 

plot(G,'XData',xp,'YData',yp, 'NodeLabel', names, 'MarkerSize', 3*degree, 'NodeFontSize', degree, 'LineWidth', G.Edges.Weight, 'NodeColor', color, 'EdgeColor', coloredge, 'EdgeAlpha',0.1);
set(gca,'xtick',[]); set(gca,'ytick',[]); set(gcf,'PaperUnits','inches','PaperPosition',[0 0 20 20]); axis off; set(gca,'LooseInset',get(gca,'TightInset'));
if flag == 3
    rectangle('Position',[30/100*B 0 17/100*B 9.8/50*A],'Curvature',0.2, 'EdgeColor', comcolor2, 'LineWidth',3); hold on
    rectangle('Position',[97/200*B 0 107/200*B 102/100*A],'Curvature',0.2, 'EdgeColor', comcolor3, 'LineWidth',3); hold on
    rectangle('Position',[0 2/10*A 45/100*B 8/10*A],'Curvature',0.2, 'EdgeColor', comcolor1, 'LineWidth',3); hold on
end
cd ..
