%%%% Modelo %%%%

n = 96;
p = 3962/(96*96);
p= [1-p,p];
N = 0:1;
MA = zeros(n,n);
for i=1:n
    MA(i,:) = randsample(N, 96, true, p);
    MA(i,i) = 0;
end

G = graph(MA, 'upper');
graus = degree(G)';
% 2 - Path
path = distances(G);
diameter = max(max(path));
nzele = nnz(path);
pathmedium = sum(sum(path))/nzele;

% 3 - Clusterization
neigh = {};
for i=1:n
    neigh{i} = neighbors(G,i)';
end
for i=1:n
    if size(neigh{i},2) ~= 0
        for j=1:size(neigh{i},2)
            inter(i,j) = size(intersect(neigh{i},neigh{neigh{i}(1,j)}),2);
        end
    else
        inter(i,1) = 0;
    end
end
%C_i = 2*E_i/(degree_i(degree_i-1))
for i=1:n
    E(i,1) = max(inter(i,:));
    C(i,1) = 2*E(i,1)/(graus(1,i)*(graus(1,i)-1));
end