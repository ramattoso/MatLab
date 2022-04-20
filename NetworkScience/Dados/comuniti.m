% Constructing comunities
%function [com, MA] = comuniti(x,nnodes, degree)
%cd ('Dados')
x = importdata('matrix.mat'); % matrix adjacency. 
nnodes = size(x,1);
for i=1:nnodes
    degree(i) = nnz(x(i,:)); %nnz(x(i,:)) -> in case of not pondered edges
end
    matrixadj = x;
    Qnew = 0; Qold = 1;
    com = [];
    for i=1:nnodes
        maisligacoes = max(x(i,:));
        quem = find(x(i,:) == maisligacoes);
        aux = [];
        for k = 1:size(quem,2)
            aux(1,k) = quem(1,k);
            aux(2,k) = degree(1,quem(1,k));
        end
        p = find(aux(2,:) == max(aux(2,:)));
        escolha = aux(1,p(1,1));
        if isempty(com)
            com(1,1) = i;
            com(1,2) = escolha;
        else
            if ismember(i,com)
                p = find(com == i);
                linha = mod(p,size(com,1));
                if linha == 0
                    linha = size(com,1);
                end
                if ~ismember(escolha,com)
                   com(linha,size(com,2)+1) = escolha;
                else 
                    pe = find(com == escolha);
                    linhae = mod(pe,size(com,1));
                    if linhae == 0
                        linhae = size(com,1);
                    end 
                    if linhae ~=linha
                        get = [];
                        get = com(linha,[1:end]);
                        cont = 1; oi = size(com,2);
                        while cont <= size(get,2)
                            com(linhae, oi+cont) = get(cont);
                            cont = cont +1;
                        end 
                        com(linha,:) = [];
                    end
                end
            elseif ismember(escolha,com)
                p = find(com == escolha);
                linha = mod(p,size(com,1));
                if linha == 0
                    linha = size(com,1);
                end
                coluna = find(com(linha,:) == 0);
                if ~isempty(coluna)
                    com(linha,coluna(1)) = i;
                else
                    com(linha, size(com,2)+1) = i;
                end
            else
                linha = size(com,1)+1;
                com(linha, 1) = i;
                com(linha, 2) = escolha;
            end
        end
    end
    com( :, all( ~any( com ), 1 ) ) = []; 
    for i=1:nnodes
        for j=1:nnodes
            gg(i,j) = degree(i)*degree(j)/(2*nnodes);
            a = find(com == i);
            b = find(com == j);
            if mod(a,16) == mod(b,16)
                Delta(i,j) = 1;
            else
                Delta(i,j) = 0;
            end
        end
    end
    Q = abs(x-gg).*Delta;
    Q = Q/(2*nnodes);
    Q = max(max(Q));
    Qold = Qnew;
    Qnew = Q;
%while Qnew < Qold
    for i=1:nnodes
        maisligacoes = max(x(i,:));
        quem = find(x(i,:) == maisligacoes);
        aux = [];
        p = find(x(i,:)==max(x(i,:)));
        linhae = mod(find(com == i),size(com,1));
        if linhae ~=0 
            falta = setdiff(p,com(linhae,:));
            if ~isempty(falta)
                escolha = falta(1,1);
                linha = mod(find(com == escolha),size(com,1));
                linhae = mod(find(com == i),size(com,1));
                if linhae ~= linha && linha ~=0
                    get = [];
                    get = com(linha,[1:end]);
                    cont = 1; oi = size(com,2);
                    while cont <= size(get,2)
                        com(linhae, oi+cont) = get(cont);
                        cont = cont +1;
                    end 
                    com(linha,:) = [];
                end
            end
        end
    end
    for i=1:nnodes
        for j=1:nnodes
            gg(i,j) = degree(i)*degree(j)/(2*nnodes);
            a = find(com == i);
            b = find(com == j);
            if mod(a,3) == mod(b,3)
                Delta(i,j) = 1;
            else
                Delta(i,j) = 0;
            end
        end
    end
    Q = abs(x-gg).*Delta;
    Q = Q/(2*nnodes);
    Q = max(max(Q));
    Qold = Qnew;
    Qnew = Q;