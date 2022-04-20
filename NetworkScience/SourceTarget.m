clear; clc; close all;
%% Builds the source-target relation
cd ('Dados')
T = readtable('LN.csv');
Detective = []; Persons = []; Murder=[]; CEP =[]; Place = []; Motive = []; Words = [];
line = size(T,1);
col = size(T,2);
%% Detective
aux=0;
for i=1:line
    id = table2array(T(i,4));
    equal = find(table2array(T(:,4))==id);
    n = size(equal,1);
    for j=1:n
         if i~=equal(j)
            Detective(j+aux,1) = i;
            Detective(j+aux,2) = equal(j);
         else
             j=j+1;
         end
    end
    aux = size(Detective,1);
end
%% Personagens
aux=0; listpers = {};
for i=1:line
    auxchar = join(erase(string(T{i,5}), "'"), '', 2);
    C = strsplit(auxchar,{' ',','},'CollapseDelimiters',true);
    for k=1:size(C,2)
        listpers{i}(1,k) = str2double(C(k));
    end
end
for i=1:line
    for m=i+1:line
        for k=1:size(listpers{i},2)
            id = listpers{i}(1,k);
            for n=1:size(listpers{m},2)
                if listpers{m}(1,n) == id
                        Persons(n+aux,1) = i;
                        Persons(n+aux,2) = m;
                end
                aux = size(Persons,1);
            end
        end
    end
end
%% Tipo de morte
aux = 0; listmurder = {};
for i=1:line
    auxchar = join(erase(string(T{i,6}), "'"), '', 2);
    listmurder{i} = strsplit(auxchar,{' ',','},'CollapseDelimiters',true);
end

for i=1:line
    for m=i+1:line
        for k=1:size(listmurder{i},2)
            id = listmurder{i}(1,k);
            for n=1:size(listmurder{m},2)
                if listmurder{m}(1,n) == id
                        Murder(n+aux,1) = i;
                        Murder(n+aux,2) = m;
                end
                aux = size(Murder,1);
            end
        end
    end
end
%% CEP
aux=0; listcep = {};
for i=1:line
    auxchar = join(erase(string(T{i,7}), "'"), '', 2);
    listcep{i} = strsplit(auxchar,{' ',','},'CollapseDelimiters',true);
end
for i=1:line
    for m=i+1:line
        for k=1:size(listcep{i},2)
            id = listcep{i}(1,k);
            for n=1:size(listcep{m},2)
                if listcep{m}(1,n) == id
                        CEP(n+aux,1) = i;
                        CEP(n+aux,2) = m;
                end
                aux = size(CEP,1);
            end
        end
    end
end
%% lugar
aux=0; listplace = {};
for i=1:line
    auxchar = join(erase(string(T{i,8}), "'"), '', 2);
    listplace{i} = strsplit(auxchar,{' ',','},'CollapseDelimiters',true);
end
for i=1:line
    for m=i+1:line
        for k=1:size(listplace{i},2)
            id = listplace{i}(1,k);
            for n=1:size(listplace{m},2)
                if listplace{m}(1,n) == id
                        Place(n+aux,1) = i;
                        Place(n+aux,2) = m;
                end
                aux = size(Place,1);
            end
        end
    end
end

%% motivo
aux=0; listmotive = {};
for i=1:line
    auxchar = join(erase(string(T{i,9}), "'"), '', 2);
    listmotive{i} = strsplit(auxchar,{' ',','},'CollapseDelimiters',true);
end
for i=1:line
    for m=i+1:line
        for k=1:size(listmotive{i},2)
            id = listmotive{i}(1,k);
            for n=1:size(listmotive{m},2)
                if listmotive{m}(1,n) == id
                        Motive(n+aux,1) = i;
                        Motive(n+aux,2) = m;
                end
                aux = size(Motive,1);
            end
        end
    end
end

%% words in common
aux = 0; listwords = {}; listaux ={}; j=1;
for i=1:line
    auxchar = join(erase(string(T{i,3}),"'"),'',2);
    listaux{i} = strsplit(auxchar,{' ',' '},'CollapseDelimiters',true);
    for cont =1:size(listaux{i},2)
        if  size(char(listaux{i}(1,cont)),2) > 2
            if lower(listaux{i}(1,cont)) ~= char('the')
                if lower(listaux{i}(1,cont)) ~= char('with')
                    listwords{i}(1,j) = listaux{i}(1,cont);
                    j = j+1;
                end
            end
        end
    end
    j =1;
end
for i=1:line
    for m=i+1:line
        for k=1:size(listwords{i},2)
            id = listwords{i}(1,k);
            for n=1:size(listwords{m},2)
                if listwords{m}(1,n) == id
                        Words(n+aux,1) = i;
                        Words(n+aux,2) = m;
                end
                aux = size(Words,1);
            end
        end
    end
end
%% Pós processamento 
Detective = Detective(any(Detective,2),:);
Persons = Persons(any(Persons,2),:);
Murder = Murder(any(Murder,2),:);
CEP = CEP(any(CEP,2),:);
Place = Place(any(Place,2),:);
Motive = Motive(any(Motive,2),:);
Words = Words(any(Words,2),:);
csvwrite('Detective.csv',Detective);
csvwrite('Persons.csv',Persons);
csvwrite('Murder.csv',Murder);
csvwrite('CEP.csv',CEP);
csvwrite('Place.csv',Place);
csvwrite('Motive.csv',Motive);
csvwrite('Words.csv',Words);
cd ..