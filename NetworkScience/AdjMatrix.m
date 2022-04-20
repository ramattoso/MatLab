clear; clc; close all;
n=96; AdjM = cell(n,1); 
%% Builds the Adjacency Matrix
cd ('Dados')
T = readtable('Detective.csv'); 
line = size(T,1); MA_D = zeros(n,n);
for i=1:line
    source = table2array(T(i,1));
    target = table2array(T(i,2));
    MA_D(source,target) = MA_D(source,target) + 1;
    MA_D(target,source) = MA_D(source,target);
end
MA_D = MA_D/2;
T = readtable('Persons.csv');
line = size(T,1); MA_Pe = zeros(n,n);
for i=1:line
    source = table2array(T(i,1));
    target = table2array(T(i,2));
    MA_Pe(source,target) = MA_Pe(source,target)+1;
    MA_Pe(target,source) = MA_Pe(source, target);
end

T = readtable('CEP.csv');
line = size(T,1); MA_CEP = zeros(n,n);
for i=1:line
    source = table2array(T(i,1));
    target = table2array(T(i,2));
    MA_CEP(source,target) = MA_CEP(source,target)+1;
    MA_CEP(target,source) = MA_CEP(source, target);
end

T = readtable('Place.csv');
line = size(T,1); MA_P = zeros(n,n);
for i=1:line
    source = table2array(T(i,1));
    target = table2array(T(i,2));
    MA_P(source,target) = MA_P(source,target)+1;
    MA_P(target,source) = MA_P(source, target);
end

T = readtable('Murder.csv');
line = size(T,1); MA_M = zeros(n,n);
for i=1:line
    source = table2array(T(i,1));
    target = table2array(T(i,2));
    MA_M(source,target) = MA_M(source,target)+1;
    MA_M(target,source) = MA_M(source, target);
end

T = readtable('Motive.csv');
line = size(T,1); MA_Mo = zeros(n,n);
for i=1:line
    source = table2array(T(i,1));
    target = table2array(T(i,2));
    MA_Mo(source,target) = MA_Mo(source,target)+1;
    MA_Mo(target,source) = MA_Mo(source, target);
end

T = readtable('Words.csv');
line = size(T,1); MA_W = zeros(n,n);
for i=1:line
    source = table2array(T(i,1));
    target = table2array(T(i,2));
    MA_W(source,target) = MA_W(source,target)+1;
    MA_W(target,source) = MA_W(source, target);
end
Matrix = MA_Pe + MA_W + MA_Mo + MA_P + MA_M + MA_CEP + MA_D;
% for i=1:size(Matrix,1)
%      for j=1:size(Matrix,2)
%          if Matrix(i,j) ~=0
%              Matrix(i,j) = 1;
%          end
%      end
% end
% l = 1; k=1;
for i=1:size(Matrix,1)
     for j=1:size(Matrix,2)
         if Matrix(i,j) < 5
             AuxMA(i,j) = 0;
         else
             AuxMA(i,j) = Matrix(i,j);
         end
     end
end

save matrix.mat Matrix
cd ..
%% printing books names on AdjM
% T = readtable('LN.csv'); 
% for i=1:n
%     AdjM{i,1} = join(erase(string(T{i,3}), "'"), "", 2);
% end
% %writecell(AdjM, 'AdjM.txt');
% save AdjM.mat AdjM;