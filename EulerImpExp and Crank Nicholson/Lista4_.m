%%% Código desenvolvido por Raquel Mattoso para a disciplina de Método de Elementos Finitos %%%
%%% Inserindo o domínio %%%
clear; clc;
a = 0;
l = 1; % comprimento da barra 
h = 0.05; % tamanho da discretização no espaço
X = a:h:l; % conjunto de pontos de X
Xfront = h:h:(l-h); % desconsiderando pontos da fronteira
n = size(X,2);
tmax = 1.0;
t = 0;
Deltat = 0.0004; % para Euler Explícito Deltat <= 1/4*h^2
cont = 1;
flag = 3; % 1 -> Implícito, 2 -> Explícito, 3-> CN, 4 -> Pós processamento
f = 0; % fonte do problema
k = 1; % capacidade térmica
u0t = 0; ult = 0; % condição de contorno dirichlet homogêneo
ux0 = 1/(pi*pi)*sin(pi*Xfront); ux0 = ux0'; % condição inicial

syms x % tratado como símbolo
phi1(x) = x/h; 
phi2(x) = (h-x)/h;
dphi1(x) = diff(phi1(x));
dphi2(x) = diff(phi2(x));
M = zeros(n-2,n-2); K = zeros(n-2,n-2);
for j=1:(n-2)
    for i=1:(n-2)
        if i==j
            M(i,j) = int(phi1*phi1,0,h) + int(phi2*phi2,0,h);
            K(i,j) = k*(int(dphi1*dphi1,0,h) + int(dphi2*dphi2,0,h));
        end
        if i==j+1 || i==j-1
            M(i,j) = int(phi1*phi2,0,h);
            K(i,j) = k*int(dphi1*dphi2,0,h);
        end
    end
end
% Euler implícito

if flag == 1
    tic;
    aux = Deltat*(int(phi1*f,0,h)+ int(phi2*f,0,h));
    while t <= tmax
        B = M*ux0 + aux;
        B(1,1) = B(1,1) - u0t; B(n-2,1) = B(n-2,1) - ult;
        ux = (Deltat*K+M)\B;
        ux0 = eval(ux);
        uh = [0,ux0',0];
        EI(cont,:) = uh(1,:);
        cont = cont + 1;
        t = t+Deltat;
    end
    save 'EI.mat' EI;
    toc;
    time = toc;
    save timeEI.mat time;    
end
% Euler Explícito 
if flag == 2
    tic;
    aux = Deltat*(int(phi1*f,0,h)+int(phi2*f,0,h));
    while t <= tmax
        ux = (eye(n-2,n-2) - Deltat*inv(M)*K)*ux0 + aux;
        ux(1,1) = ux(1,1) - u0t; ux(n-2,1) = ux(n-2,1) - ult;
        ux0 = eval(ux);
        uh = [0,ux0',0];
        EE(cont,:) = uh(1,:);
        cont = cont +1;
        t = t+Deltat;
    end
    save EE.mat EE;
    toc;
    time = toc;
    save timeEE.mat time; 
end
% Crank - Nicolson
if flag == 3
    tic;
    aux = Deltat*(int(phi1*f,0,h)+int(phi2*f,0,h));
    while t <= tmax
        B = (M - 1/2*Deltat*K)*ux0 + aux;
        B(1,1) = B(1,1) - u0t; B(n-2,1) = B(n-2,1) - ult;
        ux = (M+(Deltat/2)*K)\B;
        ux0 = eval(ux);
        uh = [0,ux0',0];
        CN(cont,:) = uh(1,:);
        cont = cont + 1;
        t = t+Deltat;
    end
    save CN.mat CN;
    toc;
    time = toc;
    save timeCN.mat time; 
end
%% 
% Pós processamento

if flag == 4
     EI = importdata('EI.mat');
     EE = importdata('EE.mat');
     CN = importdata('CN.mat');
     aux = tmax/Deltat;
    for i=1:aux
        Exata(i,:) = (1.0/(pi*pi))*sin(pi*X)*exp(-(k*pi*pi)*i*Deltat);
    end 
    
    figure(1);
    plot(X,Exata(10,:),'k'); hold on;
    plot(X,EI(10,:),'b'); hold on;
    plot(X,EE(10,:),'r'); hold on;
    plot(X,CN(10,:),'g'); legend('Exata','Euler Implícito','Euler Explícito','Crank-Nicolson'); 
    title('Comparação entre as soluções para t=0.1');ylabel('$u_h , u$','interpreter','latex'); 
    xlabel('$x$(u.m.)','interpreter','latex');
    
    figure(2);
    plot(X,Exata(aux,:),'k'); hold on;
    plot(X,EI(aux,:),'b'); hold on;
    plot(X,EE(aux,:),'r'); hold on;
    plot(X,CN(aux,:),'g'); legend('Exata','Euler Implícito','Euler Explícito','Crank-Nicolson'); 
    title('Comparação entre as soluções para t=1.0'); ylabel('$u_h , u$','interpreter','latex'); 
    xlabel('$x$(u.m.)','interpreter','latex');
        
    for i=1:(tmax/Deltat)
        ErroEI(i) = norm(Exata(i,:)-EI(i,:));
        ErroEE(i) = norm(Exata(i,:)-EE(i,:));
        ErroCN(i) = norm(Exata(i,:)-CN(i,:));
    end
    figure(3);
    tempo = 0:Deltat:tmax;
    ErroEI = [0,ErroEI];
    ErroEE = [0,ErroEE];
    ErroCN = [0,ErroCN];
    plot(tempo,ErroEI,'b'); hold on;
    plot(tempo,ErroEE,'r'); hold on;
    plot(tempo,ErroCN,'g'); legend('Euler Implícito','Clark-Nicolson');
    title('Erro Euler Explícito com condição violada'); 
    xlabel('$t$(u.t.)','interpreter','latex'); ylabel('$\|u - u_h\|$','interpreter','latex');
end