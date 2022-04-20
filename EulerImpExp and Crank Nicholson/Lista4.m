%%% Código desenvolvido por Raquel Mattoso para a disciplina de Método de Elementos Finitos %%%
% O problema parabólico escolhido na forma diferencial é o problema de calor
% $\frac{\partial u}{\partial t} - k \frac{\partial ^2 u}{\partial x^2} = f 
% \qquad 0<x<l, \quad t>0$,
% 
% $$u(0,t) = u(l,t) = 0 \qquad t>0,$$
% 
% $u(x,0) = u_0(x) \qquad 0<x<l$.

%%% Inserindo o domínio %%%
clear; clc;
tic;
a = 0;
l = 1; % comprimento da barra 
h = 0.1; % tamanho da discretização no espaço
X = a:h:l; % conjunto de pontos de X
Xfront = h:h:(l-h); % desconsiderando pontos da fronteira
n = size(X,2);
tmax = 1.0;
t = 0;
Deltat = 0.2; % para Euler Explícito Deltat <= 1/4*h^2
cont = 1;
flag = 2; % 1 -> Implícito, 2 -> Explícito, 3-> CN, 4 -> Pós processamento
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
% EULER IMPLÍCITO

if flag == 1
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
    aux = Deltat*(int(phi1*f,0,h)+int(phi2*f,0,h));
    while t <= tmax
        B = M*ux0 - Deltat*K*ux0;
        B(1,1) = B(1,1) - u0t; B(n-2,1) = B(n-2,1) - ult;
        ux = M\B; 
        ux0 = (ux);
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
     %EI = importdata('EI.mat');
     EE = importdata('EE.mat');
     %CN = importdata('CN.mat');
     aux = tmax/Deltat;
    for i=1:aux
        Exata(i,:) = (1.0/(pi*pi))*sin(pi*X)*exp(-(k*pi*pi)*i*Deltat);
    end 
    
%     figure(1);
%     plot(X,Exata(1,:),'k'); hold on;
%     plot(X,EI(1,:),'b'); hold on;
%     plot(X,EE(1,:),'r'); hold on;
%     plot(X,CN(1,:),'g'); legend('Exata','Euler Implícito','Euler Explícito','Crank-Nicolson'); 
%     title('Comparação entre as soluções para t=0.0');ylabel('$u_h$','interpreter','latex'); 
%     xlabel('$x$(u.m.)','interpreter','latex');
%     
%     figure(2);
%     plot(X,Exata(aux,:),'k'); hold on;
%     plot(X,EI(aux,:),'b'); hold on;
%     plot(X,EE(aux,:),'r'); hold on;
%     plot(X,CN(aux,:),'g'); legend('Exata','Euler Implícito','Euler Explícito','Crank-Nicolson'); 
%     title('Comparação entre as soluções para t=1.0'); ylabel('$u_h$','interpreter','latex'); 
%     xlabel('$x$(u.m.)','interpreter','latex');
        
    for i=1:(tmax/Deltat)
        %ErroEI(i) = norm(Exata(i,:)-EI(i,:));
        ErroEE(i) = norm(Exata(i,:)-EE(i,:));
        %ErroCN(i) = norm(Exata(i,:)-CN(i,:));
    end
    figure(3);
    tempo = Deltat:Deltat:tmax;
    %plot(tempo,ErroEI,'b'); hold on;
    plot(tempo,ErroEE,'r'); hold on;
    %plot(tempo,ErroCN,'g'); legend('Euler Implícito','Euler Explícito','Clark-Nicolson');
    title('Erro Euler Explícito'); 
    xlabel('$t$(u.t.)','interpreter','latex'); ylabel('$\|Exata - Numerica\|$','interpreter','latex');
end