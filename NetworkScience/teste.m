x= [0 0.5 1.125, 3.271, 22.271, 41.813, 59.250, 77.687, 82.542];
y = [0 0 9/96, 39/96, 1, 1, 1, 1, 1];
pensar numa distribui��o considerando o grau m�dio
Verificar se um modelo feito aleat�iro bate com os valores obtidos na realidade
a probabilidade vai ser a densidade de arestas =# de arestas / valor m�ximo que podia existir 
fazer gr�fico da evolu��o do v�rtice de maior grau.
Calcular a vari�ncia da distribui��o de graus
calcular gamma baseado na f�rmula do kmax, deu 2.0046 mais ou menos
plot(x,y,'*')
t = 0:1:83;
hold on
g = 1-exp(-0.15*t);
plot(g)