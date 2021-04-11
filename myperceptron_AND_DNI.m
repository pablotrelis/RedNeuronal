function myperceptron_AND_DNI
%function myperceptron_AND_DNI
% Funcion principal que realiza las funciones de
 %1) Creación de las variables para el banco de entrenamiento y banco de
    %validación
 %2) Creación de la red neuronal (perceptron)
 %3) Entrenamiento de la red neuronal con el set de valores del banco de
    %entrenamiento
 %4) Validación de la red con el banco de validación.
 %5) Calculo y representación del error cometido.
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creamos entradas y salidas de entrenamiento para una funcion AND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inicializamos las variables E1, E2 y SE ideales para entrenamiento
lg=0; %lg=0 AND lg=1 XOR
num_ent=5000;
num_test=20;
posibilidades=[0 1];
for ii=1:1:(num_ent+num_test)
    if rand()<=1/2
        E1(ii)=posibilidades(1);
    else
        E1(ii)=posibilidades(2);
    end
    if rand()<=1/2
        E2(ii)=posibilidades(1);
    else
        E2(ii)=posibilidades(2);
    end
    if lg==0
        SE(ii)=and(E1(ii),E2(ii));
    else
        SE(ii)=xor(E1(ii),E2(ii));
    end
end
E1V=E1(num_ent+1:length(E1));
E2V=E2(num_ent+1:length(E2));
SEV=SE(num_ent+1:length(SE));
E1=E1(1:num_ent);
E2=E2(1:num_ent);
SE=SE(1:num_ent);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Creamos entradas y salidas de validacion E1V, E2V, SEV para test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inicializamos un perceptron para 2 entradas %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myperceptron=initialize_perceptron(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Entrenamos el perceptron para un LR, por defecto 0.7
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LR=0.7;
myperceptronT=train_perceptron(myperceptron,LR,[E1;E2]',SE);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluamos el perceptron
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_est=useperceptron(myperceptronT,[E1V;E2V]');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculamos y representamos el error cometido
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error=mean(abs(SEV-S_est));
 
close all
figure,
plot(SEV,'ok','LineWidth',2),hold on
plot(S_est,'xr','LineWidth',2)
%set(gcf,'FontWeight','bold')
set(gca,'FontSize',12) %# Fix font size of the text in the current axes 
set(gca,'FontWeight','bold')  %# Fix Bold text in the current axes 
xlabel('Number of test','FontWeight','bold')
ylabel('Output values','FontWeight','bold')
axis([-1 length(SEV)+1 -0.1 1.3])
legend('Correct Values','Perceptron Output')
title('Evaluation of Perceptron Ouput for an AND','FontWeight','bold')
 
%%%%%%%%%%%%%
% FUNCTIONS %
%%%%%%%%%%%%%
    function [myperceptron]=initialize_perceptron(n_inputs)
    %function [myperceptron]=initialize_perceptron(n_inputs)
    %funcion que inicializa un perceptron y le asigna pesos aleatorios Esta función crea e inicializa la estructura myperceptron.weights donde se guardan los pesos del perceptron. Los valores de los pesos iniciales deben ser números aleatorios entre -1 y 1.
    % INPUTS:
        % n_inputs: numero de entradas al perceptron
    % OUPUT:
          % myperceptron:  estructura con el perceptron 
              % myperceptron.bias:       %Bias del perceptron (e.g. 1)
              % myperceptron.weights:    %Pesos del perceptron
    myperceptron.weights=rand(1,n_inputs+1)*2-1;
    myperceptron.bias=1;
    end

    function myperceptron=train_perceptron(myperceptron,LR,input,output)
    % function myperceptron=train_perceptron(myperceptron,LR,input,output)
    % funcion que modifica los pesos del perceptron para que vaya aprendiendo a
    % a partir de los valores de entrada que se le indican
    % ESTE PERCEPTRON UTILIZA:
        % Funcion sigma SIGMOIDAL
        % Entrenamiento DELTA RULE
    % INPUTS:
       % myperceptron:  estructura con el perceptron
               %myperceptron.bias:       %Bias del perceptron (e.g. 1)
               %myperceptron.weights:    %Pesos del perceptron 
       % LR: learning rate (e.g. 0.7)
       % input: matriz con valores de entrada de entrenamiento (e.g. [E1 E2])
       % output: vector con valores de salida de entrenamiento (e.g. [SE])
    % OUPUT:
          % myperceptron:  estructura con el perceptron ya entrenado
              %myperceptron.bias:       %Bias del perceptron (e.g. 1)
              %myperceptron.weights:    %Pesos del perceptron ya entrenado
        a=1;
        [N,M]=size(input);
        for i=1:1:N
            for j=1:1:M
                E(j)=input(i,j);        
            end
                v=myperceptron.bias*myperceptron.weights(1);
            for j=1:1:M
                v=v+myperceptron.weights(j+1)*E(j); %Sumatorio v
            end
            % Función SIGMOIDAL
            y=1/(1+exp(-a*v));
            
            if y~=output(i) %Resultado obtenido diferente al esperado
                e=output(i)-y; %Cálculo de error esperado-obtenido
                    myperceptron.weights(1)=myperceptron.weights(1)+...
                                                LR*e*myperceptron.bias;
                for j=1:1:M %Recalculamos los pesos
                    myperceptron.weights(j+1)=myperceptron.weights(j+1)+...
                                                                LR*e*E(j);
                end      
            end     
        end
    end %END function

    function out=useperceptron(myperceptron,input)
    % function out=useperceptron(myperceptron,input)
    % funcion que utiliza el perceptron para calcular las salidas a partir de
    % las entradas de acuerdo con lo que haya aprendido el perceptron en la
    % fase de entrenamiento
        [N,M]=size(input);
        for i=1:1:N
            for j=1:1:M
                E(j)=input(i,j);        
            end
                v=myperceptron.bias*myperceptron.weights(1);
            for j=1:1:M
                v=v+myperceptron.weights(j+1)*E(j);
            end
            if v<0
                out(i)=0;
            else
                out(i)=1;
            end

        end
    end

end %END MAIN function