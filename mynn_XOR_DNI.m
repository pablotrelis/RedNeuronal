function results=mynn_XOR_DNI
%function mynn_XOR_profesor_AND_profesor
% Funcion principal que realiza las funciones de
 %1) Creación de las variables para el banco de entrenamiento y banco de
    %validación
 %2) Creación de la red neuronal de DOS CAPAS
 %3) Entrenamiento de la red neuronal con el set de valores del banco de
    %entrenamiento
 %4) Validación de la red con el banco de validación.
 %5) Calculo y representación del error cometido.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creamos entradas y salidas de entrenamiento para una funcion AND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inicializamos las variables E1, E2 y SE ideales para entrenamiento
VT=5000;                %numero de muestras de entrada
NV=round(VT);           %aseguramos que VT sea un numero entero
E1=round(rand(NV,1));   %vector de NV valores de entrada en pin 1
E2=round(rand(NV,1));   %vector de NV valores de entrada en pin 2
SE=double(xor(E1,E2));   %vector salida ideal para las entradas E1 y E

%%%% Creamos entradas y salidas de validacion E1V, E2V, SEV para test
NVV=20;
E1V=round(rand(NVV,1));
E2V=round(rand(NVV,1));
SEV=double(xor(E1V,E2V));
 
% Inicializamos un perceptron para 2 entradas %%%
mynn=initialize_nn(2);
% Entrenamos el perceptron para un LR, por defecto 0.7
LR=0.7;
mynnT=train_nn(mynn,LR,[E1 E2],SE);
%evaluamos el perceptron
S_est=usenn(mynnT,[E1V E2V]);

%Calculamos y representamos el error cometido
results.error=mean(abs(SEV-S_est));

results.S_est=S_est;
results.SEV=SEV;
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
title('Evaluation of MULTILAYER Ouput for an XOR','FontWeight','bold')
 
    function [mynn]=initialize_nn(n_inputs)
    %function [myperceptron]=initialize_perceptron(n_inputs)
    %funcion que inicializa un perceptron y le asigna pesos aleatorios
    %funcion que inicializa una estructura de perceptron
    % INPUTS:
        % n_inputs: numero de entradas al perceptron
    % OUTPUTS:
        % myperceptron:  estructura con el perceptron
               %myperceptron.bias:       %Bias del perceptron (e.g. 1)
               %myperceptron.weights:    %Pesos del perceptron (habrá tantos
                                         %como indice n_inputs +1 )
    rand('state',sum(100*clock));  %inicializa random en función del reloj
    n_neuronas=3;
    for i=1:1:n_neuronas
        [myperceptron] = inicialize_perceptron(n_inputs);
        mynn.weights(i,:)=myperceptron.weights;
        mynn.bias(i)=myperceptron.bias;
    end
    function [myperceptron] = inicialize_perceptron(n_inputs)
        myperceptron.weights=rand(1,n_inputs+1)*2-1;
        myperceptron.bias=1;
    end
    
    end

    function mynnT=train_nn(mynn,LR,input,output)
    % funcion que modifica los pesos la red para que vaya aprendiendo 
    % ESTE PERCEPTRON UTILIZA:
        % Funcion sigma SIGMOIDAL
        % Entrenamiento BACKPROPAGATION
    % INPUTS:
       % mynn:  estructura con el perceptron
               %mynnT.bias:       %Bias del perceptron (e.g. 1)
               %mynnT.weights:    %Pesos del perceptron 
       % LR: learning rate (e.g. 0.7)
       % input: matriz con valores de entrada de entrenamiento (e.g. [E1 E2])
       % output: vector con valores de salida de entrenamiento (e.g. [SE])
    % OUPUT:
          % mynnT:  estructura con el perceptron ya entrenado
              %mynnT.bias:       %Bias del perceptron (e.g. 1)
              %mynnT.weights:    %Pesos del perceptron ya entrenado  
    [N,M]=size(input);
    for i=1:1:N
        myperceptron1.weights=mynn.weights(1,:);
        myperceptron1.bias=mynn.bias(1);
        myperceptron2.weights=mynn.weights(2,:);
        myperceptron2.bias=mynn.bias(2);
        myperceptron3.weights=mynn.weights(3,:);
        myperceptron3.bias=mynn.bias(3);
        Y(1)=useperceptron(myperceptron1,input(i));
        Y(2)=useperceptron(myperceptron2,input(i));
        Y3=useperceptron(myperceptron2,[Y(1) Y(2)]);
 
        e3=Y3*(1-Y3)*(output(i)-Y3);
        mynn.weights(3,1)=mynn.weights(3,1)+LR*e3*mynn.bias(3);
        for j=1:1:M
            mynn.weights(3,j+1)=mynn.weights(3,j+1)+LR*e3*Y(j);
        end
        e1=Y(1)*(1-Y(1))*mynn.weights(3,2)*e3;
        e2=Y(2)*(1-Y(2))*mynn.weights(3,3)*e3;
        mynn.weights(1,1)=mynn.weights(1,1)+LR*e1*mynn.bias(1);
        for j=1:1:M
            mynn.weights(1,j+1)=mynn.weights(1,j+1)+LR*e1*input(i,j);
        end
        mynn.weights(2,1)=mynn.weights(2,1)+LR*e2*mynn.bias(2);
        for j=1:1:M
            mynn.weights(2,j+1)=mynn.weights(2,j+1)+LR*e2*input(i,j);
        end  
    end
    mynnT=mynn;
    end

    function out=usenn(mynn,input)
    % function out=useperceptron(myperceptron,input)
    % funcion que utiliza el perceptron para calcular las salidas a partir de
    % las entradas de acuerdo con lo que haya aprendido el perceptron en la
    % fase de entrenamiento
    [N,M]=size(input);
    for i=1:1:N
        myperceptron1.weights=mynn.weights(1,:);
        myperceptron1.bias=mynn.bias(1);
        myperceptron2.weights=mynn.weights(2,:);
        myperceptron2.bias=mynn.bias(2);
        myperceptron3.weights=mynn.weights(3,:);
        myperceptron3.bias=mynn.bias(3);
        Y(1)=useperceptron(myperceptron1,input(i));
        Y(2)=useperceptron(myperceptron2,input(i));
        Y3=useperceptron(myperceptron2,[Y(1) Y(2)]);
        out(i)=Y3;
    end
    end

    function [out]=useperceptron(myperceptron,input)
    [N,M]=size(input);
    for i=1:1:N
        for j=1:1:M
            E(j)=input(i,j);        
        end
            v=myperceptron.bias*myperceptron.weights(1);
        for j=1:1:M
            v=v+myperceptron.weights(j+1)*E(j);
        end
        out(i)=1/(1+exp(-v));   
    end
    end

end