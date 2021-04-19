function results=mynn_XOR_DNI(num_ent,num_test)
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
if nargin<1
    num_ent=5000;
end
if nargin<2
    num_test=20;
end
% Inicializamos las variables E1, E2 y SE ideales para entrenamiento
VT=num_ent;                %numero de muestras de entrada
NV=round(VT);           %aseguramos que VT sea un numero entero
E1=round(rand(NV,1));   %vector de NV valores de entrada en pin 1
E2=round(rand(NV,1));   %vector de NV valores de entrada en pin 2
SE=double(xor(E1,E2));   %vector salida ideal para las entradas E1 y E

%%%% Creamos entradas y salidas de validacion E1V, E2V, SEV para test
NVV=num_test;
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
        mynn.weights(i,:)=rand(1,n_inputs+1)*2-1;
        mynn.bias(i)=1;
    end  
    end

    function [mynnT]=train_nn(mynn,LR,input,output)
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
        [N,~]=size(input);
            myp1.bias=mynn.bias(1);
            myp2.bias=mynn.bias(2);
            myp3.bias=mynn.bias(3);
        for i=1:1:N
            myp1.weights=mynn.weights(1,:);
            myp2.weights=mynn.weights(2,:);
            myp3.weights=mynn.weights(3,:);

            in=input(i,:);
            out=output(i);

            Y(1)=usep(myp1,in);
            Y(2)=usep(myp2,in);
            Y(3)=usep(myp3,[Y(1) Y(2)]);

            e3=Y(3)*(1-Y(3))*(out-Y(3));
            e1=Y(1)*(1-Y(1))*mynn.weights(3,2)*e3;
            e2=Y(2)*(1-Y(2))*mynn.weights(3,3)*e3;

            mynn.weights(3,1)=mynn.weights(3,1)+LR*e3*mynn.bias(3);
            mynn.weights(3,2)=mynn.weights(3,2)+LR*e3*Y(1);
            mynn.weights(3,3)=mynn.weights(3,3)+LR*e3*Y(2);

            mynn.weights(1,1)=mynn.weights(1,1)+LR*e1*mynn.bias(1);
            mynn.weights(1,2)=mynn.weights(1,2)+LR*e1*in(1);
            mynn.weights(1,3)=mynn.weights(1,3)+LR*e1*in(2);

            mynn.weights(2,1)=mynn.weights(2,1)+LR*e2*mynn.bias(2);
            mynn.weights(2,2)=mynn.weights(2,2)+LR*e2*in(1);
            mynn.weights(2,3)=mynn.weights(2,3)+LR*e2*in(2);

        end
        mynnT=mynn;

    end

    function [out]=usenn(mynn,input)
    [N,~]=size(input);
        for i=1:1:N
            in=input(i,:);
            myp1.weights=mynn.weights(1,:);
            myp1.bias=mynn.bias(1);
            myp2.weights=mynn.weights(2,:);
            myp2.bias=mynn.bias(2);
            myp3.weights=mynn.weights(3,:);
            myp3.bias=mynn.bias(3);

            Y(1)=usep(myp1,in);
            Y(2)=usep(myp2,in);
            Y3=usep(myp3,[Y(1) Y(2)]);
            out(i)=Y3;
        end


    end

    function [out]=usep(myperceptron,input)
        M=length(input);
            for j=1:1:M
                E(j)=input(j);        
            end
                v=myperceptron.bias*myperceptron.weights(1);
            for j=1:1:M
                v=v+myperceptron.weights(j+1)*E(j);
            end
            %SIGMOIDAL
            out=1/(1+exp(-v));   

    end

end