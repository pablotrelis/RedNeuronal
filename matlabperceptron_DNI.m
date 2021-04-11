function results=matlabperceptron_DNI(VT,NVV)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Creamos entradas y salidas de entrenamiento para una funcion AND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inicializamos las variables E1, E2 y SE ideales para entrenamiento
if nargin<1
    VT=5000;
end                %numero de muestras de entrada
NV=round(VT);           %aseguramos que VT sea un numero entero
E1=round(rand(NV,1));   %vector de NV valores de entrada en pin 1
E2=round(rand(NV,1));   %vector de NV valores de entrada en pin 2
SE=double(xor(E1,E2));   %vector salida ideal para las entradas E1 y E

%%%% Creamos entradas y salidas de validacion E1V, E2V, SEV para test
if nargin<2
    NVV=20;
end
E1V=round(rand(NVV,1));
E2V=round(rand(NVV,1));
SEV=double(xor(E1V,E2V));
 
% Inicializamos un perceptron para 2 entradas %%%
net=feedforwardnet(3);
% Entrenamos el perceptron para un LR, por defecto 0.7
LR=0.7;
in=[E1 E2];
net=train(net,in',SE');
%evaluamos el perceptron
S_est=net([E1V E2V]');

%Calculamos y representamos el error cometido
results.error=mean(abs(SEV-S_est));
results.S_est=S_est;
results.SEV=SEV;

end %END MAIN function