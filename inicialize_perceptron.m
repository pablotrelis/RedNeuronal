function [myperceptron] = inicialize_perceptron(n_inputs)
%function [myperceptron]=initialize_perceptron(n_inputs)
%funcion que inicializa un perceptron y le asigna pesos aleatorios. Esta 
%función crea e inicializa la estructura myperceptron.weights donde se 
%guardan los pesos del perceptron. Los valores de los pesos iniciales 
%deben ser números aleatorios entre -1 y 1.
    % INPUTS:
        % n_inputs: numero de entradas al perceptron
    % OUPUT:
        % myperceptron: estructura con el perceptron ya entrenado
        % myperceptron.bias: %Bias del perceptron
        % myperceptron.weights: %Pesos del perceptron
myperceptron.weights=rand(1,n_inputs+1)*2-1;
myperceptron.bias=1;
end

