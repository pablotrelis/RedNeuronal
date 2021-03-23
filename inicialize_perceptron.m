function [myperceptron] = inicialize_perceptron(n_inputs)
myperceptron.weights=rand(1,n_inputs+1)*2-1;
myperceptron.bias=1;
end

