function [mynn]=initialize_nn(n_inputs)
    n_neuronas=3;
    for i=1:1:n_neuronas
        [myperceptron] = inicialize_perceptron(n_inputs);
        mynn.weights(i,:)=myperceptron.weights;
        mynn.bias(i)=myperceptron.bias;
    end
    
end