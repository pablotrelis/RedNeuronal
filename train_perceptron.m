function [myperceptron]=train_perceptron(myperceptron,LR,input,output)
    
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
            y=0;
        else
            y=1;
        end
        
        if y~=output(i)
            e=output(i)-y;
                myperceptron.weights(1)=myperceptron.weights(1)+LR*e*myperceptron.bias;
            for j=1:1:M
                myperceptron.weights(j+1)=myperceptron.weights(j+1)+LR*e*E(j);
            end
        
        end
        
    end



end