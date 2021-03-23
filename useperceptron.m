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
        if v<0
            out(i)=0;
        else
            out(i)=1;
        end
                
    end

end