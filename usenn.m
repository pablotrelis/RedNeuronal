function [out]=usenn(mynn,input)
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