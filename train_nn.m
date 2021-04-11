function [mynn]=train_nn(mynn,LR,input,output)
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

        if Y3~=output
            e3=Y3*(1-Y3)*(output(i)-Y3);
            mynn.weights(3,1)=mynn.weights(3,1)+LR*e3*mynn.bias(3);
            for j=1:1:M
                mynn.weights(3,j+1)=mynn.weights(3,j+1)+LR*e3*Y(j);
            end
        else
            e3=0;
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
    
end