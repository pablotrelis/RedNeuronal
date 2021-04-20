classdef main_P4_SCBIO_DNI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                  matlab.ui.Figure
        GridLayout                matlab.ui.container.GridLayout
        GridLayout4               matlab.ui.container.GridLayout
        GridLayout13              matlab.ui.container.GridLayout
        ShowInputsButton          matlab.ui.control.Button
        GridLayout9               matlab.ui.container.GridLayout
        matlabperceptronButton    matlab.ui.control.StateButton
        GridLayout8               matlab.ui.container.GridLayout
        mynn_XORButton            matlab.ui.control.StateButton
        GridLayout7               matlab.ui.container.GridLayout
        myperceptron_XORButton    matlab.ui.control.StateButton
        GridLayout5               matlab.ui.container.GridLayout
        myperceptron_ANDButton    matlab.ui.control.StateButton
        GridLayout2               matlab.ui.container.GridLayout
        GridLayout6               matlab.ui.container.GridLayout
        GridLayout12              matlab.ui.container.GridLayout
        RUNButton                 matlab.ui.control.Button
        GridLayout3               matlab.ui.container.GridLayout
        GridLayout11              matlab.ui.container.GridLayout
        num_testEditField         matlab.ui.control.NumericEditField
        NmerotestsLabel           matlab.ui.control.Label
        GridLayout10              matlab.ui.container.GridLayout
        num_entEditField          matlab.ui.control.NumericEditField
        NmeroentrenamientosLabel  matlab.ui.control.Label
        num_test                  matlab.ui.control.Slider
        num_ent                   matlab.ui.control.Slider
        UIAxes                    matlab.ui.control.UIAxes
    end

    
    methods (Access = private) 
        function graf(app,results)
            SEV=results.SEV;
            S_est=results.S_est;
            plot(app.UIAxes,SEV,'ok','LineWidth',2)
            hold(app.UIAxes,'on')
            plot(app.UIAxes,S_est,'xr','LineWidth',2)
            set(app.UIAxes,'FontSize',12) %# Fix font size of the text in the current axes 
            set(app.UIAxes,'FontWeight','bold')  %# Fix Bold text in the current axes 
            xlabel(app.UIAxes,'Number of test','FontWeight','bold')
            ylabel(app.UIAxes,'Output values','FontWeight','bold')
            axis(app.UIAxes,[-1 length(SEV)+1 -0.1 1.3])
            legend(app.UIAxes,'Correct Values','Perceptron Output')
            %title(app.UIAxes,'Evaluation of Perceptron Ouput for an AND','FontWeight','bold')
            hold(app.UIAxes,'off')
        end
        
        function [results,Vinputs]=myperceptron_AND_DNI(~,lg,num_ent,num_test)
        %function myperceptron_AND_DNI
        % Funcion principal que realiza las funciones de
         %1) Creación de las variables para el banco de entrenamiento y banco de
            %validación
         %2) Creación de la red neuronal (perceptron)
         %3) Entrenamiento de la red neuronal con el set de valores del banco de
            %entrenamiento
         %4) Validación de la red con el banco de validación.
         %5) Calculo y representación del error cometido.
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Creamos entradas y salidas de entrenamiento para una funcion AND
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Inicializamos las variables E1, E2 y SE ideales para entrenamiento
        if nargin<1
            lg=0; %lg=0 AND lg=1 XOR
        end
        if nargin<2
            num_ent=5000;
        end
        if nargin<3
            num_test=20;
        end
        % Inicializamos las variables E1, E2 y SE ideales para entrenamiento
        VT=num_ent;                %numero de muestras de entrada
        NV=round(VT);           %aseguramos que VT sea un numero entero
        E1=round(rand(NV,1));   %vector de NV valores de entrada en pin 1
        E2=round(rand(NV,1));   %vector de NV valores de entrada en pin 2
        if lg==1
        SE=double(xor(E1,E2));   %vector salida ideal para las entradas E1 y E
        else
        SE=double(and(E1,E2));
        end
        %%%% Creamos entradas y salidas de validacion E1V, E2V, SEV para test
        NVV=num_test;
        E1V=round(rand(NVV,1));
        E2V=round(rand(NVV,1));
        if lg==1
        SEV=double(xor(E1V,E2V));
        else
        SEV=double(and(E1V,E2V));  
        end
        Vinputs=[E1V E2V];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Creamos entradas y salidas de validacion E1V, E2V, SEV para test
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Inicializamos un perceptron para 2 entradas %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        myperceptron=initialize_perceptron(2);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Entrenamos el perceptron para un LR, por defecto 0.7
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LR=0.7;
        myperceptronT=train_perceptron(myperceptron,LR,[E1 E2],SE);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Evaluamos el perceptron
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S_est=useperceptron(myperceptronT,[E1V E2V]);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Calculamos y representamos el error cometido
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        results.error=mean(abs(SEV-S_est));
        results.SEV=SEV;
        results.S_est=S_est;
        
        % close all
        % figure,
        % plot(SEV,'ok','LineWidth',2),hold on
        % plot(S_est,'xr','LineWidth',2)
        % %set(gcf,'FontWeight','bold')
        % set(gca,'FontSize',12) %# Fix font size of the text in the current axes 
        % set(gca,'FontWeight','bold')  %# Fix Bold text in the current axes 
        % xlabel('Number of test','FontWeight','bold')
        % ylabel('Output values','FontWeight','bold')
        % axis([-1 length(SEV)+1 -0.1 1.3])
        % legend('Correct Values','Perceptron Output')
        % title('Evaluation of Perceptron Ouput for an AND','FontWeight','bold')
         
        %%%%%%%%%%%%%
        % FUNCTIONS %
        %%%%%%%%%%%%%
            function [myperceptron]=initialize_perceptron(n_inputs)
            %function [myperceptron]=initialize_perceptron(n_inputs)
            %funcion que inicializa un perceptron y le asigna pesos aleatorios Esta función crea e inicializa la estructura myperceptron.weights donde se guardan los pesos del perceptron. Los valores de los pesos iniciales deben ser números aleatorios entre -1 y 1.
            % INPUTS:
                % n_inputs: numero de entradas al perceptron
            % OUPUT:
                  % myperceptron:  estructura con el perceptron 
                      % myperceptron.bias:       %Bias del perceptron (e.g. 1)
                      % myperceptron.weights:    %Pesos del perceptron
            myperceptron.weights=rand(1,n_inputs+1)*2-1;
            myperceptron.bias=1;
            end
        
            function myperceptron=train_perceptron(myperceptron,LR,input,output)
            % function myperceptron=train_perceptron(myperceptron,LR,input,output)
            % funcion que modifica los pesos del perceptron para que vaya aprendiendo a
            % a partir de los valores de entrada que se le indican
            % ESTE PERCEPTRON UTILIZA:
                % Funcion sigma SIGMOIDAL
                % Entrenamiento DELTA RULE
            % INPUTS:
               % myperceptron:  estructura con el perceptron
                       %myperceptron.bias:       %Bias del perceptron (e.g. 1)
                       %myperceptron.weights:    %Pesos del perceptron 
               % LR: learning rate (e.g. 0.7)
               % input: matriz con valores de entrada de entrenamiento (e.g. [E1 E2])
               % output: vector con valores de salida de entrenamiento (e.g. [SE])
            % OUPUT:
                  % myperceptron:  estructura con el perceptron ya entrenado
                      %myperceptron.bias:       %Bias del perceptron (e.g. 1)
                      %myperceptron.weights:    %Pesos del perceptron ya entrenado
                a=1;
                [N,M]=size(input);
                for i=1:1:N
                    for j=1:1:M
                        E(j)=input(i,j);        
                    end
                        v=myperceptron.bias*myperceptron.weights(1);
                    for j=1:1:M
                        v=v+myperceptron.weights(j+1)*E(j); %Sumatorio v
                    end
                    % Función SIGMOIDAL
                    y=1/(1+exp(-a*v));
                    
                    if y~=output(i) %Resultado obtenido diferente al esperado
                        e=output(i)-y; %Cálculo de error esperado-obtenido
                            myperceptron.weights(1)=myperceptron.weights(1)+...
                                                        LR*e*myperceptron.bias;
                        for j=1:1:M %Recalculamos los pesos
                            myperceptron.weights(j+1)=myperceptron.weights(j+1)+...
                                                                        LR*e*E(j);
                        end      
                    end     
                end
            end %END function
        
            function out=useperceptron(myperceptron,input)
            % function out=useperceptron(myperceptron,input)
            % funcion que utiliza el perceptron para calcular las salidas a partir de
            % las entradas de acuerdo con lo que haya aprendido el perceptron en la
            % fase de entrenamiento
                [N,M]=size(input);
                for i=1:1:N
                    for j=1:1:M
                        E(j)=input(i,j);        
                    end
                        v=myperceptron.bias*myperceptron.weights(1);
                    for j=1:1:M
                        v=v+myperceptron.weights(j+1)*E(j);
                    end
        %             if v<0
        %                 out(i)=0;
        %             else
        %                 out(i)=1;
        %             end
                    out(i)=1/(1+exp(-v)); 
                end
            end
        
        end %END MAIN function
    
    
        function [results,Vinputs]=mynn_XOR_DNI(~,num_ent,num_test)
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
        Vinputs=[E1V E2V];
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
        
        end %XOR function
       
        
        function results=matlabperceptron_DNI(~,VT,NVV)
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
        %LR=0.7;
        in=[E1 E2];
        net=train(net,in',SE');
        %evaluamos el perceptron
        S_est=net([E1V E2V]');
        
        %Calculamos y representamos el error cometido
        results.error=mean(abs(SEV-S_est));
        results.S_est=S_est;
        results.SEV=SEV;
        
        end %END MAIN function
        
    end %methods
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            global modo
            global inputs
            close all
            modo=0;
            inputs=[];
        end

        % Value changed function: myperceptron_ANDButton
        function myperceptron_ANDButtonValueChanged(app, event)
            value = app.myperceptron_ANDButton.Value;
            global modo
            if value==0
                modo=0;
            else
                app.num_ent.Value=5000;
                app.num_entEditField.Value=5000;
                modo=1;
                app.myperceptron_XORButton.Value=0;
                app.mynn_XORButton.Value=0;
                app.matlabperceptronButton.Value=0;
            end
        end

        % Value changed function: myperceptron_XORButton
        function myperceptron_XORButtonValueChanged(app, event)
            value = app.myperceptron_XORButton.Value;
            global modo
            if value==0
                modo=0;
            else
                app.num_ent.Value=5000;
                app.num_entEditField.Value=5000;
                modo=2;
                app.myperceptron_ANDButton.Value=0;
                app.mynn_XORButton.Value=0;
                app.matlabperceptronButton.Value=0;
            end
        end

        % Value changed function: mynn_XORButton
        function mynn_XORButtonValueChanged(app, event)
            value = app.mynn_XORButton.Value;
            global modo
            if value==0
                modo=0;
            else
                modo=3;
                app.myperceptron_XORButton.Value=0;
                app.myperceptron_ANDButton.Value=0;
                app.matlabperceptronButton.Value=0;
            end
        end

        % Value changed function: matlabperceptronButton
        function matlabperceptronButtonValueChanged(app, event)
            value = app.matlabperceptronButton.Value;
            global modo
            if value==0
                modo=0;
            else
                app.num_ent.Value=1000;
                app.num_entEditField.Value=1000;
                modo=4;
                app.myperceptron_XORButton.Value=0;
                app.myperceptron_ANDButton.Value=0;
                app.mynn_XORButton.Value=0;
            end
        end

        % Button pushed function: RUNButton
        function RUNButtonPushed(app, event)
            global modo
            global inputs
            cla(app.UIAxes,'reset')
            ent=app.num_entEditField.Value;
            test=app.num_testEditField.Value;
            
            if modo==0
            close all   
            elseif modo==1
                [results,inputs]=myperceptron_AND_DNI(app,0,ent,test);
                app.graf(results);
                title(app.UIAxes,'Evaluation of Perceptron Ouput for an AND',...
                    'FontWeight','bold')
            elseif modo==2
                [results,inputs]=myperceptron_AND_DNI(app,1,ent,test);
                app.graf(results);
                title(app.UIAxes,'Evaluation of Perceptron Ouput for an XOR',...
                    'FontWeight','bold')
            elseif modo==3
                [results,inputs]=mynn_XOR_DNI(app,ent,test);
                %uitable('Data',inputs);
                app.graf(results);
                title(app.UIAxes,'Evaluation of Perceptron Ouput for an XOR',...
                    'FontWeight','bold')
            elseif modo==4
                results=matlabperceptron_DNI(app,ent,test);
                app.graf(results);
                title(app.UIAxes,'Evaluation of Perceptron Ouput for an XOR',...
                    'FontWeight','bold')
            end
        end

        % Value changed function: num_testEditField
        function num_testEditFieldValueChanged(app, event)
            value = app.num_testEditField.Value;
            app.num_test.Value=value;
        end

        % Value changed function: num_entEditField
        function num_entEditFieldValueChanged(app, event)
            value = app.num_entEditField.Value;
            if value>10000, value=10000; end
            app.num_ent.Value=value;
        end

        % Value changing function: num_test
        function num_testValueChanging(app, event)
            changingValue = event.Value;
            app.num_testEditField.Value=round(changingValue);
        end

        % Value changing function: num_ent
        function num_entValueChanging(app, event)
            changingValue = event.Value;
            app.num_entEditField.Value=round(changingValue);
        end

        % Button pushed function: ShowInputsButton
        function ShowInputsButtonPushed(app, event)
            global inputs
            close all
            fig_table=uifigure("HandleVisibility",'on','Position',[850,100,400,580]);
            uitable(fig_table,'Data',inputs,'Position',[10,10,380,560]);
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Color = [0.302 0.7451 0.9333];
            app.UIFigure.Position = [100 100 740 580];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'1x', '2x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.BackgroundColor = [0.6588 0.851 0.9294];

            % Create GridLayout2
            app.GridLayout2 = uigridlayout(app.GridLayout);
            app.GridLayout2.ColumnWidth = {'1x'};
            app.GridLayout2.RowHeight = {'2x', '1x'};
            app.GridLayout2.Padding = [0 0 0 0];
            app.GridLayout2.Layout.Row = 1;
            app.GridLayout2.Layout.Column = 2;
            app.GridLayout2.BackgroundColor = [0.4 0.702 0.8];

            % Create UIAxes
            app.UIAxes = uiaxes(app.GridLayout2);
            app.UIAxes.Layout.Row = 1;
            app.UIAxes.Layout.Column = 1;

            % Create GridLayout6
            app.GridLayout6 = uigridlayout(app.GridLayout2);
            app.GridLayout6.ColumnWidth = {'1x'};
            app.GridLayout6.RowHeight = {'3x', '1x'};
            app.GridLayout6.Padding = [0 0 0 0];
            app.GridLayout6.Layout.Row = 2;
            app.GridLayout6.Layout.Column = 1;
            app.GridLayout6.BackgroundColor = [0.4 0.702 0.8];

            % Create GridLayout3
            app.GridLayout3 = uigridlayout(app.GridLayout6);
            app.GridLayout3.RowHeight = {'1x', '3x'};
            app.GridLayout3.ColumnSpacing = 20;
            app.GridLayout3.Padding = [10 0 10 0];
            app.GridLayout3.Layout.Row = 1;
            app.GridLayout3.Layout.Column = 1;
            app.GridLayout3.BackgroundColor = [0.4 0.702 0.8];

            % Create num_ent
            app.num_ent = uislider(app.GridLayout3);
            app.num_ent.Limits = [0 10000];
            app.num_ent.MajorTicks = [0 2500 5000 7500 10000];
            app.num_ent.MajorTickLabels = {'0', '2500', '5000', '7500', '10000'};
            app.num_ent.ValueChangingFcn = createCallbackFcn(app, @num_entValueChanging, true);
            app.num_ent.FontSize = 11;
            app.num_ent.FontColor = [1 1 1];
            app.num_ent.Layout.Row = 2;
            app.num_ent.Layout.Column = 1;
            app.num_ent.Value = 5000;

            % Create num_test
            app.num_test = uislider(app.GridLayout3);
            app.num_test.ValueChangingFcn = createCallbackFcn(app, @num_testValueChanging, true);
            app.num_test.FontSize = 11;
            app.num_test.FontColor = [1 1 1];
            app.num_test.Layout.Row = 2;
            app.num_test.Layout.Column = 2;
            app.num_test.Value = 20;

            % Create GridLayout10
            app.GridLayout10 = uigridlayout(app.GridLayout3);
            app.GridLayout10.ColumnWidth = {'2x', '1x'};
            app.GridLayout10.RowHeight = {'1x'};
            app.GridLayout10.Padding = [0 0 0 0];
            app.GridLayout10.Layout.Row = 1;
            app.GridLayout10.Layout.Column = 1;
            app.GridLayout10.BackgroundColor = [0.4 0.702 0.8];

            % Create NmeroentrenamientosLabel
            app.NmeroentrenamientosLabel = uilabel(app.GridLayout10);
            app.NmeroentrenamientosLabel.BackgroundColor = [0.8118 0.9882 0.9686];
            app.NmeroentrenamientosLabel.Layout.Row = 1;
            app.NmeroentrenamientosLabel.Layout.Column = 1;
            app.NmeroentrenamientosLabel.Text = ' Número entrenamientos';

            % Create num_entEditField
            app.num_entEditField = uieditfield(app.GridLayout10, 'numeric');
            app.num_entEditField.Limits = [0 100000];
            app.num_entEditField.ValueChangedFcn = createCallbackFcn(app, @num_entEditFieldValueChanged, true);
            app.num_entEditField.Layout.Row = 1;
            app.num_entEditField.Layout.Column = 2;
            app.num_entEditField.Value = 5000;

            % Create GridLayout11
            app.GridLayout11 = uigridlayout(app.GridLayout3);
            app.GridLayout11.ColumnWidth = {'2x', '1x'};
            app.GridLayout11.RowHeight = {'1x'};
            app.GridLayout11.Padding = [0 0 0 0];
            app.GridLayout11.Layout.Row = 1;
            app.GridLayout11.Layout.Column = 2;
            app.GridLayout11.BackgroundColor = [0.4 0.702 0.8];

            % Create NmerotestsLabel
            app.NmerotestsLabel = uilabel(app.GridLayout11);
            app.NmerotestsLabel.BackgroundColor = [0.8118 0.9882 0.9686];
            app.NmerotestsLabel.Layout.Row = 1;
            app.NmerotestsLabel.Layout.Column = 1;
            app.NmerotestsLabel.Text = ' Número tests';

            % Create num_testEditField
            app.num_testEditField = uieditfield(app.GridLayout11, 'numeric');
            app.num_testEditField.Limits = [0 100];
            app.num_testEditField.ValueChangedFcn = createCallbackFcn(app, @num_testEditFieldValueChanged, true);
            app.num_testEditField.Layout.Row = 1;
            app.num_testEditField.Layout.Column = 2;
            app.num_testEditField.Value = 20;

            % Create GridLayout12
            app.GridLayout12 = uigridlayout(app.GridLayout6);
            app.GridLayout12.ColumnWidth = {'1x'};
            app.GridLayout12.RowHeight = {'1x'};
            app.GridLayout12.Padding = [10 10 10 0];
            app.GridLayout12.Layout.Row = 2;
            app.GridLayout12.Layout.Column = 1;
            app.GridLayout12.BackgroundColor = [0.4 0.702 0.8];

            % Create RUNButton
            app.RUNButton = uibutton(app.GridLayout12, 'push');
            app.RUNButton.ButtonPushedFcn = createCallbackFcn(app, @RUNButtonPushed, true);
            app.RUNButton.BackgroundColor = [0.149 0.9294 0.6549];
            app.RUNButton.Layout.Row = 1;
            app.RUNButton.Layout.Column = 1;
            app.RUNButton.Text = 'RUN';

            % Create GridLayout4
            app.GridLayout4 = uigridlayout(app.GridLayout);
            app.GridLayout4.ColumnWidth = {'1x'};
            app.GridLayout4.RowHeight = {'2x', '2x', '2x', '2x', '1x'};
            app.GridLayout4.Padding = [0 0 0 0];
            app.GridLayout4.Layout.Row = 1;
            app.GridLayout4.Layout.Column = 1;
            app.GridLayout4.BackgroundColor = [0.4 0.702 0.8];

            % Create GridLayout5
            app.GridLayout5 = uigridlayout(app.GridLayout4);
            app.GridLayout5.ColumnWidth = {'1x'};
            app.GridLayout5.Layout.Row = 1;
            app.GridLayout5.Layout.Column = 1;
            app.GridLayout5.BackgroundColor = [0.4 0.702 0.8];

            % Create myperceptron_ANDButton
            app.myperceptron_ANDButton = uibutton(app.GridLayout5, 'state');
            app.myperceptron_ANDButton.ValueChangedFcn = createCallbackFcn(app, @myperceptron_ANDButtonValueChanged, true);
            app.myperceptron_ANDButton.Text = 'myperceptron_AND';
            app.myperceptron_ANDButton.BackgroundColor = [0.8118 0.9882 0.9725];
            app.myperceptron_ANDButton.Layout.Row = 1;
            app.myperceptron_ANDButton.Layout.Column = 1;

            % Create GridLayout7
            app.GridLayout7 = uigridlayout(app.GridLayout4);
            app.GridLayout7.ColumnWidth = {'1x'};
            app.GridLayout7.Layout.Row = 2;
            app.GridLayout7.Layout.Column = 1;
            app.GridLayout7.BackgroundColor = [0.4 0.702 0.8];

            % Create myperceptron_XORButton
            app.myperceptron_XORButton = uibutton(app.GridLayout7, 'state');
            app.myperceptron_XORButton.ValueChangedFcn = createCallbackFcn(app, @myperceptron_XORButtonValueChanged, true);
            app.myperceptron_XORButton.Text = 'myperceptron_XOR';
            app.myperceptron_XORButton.BackgroundColor = [0.8118 0.9882 0.9686];
            app.myperceptron_XORButton.Layout.Row = 1;
            app.myperceptron_XORButton.Layout.Column = 1;

            % Create GridLayout8
            app.GridLayout8 = uigridlayout(app.GridLayout4);
            app.GridLayout8.ColumnWidth = {'1x'};
            app.GridLayout8.Layout.Row = 3;
            app.GridLayout8.Layout.Column = 1;
            app.GridLayout8.BackgroundColor = [0.4 0.702 0.8];

            % Create mynn_XORButton
            app.mynn_XORButton = uibutton(app.GridLayout8, 'state');
            app.mynn_XORButton.ValueChangedFcn = createCallbackFcn(app, @mynn_XORButtonValueChanged, true);
            app.mynn_XORButton.Text = 'mynn_XOR';
            app.mynn_XORButton.BackgroundColor = [0.8118 0.9882 0.9686];
            app.mynn_XORButton.Layout.Row = 1;
            app.mynn_XORButton.Layout.Column = 1;

            % Create GridLayout9
            app.GridLayout9 = uigridlayout(app.GridLayout4);
            app.GridLayout9.ColumnWidth = {'1x'};
            app.GridLayout9.Layout.Row = 4;
            app.GridLayout9.Layout.Column = 1;
            app.GridLayout9.BackgroundColor = [0.4 0.702 0.8];

            % Create matlabperceptronButton
            app.matlabperceptronButton = uibutton(app.GridLayout9, 'state');
            app.matlabperceptronButton.ValueChangedFcn = createCallbackFcn(app, @matlabperceptronButtonValueChanged, true);
            app.matlabperceptronButton.Text = 'matlabperceptron';
            app.matlabperceptronButton.BackgroundColor = [0.8118 0.9882 0.9686];
            app.matlabperceptronButton.Layout.Row = 1;
            app.matlabperceptronButton.Layout.Column = 1;

            % Create GridLayout13
            app.GridLayout13 = uigridlayout(app.GridLayout4);
            app.GridLayout13.ColumnWidth = {'1x', '2x', '1x'};
            app.GridLayout13.RowHeight = {'1x'};
            app.GridLayout13.Layout.Row = 5;
            app.GridLayout13.Layout.Column = 1;
            app.GridLayout13.BackgroundColor = [0.4 0.702 0.8];

            % Create ShowInputsButton
            app.ShowInputsButton = uibutton(app.GridLayout13, 'push');
            app.ShowInputsButton.ButtonPushedFcn = createCallbackFcn(app, @ShowInputsButtonPushed, true);
            app.ShowInputsButton.BackgroundColor = [0.9804 0.7882 1];
            app.ShowInputsButton.Layout.Row = 1;
            app.ShowInputsButton.Layout.Column = 2;
            app.ShowInputsButton.Text = 'Show Inputs';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = main_P4_SCBIO_DNI

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end