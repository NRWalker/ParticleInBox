%Nick Walker
%17 June 2014

classdef ParticleInBoxInterface < handle
    
    properties (SetAccess = private)

        f = NaN;
        F = NaN;
        
    end
    
    methods
        
        function I = ParticleInBoxInterface()
             I.f = 'x .* (1 - x)';
             g = @(x) eval(I.f);
             I.F = ParticleInBox(g);             
        end
        
        function Reset(I)
            clc;
            I.Instructions();
            I.CommandInput();
        end
        
        function Instructions(I)
            disp('======================================================');
            disp('A Particle In A Box');
            disp('======================================================');
            disp('------------------------------------------------------');
            disp('Nicholas Walker');
            disp('------------------------------------------------------');
            disp('0). Input Initial Wave');
            disp(['    f = ' I.f]);
            disp(['    n = ', mat2str(I.F.n())]);
            disp(['    C = ', mat2str(floor(100000 * I.F.C()) / 100000)]);
            disp('1). Plot Real Wave');
            disp('2). Plot Imaginary Wave');
            disp('3). Plot Probability Density Function');
            disp('4). Plot All');
            disp('5). Quit');
        end
        
        function CommandInput(I)
            C = input('[in]: ', 's');
            
            if C == '0'
                disp('[out]: Remember to place "." before operators');
                f = input('[out]: f = ', 's');
                g = @(x) eval(f);
                if abs(g(0)) > 1e-6 || abs(g(1)) > 1e-6
                    disp('[out]: ValueError: Invalid Function Entry');
                    I.CommandInput();
                else
                    I.f = f;
                    I.F.InitialWave(g);
                    disp(['[out]: ' num2str(I.f)]);
                    I.Reset();
                end
            
            elseif C == '1'
                disp('[out]: press any key to close');
                I.F.PlotRealTimeEvolution();
                I.CommandInput();
            
            elseif C == '2'
                disp('[out]: press any key to close');
                I.F.PlotImaginaryTimeEvolution();
                I.CommandInput();
            
            elseif C == '3'
                disp('[out]: press any key to close');
                I.F.PlotDensityTimeEvolution();
                I.CommandInput();
            
            elseif C == '4'
                disp('[out]: press any key to close');
                I.F.PlotTimeEvolution();
                I.CommandInput();
            
            elseif C == '5'
                return;
            
            else
                disp('[out]: ValueError: Invalid Input');
                I.CommandInput();
            end
        end
    end
end            