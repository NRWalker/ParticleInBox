%Nick Walker
%10 June 2014

classdef ParticleInBox < handle
    
    properties (SetAccess = private)
        
        dx = 1e-4;
        x = 0:1e-4:1;
        S = @(n) sin((0:1e-4:1) .* n .* pi);
        dt = NaN;
        f = NaN;
        X = NaN;
        RT = NaN;
        IT = NaN;
        E = NaN;
        c = NaN;
        C = NaN;
        n = NaN;
        e = 1e-6;

    end
    
    methods
        
        function F = ParticleInBox(f) 
            F.f = Normalize(F.x, f(F.x));
            if abs(f(0)) > F.e || abs(f(1)) > F.e
                disp('Invalid Function Entry, Default Value Invoked');
                f = @(x) x .* (1 - x);
                F.f = Normalize(F.x, f(F.x));
            end
            F.X = @(n) sqrt(2) * F.S(n);
            F.RT = @(n, t) cos(F.E(n) .* t);
            F.IT = @(n, t) -sin(F.E(n) .* t);
            F.E = @(n) (n .* pi) ^ 2 ./ (2);        
            F.c = @(n) trapz(F.x, F.f .* F.X(n));
            NC = FourierConstants(F.c, F.e);
            F.C = NC(2, :);
            F.n = NC(1, :);
            if length(F.n) == 1
                w = F.n(1);
            else
                w = sum(sqrt(F.n) .* F.C);
            end
            t = 4 / (w ^ 2 * pi);
            F.dt = t * 1e-2; 
        end
        
        function InitialWave(F, f) 
            if abs(f(0)) > F.e || abs(f(1)) > F.e
                disp('Invalid Function Entry, Prior Value Retained');
                return
            end
            F.f = Normalize(F.x, f(F.x));
            F.c = @(n) trapz(F.x, F.f .* F.X(n));
            NC = FourierConstants(F.c, F.e);
            F.C = NC(2, :);
            F.n = NC(1, :);
            if length(F.n) == 1
                w = F.n(1);
            else
                w = sum(sqrt(F.n) .* F.C);
            end
            t = 4 / (w ^ 2 * pi);
            F.dt = t * 1e-2;
        end
        
        function [X] = SpatialArray(F, n)
            X = F.X(n);
        end
        
        function [RT] = RealTimeAmplitude(F, n, t)
            RT = F.RT(n, t); 
        end
        
        function [IT] = ImaginaryTimeAmplitude(F, n, t)
            IT = F.IT(n, t);
        end
        
        function c = FourierConstant(F, n)
            c = F.c(n);
        end
        
        function [C] = FourierConstantArray(F)
            C = F.C;
        end
        
        function [n] = WaveNumberArray(F)
            n = F.n;
        end
        
        function E = WaveEnergy(F, n)
            E = F.E(n);
        end
        
        function [RF] = RealWaveArray(F, n, t)
            RF = F.c(n) * F.X(n) .* F.RT(n, t);
        end
           
        function [IF] = ImaginaryWaveArray(F, n, t)
            IF = F.c(n) * F.X(n) .* F.IT(n, t);
        end
        
        function [P] = DensityArray(F, n, t)
            R = F.RealWaveArray(n, t);
            I = F.ImaginaryWaveArray(n, t);
            P = R .^ 2 + I .^ 2;
        end
        
        function [RF] = RealSuperPositionArray(F, t)
            A = zeros(length(F.n), length(F.x));
            for m = 1:length(F.n)
                A(m, :) = F.RealWaveArray(F.n(m), t);
            end
            RF = sum(A, 1);
        end
        
        function [IF] = ImaginarySuperPositionArray(F, t)
            A = zeros(length(F.n), length(F.x));
            for m = 1:length(F.n)
                A(m, :) = F.ImaginaryWaveArray(F.n(m), t);
            end
            IF = sum(A, 1);
        end
        
        function [P] = DensitySuperPositionArray(F, t)
            R = F.RealSuperPositionArray(t);
            I = F.ImaginarySuperPositionArray(t);
            P = R .^ 2 + I .^ 2;
        end
        
        function E = ExpectationEnergy(F)
            A = zeros(1, length(F.n));
            for m = 1:length(F.n)
                A(m) = F.c(m) .^ 2 .* F.E(m);
            end
            E = sum(A);
        end
        
        function x = ExpectationPosition(F, t)
            x = trapz(F.x, F.x .* F.DensitySuperPositionArray(t));
        end
        
        function p = ExpectationMomentum(F, t)
            xr = F.ExpectationPosition(t + F.e / 2);
            xl = F.ExpectationPosition(t - F.e / 2);
            p = (xr - xl) / F.e;
        end
        
        function x = ExpectationSquarePosition(F, t)
            x = trapz(F.x, F.x .^ 2 .* F.DensitySuperPositionArray(t));
        end
        
        function p = ExpectationSquareMomentum(F)
            p = 2 * F.ExpectationEnergy();
        end
        
        function s = StandardDeviationPosition(F, t)
            s = sqrt(F.ExpectationSquarePosition(t) - F.ExpectationPosition(t) .^ 2);
        end
        
        function s = StandardDeviationMomentum(F, t)
            s = sqrt(F.ExpectationSquareMomentum() - F.ExpectationMomentum(t) .^ 2);
        end
        
        function H = HeisenbergUncertainty(F, t)
            H = F.StandardDeviationPosition(t) * F.StandardDeviationMomentum(t);
        end
        
        function PlotInitialWave(F)
            x = F.x;
            y = F.f(t);
            figure;
            plot(x, y);
            title('f(x)');
            ylabel('f(x)');
            xlabel('x');
            y1 = 1.2 * min(y);
            y2 = 1.2 * max(y);
            axis([0 1 y1 y2]);
        end
        
        function PlotRealWave(F, t)
            x = F.x;
            y = F.RealSuperPositionArray(t);
            figure;
            plot(x, y);
            title('\Re(\Psi)');
            ylabel('\Re(\Psi)');
            xlabel('x');
            y1 = 1.2 * min(y);
            y2 = 1.2 * max(y);
            axis([0 1 y1 y2]);
        end
        
        function PlotImaginaryWave(F, t)
            x = F.x;
            y = F.ImaginarySuperPositionArray(t);
            figure;
            plot(x, y);
            title('\Im(\Psi)');
            ylabel('\Im(\Psi)');
            xlabel('x');
            y1 = 1.2 * min(y);
            y2 = 1.2 * max(y);
            axis([0 1 y1 y2]);
        end
        
        function PlotDensity(F, t)
            x = F.x;
            y = F.DensitySuperPositionArray(t);
            figure;
            plot(x, y);
            title('\Psi^* \Psi');
            ylabel('\Psi^* \Psi');
            xlabel('x');
            y1 = 1.2 * min(y);
            y2 = 1.2 * max(y);
            axis([0 1 y1 y2]);
        end
        
        function PlotRealTimeEvolution(F)
            dt = F.dt;
            x = F.x;
            y = F.DensitySuperPositionArray(0);
            h = max(y);
            global Key_Press;
            Key_Press = 0;
            set(gcf, 'KeyPressFcn', @GetKeyPress);
            n = 0;
            while ~Key_Press
                y = F.RealSuperPositionArray(n * dt);
                plot(x, y);
                title('\Re(\Psi)');
                ylabel('\Re(\Psi)');
                xlabel('x');
                axis([0 1 -h h]);
                n = n + 1;
                drawnow;
            end
            close
        end
        
        function PlotImaginaryTimeEvolution(F)
            dt = F.dt;
            x = F.x;
            y = F.DensitySuperPositionArray(0);
            h = max(y);
            global Key_Press;
            Key_Press = 0;
            set(gcf, 'KeyPressFcn', @GetKeyPress);
            n = 0;
            while ~Key_Press
                y = F.ImaginarySuperPositionArray(n * dt);
                plot(x, y);
                title('\Im(\Psi)');
                ylabel('\Im(\Psi)');
                xlabel('x');
                axis([0 1 -h h]);
                n = n + 1;
                drawnow;
            end
            close
        end
        
        function PlotDensityTimeEvolution(F)
            dt = F.dt;
            x = F.x;
            y = F.DensitySuperPositionArray(0);
            h = 1.2 * max(y);
            global Key_Press;
            Key_Press = 0;
            set(gcf, 'KeyPressFcn', @GetKeyPress);
            n = 0;
            while ~Key_Press
                y = F.DensitySuperPositionArray(n * dt);
                plot(x, y);
                title('\Psi^* \Psi');
                ylabel('\Psi^* \Psi');
                xlabel('x');
                axis([0 1 0 h]);
                n = n + 1;
                drawnow;
            end
            close(gcf);
        end
        
        function PlotTimeEvolution(F)
            dt = F.dt;
            x = F.x;
            y = F.DensitySuperPositionArray(0);
            h = max(y);
            suptitle('Particle Time Evolution');
            global Key_Press;
            Key_Press = 0;
            set(gcf, 'KeyPressFcn', @GetKeyPress);
            n = 0;
            while ~Key_Press
                R = F.RealSuperPositionArray(n * dt);
                I = F.ImaginarySuperPositionArray(n * dt);
                P = F.DensitySuperPositionArray(n * dt);
                
                subplot(1, 2, 1);
                plot(x, P);
                title('\Psi^* \Psi');
                ylabel('\Psi^* \Psi');
                xlabel('x');
                axis([0, 1, 0, 1.2 * h]);
                
                subplot(2, 2, 2);
                plot(x, R);
                title('\Re(\Psi)');
                ylabel('\Re(\Psi)');
                xlabel('x');
                axis([0 1 -h h]);
                
                subplot(2, 2, 4);
                plot(x, I);
                title('\Im(\Psi)');
                ylabel('\Im(\Psi)');
                xlabel('x');
                axis([0 1 -h h]);
                n = n + 1;
                
                drawnow;
            end
            close
        end
        
        function VideoRealTimeEvolution(F)
            dt = F.dt;
            x = F.x;
            y = F.DensitySuperPositionArray(0);
            h = max(y);
            global Key_Press;
            Key_Press = 0;
            set(gcf, 'KeyPressFcn', @GetKeyPress);
            n = 0;
            vidObj = VideoWriter('Real.avi');
            open(vidObj);
            while ~Key_Press
                y = F.RealSuperPositionArray(n * dt);
                plot(x, y);
                title('\Re(\Psi)');
                ylabel('\Re(\Psi)');
                xlabel('x');
                axis([0 1 -h h]);
                n = n + 1;
                writeVideo(vidObj, getframe(gcf));
                drawnow;
            end
            close(gcf);
            close(vidObj);
            winopen('Real.avi');
        end
        
        function VideoImaginaryTimeEvolution(F)
            dt = F.dt;
            x = F.x;
            y = F.DensitySuperPositionArray(0);
            h = max(y);
            global Key_Press;
            Key_Press = 0;
            set(gcf, 'KeyPressFcn', @GetKeyPress);
            n = 0;
            vidObj = VideoWriter('Imaginary.avi');
            open(vidObj);
            while ~Key_Press
                y = F.ImaginarySuperPositionArray(n * dt);
                plot(x, y);
                title('\Im(\Psi)');
                ylabel('\Im(\Psi)');
                xlabel('x');
                axis([0 1 -h h]);
                n = n + 1;
                writeVideo(vidObj, getframe(gcf));
                drawnow;
            end
            close(gcf);
            close(vidObj);
            winopen('Imaginary.avi');
        end
        
        function VideoDensityTimeEvolution(F)
            dt = F.dt;
            x = F.x;
            y = F.DensitySuperPositionArray(0);
            h = 1.2 * max(y);
            global Key_Press;
            Key_Press = 0;
            set(gcf, 'KeyPressFcn', @GetKeyPress);
            n = 0;
            vidObj = VideoWriter('Density.avi');
            open(vidObj);
            while ~Key_Press
                y = F.DensitySuperPositionArray(n * dt);
                plot(x, y);
                title('\Psi^* \Psi');
                ylabel('\Psi^* \Psi');
                xlabel('x');
                axis([0 1 0 h]);
                n = n + 1;
                writeVideo(vidObj, getframe(gcf));
                drawnow;
            end
            close(gcf);
            close(vidObj);
            winopen('Density.avi');
        end
        
        function VideoTimeEvolution(F)
            dt = F.dt;
            x = F.x;
            y = F.DensitySuperPositionArray(0);
            h = max(y);
            suptitle('Particle Time Evolution');
            global Key_Press;
            Key_Press = 0;
            set(gcf, 'KeyPressFcn', @GetKeyPress);
            n = 0;
            vidObj = VideoWriter('All.avi');
            open(vidObj);
            while ~Key_Press
                R = F.RealSuperPositionArray(n * dt);
                I = F.ImaginarySuperPositionArray(n * dt);
                P = F.DensitySuperPositionArray(n * dt);
                
                subplot(1, 2, 1);
                plot(x, P);
                title('\Psi^* \Psi');
                ylabel('\Psi^* \Psi');
                xlabel('x');
                axis([0, 1, 0, 1.2 * h]);
                
                subplot(2, 2, 2);
                plot(x, R);
                title('\Re(\Psi)');
                ylabel('\Re(\Psi)');
                xlabel('x');
                axis([0 1 -h h]);
                
                subplot(2, 2, 4);
                plot(x, I);
                title('\Im(\Psi)');
                ylabel('\Im(\Psi)');
                xlabel('x');
                axis([0 1 -h h]);
                n = n + 1;
                
                writeVideo(vidObj, getframe(gcf));
                drawnow;
            end
            close(gcf);
            close(vidObj);
            winopen('All.avi');
        end
        
    end
end
            