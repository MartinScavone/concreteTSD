function y = Make_Signal(Type, Size)
%
% Generates the four test signal used by:
% Donoho and Johnstone(1994) "Ideal spatial adaptation by wavelet
% shrinkage" Biometrika, 81, 3, pp. 425-455.
%
% Inputs:
% - Type: Signal type ('Blocks', 'Bumps', 'Heavisine', 'Doppler')
% - Size: Number of data points sampled from the signal
%
% Output:
% - y: The signal
%
% Written by Samer Katicha 04/01/2013
% skaticha@vtti.vt.edu

x = linspace(0,1,Size);

switch Type
    
    case 'Blocks'
        t = [0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81];
        h = [4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2,];
        y = zeros(size(x));
        for i=1:length(t)
            y = y+h(i)*(1+sign(x-t(i)))/2;
        end
        
    case 'Bumps'
        t = [0.1, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81];
        h = [4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 2.1, 4.2];
        w = [.005, .005, .006, .01, .01, .03, .01, .01, .005, .008, .005];
        y = zeros(size(x));
        for i=1:length(t)
            y = y+h(i)*(1+abs((x-t(i))/w(i))).^(-4);
        end
        
    case 'Heavisine'
        y = 4*sin(4*pi*x)-sign(x-0.3)-sign(.72-x);
        
    case 'Doppler'
        e = 0.05;
        y = sqrt(x.*(1-x)).*sin(2*pi*(1+e)./(x+e));
        
    case 'Spikes'
        p = randperm(Size);
        t = x(p(1:100));
        h = 8*ones(size(t));
        y = zeros(size(linspace(0,1,Size+1)));
        for i=1:length(t)
            y = y+h(i)*(1+sign(linspace(0,1,Size+1)-t(i)))/2;
        end
        y = diff(y);
        
    case 'Hat'
        y = -1+1.5*x+0.2*normpdf(x,0.5,0.02);
        
    case 'Harmonic'
        y = exp(-30*(1-x)).*cos(30*pi*(1-x));
        
    case 'Rapid'
        y = 0.8./(1+exp(-75*(x-0.8)));
        
    otherwise
        warning('Invalid input for signal Type; Zero signal returned')
        y = zeros(size(x));
        
end
y = y(:);
if ~strcmp(Type,'Hat')
    if ~strcmp(Type,'Harmonic')
        if ~strcmp(Type,'Rapid')
            if std(y)~=0
                y = 7*y/std(y);
            end
        end
    end
end