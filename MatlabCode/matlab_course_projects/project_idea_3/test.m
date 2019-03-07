% Program 5d: Animated bifurcation diagram for a SFR resonator.
% The parameter kappa (related to fibre coupling) is increased. 
clear
halfN=7999;
N=2*halfN+1;
N1=1+halfN;
E1=zeros(1,N);
E2=zeros(1,N); 
Esqr=zeros(1,N);
Esqr1=zeros(1,N);
ptsup=zeros(1,N); 
for j = 1:60 % amount of frames, typically 60
    F(j) = getframe;
    format long;
    C=.345913; % change from .345913 to .5, rapid bifurcation
    E1(1)=0;
    kappa=.001*j; % .001 keeps the plot in the center, but higher values make the plot hard to see, lower ...
    % ... values make bifurcation harder to produce on the frame scale
    Pmax=120; % change from 120 to 70... the lower the value the higher off the plot the bifurcation starts
    phi=0; % change phi from 0 to pi/4
    % Ramp the power up 
    for n=1:halfN
        E2(n+1)=E1(n)*exp(1i*(abs(C*E1(n))^2-phi)); 
        E1(n+1)=1i*sqrt(1-kappa)*sqrt(n*Pmax/N1)+sqrt(kappa)*E2(n+1); 
        Esqr(n+1)=abs(E1(n+1))^2;
    end
    % Ramp the power down 
    for n=N1:N
        E2(n+1)=E1(n)*exp(1i*(abs(C*E1(n))^2-phi)); 
        E1(n+1)=1i*sqrt(1-kappa)*sqrt(2*Pmax-n*Pmax/N1)+sqrt(kappa)*E2(n+1);
        Esqr(n+1)=abs(E1(n+1))^2; 
    end
    for n=1:halfN
        Esqr1(n)=Esqr(N+1-n); 
        ptsup(n)=n*Pmax/N1;
    end
    % Plot the bifurcation diagrams 
    fsize=15;
    hold 
    plot(ptsup(1:halfN), Esqr(1:halfN), '.', 'MarkerSize', 1, 'Color', 'r');
    plot(ptsup(1:halfN), Esqr1(1:halfN), '.', 'MarkerSize', 1, 'Color', 'b');
    xlabel("Input power");
    ylabel("Output power");
end
movie(F,5)