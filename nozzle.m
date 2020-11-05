% Main program
clear all;close all;clc;
n  = 31; % given set of number of nodes  
x = linspace(0,3,n);
mdot_exact = 0.579*ones(1,n);
 
% This 'for loop' ensures to call the function for the length of n(i.e, 4
% times) to calculate velocity for different grid spacing
tic
[rho_non,V_non,T_non,P_non,M_non,mdot_nonconser,rho_th_nonconser,P_th_nonconser,V_th_nonconser,T_th_nonconser,nt,rho_nonconser_throat,P_nonconser_throat,T_nonconser_throat,M_nonconser_throat] = nozzleweek7nonconser(n,x);
nozzleweek7nonconser = toc;
 
tic
[rho,V,T,P,M,mdot_conser,rho_th_con,P_th_con,V_th_con,T_th_con,nt,rho_conser_throat,P_conser_throat,T_conser_throat,M_conser_throat] = nozzleweek7conser(n,x);
nozzleweek7conser = toc;
 
% Exact solution at Throat
rho_exact = 0.639;
T_exact = 0.833;
P_exact = 0.528;
M_exact = 1;
 
% Figure showing the comparison of both the governing equation forms
figure(7)
subplot(5,1,1)
plot(x,P_non)    
title_text = sprintf('Variation of primitive variables w.r.t nozzle length in both the forms of governing equations');
title(title_text,'Fontsize',12)
hold on
grid on
plot(x,P,'*')
legend({'\bf Non-conservative','\bf Conservative'},'Fontsize',12)
xlabel('\bf x/L','Fontsize',12)
ylabel('\bf P/P_0','Fontsize',12)
 
subplot(5,1,2)
plot(x,T_non)
hold on
grid on
plot(x,T,'*')
legend({'\bf Non-conservative','\bf Conservative'},'Fontsize',12)
xlabel('\bf x/L','Fontsize',12)
ylabel('\bf T/T_0','Fontsize',12)
 
subplot(5,1,3)
plot(x,rho_non)
hold on
grid on
plot(x,rho,'-*')
legend({'\bf Non-conservative','\bf Conservative'},'Fontsize',12)
xlabel('\bf x/L','Fontsize',12)
ylabel('\bf rho/rho_0','Fontsize',12)
 
subplot(5,1,4)
plot(x,M_non)
hold on
grid on
plot(x,M,'*')
legend({'\bf Non-conservative','\bf Conservative'},'Location','northwest','Fontsize',12)
xlabel('\bf x/L','Fontsize',12)
ylabel('\bf M','Fontsize',12)
 
subplot(5,1,5)
plot(x,V_non)
hold on
grid on
plot(x,V,'*')
legend({'\bf Non-conservative','\bf Conservative'},'Location','northwest','Fontsize',12)
xlabel('\bf x/L','Fontsize',12)
ylabel('\bf V/a_0','Fontsize',12)
 
% Figure showing the mass flow rate variation
figure(8)
plot(x,mdot_nonconser)
hold on
grid on
plot(x,mdot_conser,'*')
plot(x,mdot_exact,'k','linewidth',1.5)
legend({'\bf Non-conservative','\bf Conservative'},'Location','north','Fontsize',12)
xlabel('\bf x/L','Fontsize',12)
ylabel('\bf mdot','Fontsize',12)
title_text = sprintf('Variation of mass flow rate w.r.t nozzle length in both the forms of governing equations');
title(title_text,'Fontsize',12)
 
% Print values
fprintf('Number of nodes = %f',n)
fprintf('\n');
fprintf('rho_throat_nonconservative = %g',rho_nonconser_throat)
fprintf('\n');
fprintf('P_throat_nonconservative = %g',P_nonconser_throat)
fprintf('\n');
fprintf('T_throat_nonconservative = %g',T_nonconser_throat)
fprintf('\n');
fprintf('M_throat_nonconservative = %g',M_nonconser_throat)
fprintf('\n');
fprintf('\n');
 
fprintf('rho_throat_conservative = %g',rho_conser_throat)
fprintf('\n');
fprintf('P_throat_conservative = %g',P_conser_throat)
fprintf('\n');
fprintf('T_throat_conservative = %g',T_conser_throat)
fprintf('\n');
fprintf('M_throat_conservative = %g',M_conser_throat)
fprintf('\n');
fprintf('\n');
 
fprintf('rho_analytical = %g',rho_exact)
fprintf('\n');
fprintf('P_analytical = %g',P_exact)
fprintf('\n');
fprintf('T_analytical = %g',T_exact)
fprintf('\n');
fprintf('M_analytical = %g',M_exact)
fprintf('\n');



% Program to simulate isentropic flow through quasi 1D subsonic-supersonic nozzle using Maccormack method for 
% conservative form of governing equations
 
function [rho,V,T,P,M,mdot_conser,rho_th_con,P_th_con,V_th_con,T_th_con,nt,rho_conser_throat,P_conser_throat,T_conser_throat,M_conser_throat] = nozzleweek7conser(n,x);
 
    % Inputs for mesh
    x = linspace(0,3,n);
    dx = x(2)-x(1);
 
    % Initial conditions
        for i = 1:length(x)
            if x(i)>=0 && x(i)<=0.5
                    rho(i) = 1; T(i) = 1;
                elseif x(i)>=0.5 && x(i)<=1.5
                    rho(i) = 1-0.366*(x(i)-0.5); T(i) = 1-0.167*(x(i)-0.5);
                elseif x(i)>=1.5 && x(i)<=3
                    rho(i) = 0.634-0.3879*(x(i)-1.5); T(i) = 0.833-0.3507*(x(i)-1.5);
            end
        end
 
    A = 1+2.2*(x-1.5).^2; % Area
    V = 0.59./(rho.*A); % Velocity
    P  = rho.*T; % Pressure
    gamma = 1.4; % Ratio of specific heats for Air
    mdot_conser = rho.*A.*V; % Mass flow rate
    
    % Throat
    throat = (n+1)*0.5;
    
    % Solution vector
    U1 = rho.*A;
    U2 = rho.*A.*V;
    U3 = rho.*A.*((T/(gamma-1)) + ((gamma/2)*V.^2));
 
    % Time steps
    nt = 1400;
    dt = 0.0206;
 
    % Outer Time loop. This loop is for time-marching. Since the governing
    % equations are unsteady, this loop is justified.
    for k = 1:nt
 
        % Solution update
        U1_old = U1;
        U2_old = U2;
        U3_old = U3;
 
        % Predictor step
        for j = 2:n-1 
            dAdx = (A(j+1)-A(j))/dx;
            F1 = U2;
            F2 = ((U2.^2./U1) + ((gamma-1)/gamma)*(U3-(gamma/2)*(U2.^2./U1)));
            F3 = (gamma*(U3.*U2./U1)) - ((gamma*(gamma-1)/2)*(U2.^3./U1.^2));
            J2(j) = (1/gamma)*(rho(j).*T(j)).*dAdx;
 
            % Continuity Equation
            dU1_dt_p(j) = -(F1(j+1) - F1(j))/dx;
 
            % Momentum Equation
            dU2_dt_p(j) = -((F2(j+1) - F2(j))/dx) + J2(j);
 
            % Energy Equation
            dU3_dt_p(j) = -(F3(j+1) - F3(j))/dx;
 
            % Solution Update after Predictor step
            U1(j) = U1(j) + dU1_dt_p(j)*dt;
            U2(j) = U2(j) + dU2_dt_p(j)*dt;
            U3(j) = U3(j) + dU3_dt_p(j)*dt;            
        end 
 
        % Corrector step
        for j = 2:n-1 
            dAdx = (A(j)-A(j-1))/dx;
            F1 = U2;
            F2 = ((U2.^2./U1) + ((gamma-1)/gamma)*(U3-(gamma/2)*(U2.^2./U1)));
            F3 = (gamma*(U3.*U2./U1)) - ((gamma*(gamma-1)/2)*(U2.^3./U1.^2));        
            J2(j) = (1/gamma)*(rho(j).*T(j)).*dAdx;
 
            % Continuity Equation
            dU1_dt_c(j) = -(F1(j) - F1(j-1))/dx;
 
            % Momentum Equation
            dU2_dt_c(j) = -((F2(j) - F2(j-1))/dx) + J2(j);
 
            % Energy Equation
            dU3_dt_c(j) = -(F3(j) - F3(j-1))/dx;                   
        end
 
        % Compute the average values
        dU1dt = 0.5*(dU1_dt_p + dU1_dt_c);
        dU2dt = 0.5*(dU2_dt_p + dU2_dt_c);
        dU3dt = 0.5*(dU3_dt_p + dU3_dt_c);
 
 
        % Final solution update
        for i = 2:n-1
            U1(i) = U1_old(i) + dU1dt(i)*dt;
            U2(i) = U2_old(i) + dU2dt(i)*dt;
            U3(i) = U3_old(i) + dU3dt(i)*dt;
        end
 
        % Applying BCS
        % Inlet
        U1(1) = rho(1)*A(1);
        U2(1) = 2*U2(2) - U2(3);
        V(1) = U2(1)/U1(1);
        U3(1) = (U1(1).*((T(1)/(gamma-1)) + ((gamma/2)*V(1).^2))); 
 
        % Outlet
        U1(n) = 2*U1(n-1) - U1(n-2);
        U2(n) = 2*U2(n-1) - U2(n-2);
        U3(n) = 2*U3(n-1) - U3(n-2);
 
        % Primitive variables
        rho = U1./A;
        V = U2./U1;
        T = ((gamma-1)*((U3./U1)-(gamma/2)*(U2./U1).^2));
 
        % Other varaiables
        P = rho.*T;
        M = V./sqrt(T); % Mach number
        mdot_conser = rho.*A.*V; % Mass flow rate
         
        % Figure showing Time-wise variation of Mass flow rate
        figure(4)
        if k == 100
                plot(x,mdot_conser,'r','linewidth',2)
                hold on
                grid on
            elseif k == 200
                plot(x,mdot_conser,'b','linewidth',2)
            elseif k == 500
                plot(x,mdot_conser,':bs','linewidth',2)
            elseif k == 700
                plot(x,mdot_conser,'k','linewidth',3)
            elseif k == 1000
                plot(x,mdot_conser,'--yo','linewidth',2)
            elseif k == 1400
                plot(x,mdot_conser,'-.m*','linewidth',1)
        end
        xlabel('\bf Distance','Fontsize',12)
        ylabel('\bf Mass flow rate','Fontsize',12)
        title_text = sprintf('Time-wise variation of Mass flow rate \n Conservative Governing Equation \n No. of nodes = %f',n);
        title(title_text)
        legend({'100 timesteps','200 timesteps','500 timesteps','700 timesteps','1000 timesteps','1400 timesteps'},'Fontsize',12)
 
        % Values at Throat
        rho_th_con(k) = rho(throat);
        T_th_con(k) = T(throat);
        P_th_con(k) = P(throat);
        V_th_con(k) = V(throat);
    end
 
    % Exact solution at Throat
    rho_exact = 0.639*ones(1,k);
    T_exact = 0.833*ones(1,k);
    P_exact = 0.528*ones(1,k);
    M_exact = ones(1,n);
    
    % Throat Values
    rho_conser_throat = rho(throat);
    P_conser_throat = P(throat);
    T_conser_throat = T(throat);
    M_conser_throat = M(throat);
 
    % Figure showing the variation of primitive variables w.r.t nozzle
    % length 
    figure(5)
    plot(x,P,'r','linewidth',2)
    hold on
    grid on
    plot(x,rho,'b','linewidth',2)
    plot(x,T,'k','linewidth',2)
    plot(x,M,':ms','linewidth',2)
    plot(x,V,'--bo','linewidth',2)
    plot(x,mdot_conser,'*','linewidth',2)
    xlabel('\bf x/L','Fontsize',12)
    ylabel('\bf Non-dimensional variables','Fontsize',12)
    legend({'\bf P/P_0','\bf rho/rho_0 ','\bf T/T_0','\bf M','\bf V/a_0','\bf mdot'},'Location','northwest')
    title_text = sprintf('Conservative Governing Equations \n MacCormack Method \n Iterations = %d \n No. of nodes = %f', k,n);
    title(title_text,'Fontsize',12)
 
    % Timewise variation of primitive variables
    figure(6)
    subplot(2,2,1)
    plot(P_th_con,'r','linewidth',2)
    grid on
    xlabel('\bf Number of Itertions','Fontsize',12)
    ylabel('\bf P_{th}/P_0','Fontsize',12)
    title_text = sprintf('Conservative Governing Equations \n Time-wise Variation');
    title(title_text,'Fontsize',12)
    
    subplot(2,2,2)
    plot(rho_th_con,'b','linewidth',2)
    grid on
    xlabel('\bf Number of Itertions','Fontsize',12)
    ylabel('\bf rho_{th}/rho_0','Fontsize',12)
    title_text = sprintf('Conservative Governing Equations \n Time-wise Variation');
    title(title_text,'Fontsize',12)
    
    subplot(2,2,3)
    plot(T_th_con,'k','linewidth',2)
    grid on
    xlabel('\bf Number of Itertions','Fontsize',12)
    ylabel('\bf T_{th}/T_0','Fontsize',12)
    
    subplot(2,2,4)
    plot(V_th_con,':m','linewidth',2)
    grid on
    xlabel('\bf Number of Itertions','Fontsize',12)
    ylabel('\bf V_{th}/a_0','Fontsize',12)
end


% Program to simulate isentropic flow through quasi 1D subsonic-supersonic nozzle using Maccormack method for 
% non-conservative form of governing equations
 
function [rho_non,V_non,T_non,P_non,M_non,mdot_nonconser,rho_th_nonconser,P_th_nonconser,V_th_nonconser,T_th_nonconser,nt,rho_nonconser_throat,P_nonconser_throat,T_nonconser_throat,M_nonconser_throat] = nozzleweek7nonconser(n,x);
   
    % Inputs for mesh
    x = linspace(0,3,n);
    dx = x(2)-x(1);
 
    % Initial conditions
    rho_non = 1-0.3146*x; % Density
    T_non = 1-0.2314*x; % Temperature
    V_non = (0.1+1.09*x).*sqrt(T_non); % Velocity
    gamma = 1.4; % Ratio of specific heats for Air
 
    % Area of the nozzle
    A = 1+2.2*(x-1.5).^2;
    
    % Throat
    throat = (n+1)*0.5;
    
    % Time steps
    nt = 1400; % Number of timesteps
    dt = 0.0206; % Tiem step size
 
    % Outer Time loop. This loop is for time-marching. Since the governing
    % equations are unsteady, this loop is justified.
    for k = 1:nt
 
        % This is to update each variable after each timestep
        rho_old = rho_non;
        V_old = V_non;
        T_old = T_non;
 
        % Predictor step
        for j = 2:n-1         
           dVdx = (V_non(j+1) - V_non(j))/dx;
           drhodx = (rho_non(j+1)-rho_non(j))/dx;
           dlnAdx = (log(A(j+1))-log(A(j)))/dx;
           dTdx = (T_non(j+1)-T_non(j))/dx;
 
           % Continuity Equation
           drho_dt_p(j) = -rho_non(j)*dVdx - rho_non(j)*V_non(j)*dlnAdx - V_non(j)*drhodx;
 
           % Momentum Equation
           dV_dt_p(j) = -V_non(j)*dVdx - (1/gamma)*(dTdx+(T_non(j)/(rho_non(j)))*drhodx);
 
           % Energy Equation
           dT_dt_p(j) = -V_non(j)*dTdx - (gamma-1)*T_non(j)*(dVdx + V_non(j)*dlnAdx);   
 
           % Solution Update after Predictor step
           rho_non(j) = rho_non(j) + drho_dt_p(j)*dt;
           V_non(j) = V_non(j) + dV_dt_p(j)*dt;
           T_non(j) = T_non(j) + dT_dt_p(j)*dt;            
        end 
 
        % Corrector step
        for j = 2:n-1
           dVdx = (V_non(j) - V_non(j-1))/dx;
           drhodx = (rho_non(j)-rho_non(j-1))/dx;
           dlnAdx = (log(A(j))-log(A(j-1)))/dx;
           dTdx = (T_non(j)-T_non(j-1))/dx;
 
           % Continuity Equation
           drho_dt_c(j) = -rho_non(j)*dVdx - rho_non(j)*V_non(j)*dlnAdx - V_non(j)*drhodx;
 
           % Momentum Equation
           dV_dt_c(j) = -V_non(j)*dVdx - (1/gamma)*(dTdx+(T_non(j)/(rho_non(j)))*drhodx);
 
           % Energy Equation
           dT_dt_c(j) = -V_non(j)*dTdx - (gamma-1)*T_non(j)*(dVdx + V_non(j)*dlnAdx);   
        end 
 
        % Compute the average values
        drhodt = 0.5*(drho_dt_p + drho_dt_c);
        dVdt = 0.5*(dV_dt_p + dV_dt_c);
        dTdt = 0.5*(dT_dt_p + dT_dt_c);
 
        % Final solution update
        for j = 2:n-1
            rho_non(j) = rho_old(j) + drhodt(j)*dt;
            V_non(j) = V_old(j) + dVdt(j)*dt;
            T_non(j) = T_old(j) + dTdt(j)*dt;
        end
 
        % Applying BCS
        % Inlet
        V_non(1) = 2*V_non(2) - V_non(3);
        % Outlet
        rho_non(n) = 2*rho_non(n-1) - rho_non(n-2);
        V_non(n) = 2*V_non(n-1) - V_non(n-2);
        T_non(n) = 2*T_non(n-1) - T_non(n-2);
 
        % Other quantities
        M_non = V_non./sqrt(T_non); % Mach number
        mdot_nonconser = rho_non.*A.*V_non; % Mass flow rate
        P_non = rho_non.*T_non; % Pressure
 
        % Figure showing Time-wise variation of Mass flow rate
        figure(1)
        if k == 100
                plot(x,mdot_nonconser,'r','linewidth',2)
                hold on
                grid on
            elseif k == 200
                plot(x,mdot_nonconser,'b','linewidth',2)
            elseif k == 500
                plot(x,mdot_nonconser,':bs','linewidth',2)
            elseif k == 700
                plot(x,mdot_nonconser,'k','linewidth',3)
            elseif k == 1000
                plot(x,mdot_nonconser,'--yo','linewidth',2)
            elseif k == 1400
                plot(x,mdot_nonconser,'-.m*','linewidth',1)
        end
        xlabel('\bf Distance','Fontsize',12)
        ylabel('\bf Mass flow rate','Fontsize',12)
        title_text = sprintf('Time-wise variation of Mass flow rate \n Non-conservative Governing Equation \n No. of nodes = %f',n);
        title(title_text)
        legend({'100 timesteps','200 timesteps','500 timesteps','700 timesteps','1000 timesteps','1400 timesteps'},'Fontsize',12)   
 
        % Values at Throat
        rho_th_nonconser(k) = rho_non(throat);
        T_th_nonconser(k) = T_non(throat);
        P_th_nonconser(k) = P_non(throat);
        V_th_nonconser(k) = V_non(throat);
    end
 
    % Exact solution at Throat
    rho_exact = 0.639*ones(1,k);
    T_exact = 0.833*ones(1,k);
    P_exact = 0.528*ones(1,k);
    M_exact = ones(1,k);
    
    % Throat Values
    rho_nonconser_throat = rho_non(throat);
    P_nonconser_throat = P_non(throat);
    T_nonconser_throat = T_non(throat);
    M_nonconser_throat = M_non(throat);
 
    % Figure showing the variation of primitive variables w.r.t nozzle
    % length 
    figure(2)
    plot(x,P_non,'r','linewidth',2)
    hold on
    grid on
    plot(x,rho_non,'b','linewidth',2)
    plot(x,T_non,'k','linewidth',2)
    plot(x,M_non,':ms','linewidth',2)
    plot(x,V_non,'--bo','linewidth',2)
    plot(x,mdot_nonconser,'*','linewidth',2)
    xlabel('\bf x/L','Fontsize',12)
    ylabel('\bf Non-dimensional variables','Fontsize',12)
    legend({'\bf P/P_0','\bf rho/rho_0 ','\bf T/T_0','\bf M','\bf V/a_0','\bf mdot'},'Location','northwest','Fontsize',12)
    title_text = sprintf('Non-conservative Governing Equations \n MacCormack Method \n Iterations = %d \n No. of nodes = %f', k,n);
    title(title_text,'Fontsize',12)
    
    % Timewise variation of primitive variables
    figure(3)
    subplot(2,2,1)
    plot(P_th_nonconser,'r','linewidth',2)
    grid on
    xlabel('\bf Number of Itertions')
    ylabel('\bf P_{th}/P_0')
    title_text = sprintf('Non-conservative Governing Equations \n Time-wise Variation');
    title(title_text,'Fontsize',12)
    
    subplot(2,2,2)
    plot(rho_th_nonconser,'b','linewidth',2)
    grid on
    xlabel('\bf Number of Itertions','Fontsize',12)
    ylabel('\bf rho_{th}/rho_0','Fontsize',12)
    title_text = sprintf('Non-conservative Governing Equations \n Time-wise Variation');
    title(title_text,'Fontsize',12)
 
    subplot(2,2,3)
    plot(T_th_nonconser,'k','linewidth',2)
    grid on
    xlabel('\bf Number of Itertions','Fontsize',12)
    ylabel('\bf T_{th}/T_0','Fontsize',12)
 
    subplot(2,2,4)
    plot(V_th_nonconser,':m','linewidth',2)
    grid on
    xlabel('\bf Number of Itertions','Fontsize',12)
    ylabel('V_{th}/a_0','Fontsize',12) 
end
