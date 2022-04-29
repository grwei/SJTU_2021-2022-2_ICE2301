%% ch4.m
% Description: MATLAB code for Homework (Chapter 4) (ICE2301, 2022 Spring)
% Author: Guorui Wei (危国锐) (313017602@qq.com; weiguorui@sjtu.edu.cn)
% Student ID: 516021910080
% Created: 2022-04-17
% Last modified: 2022-04-24

%% Initialize project

clc; clear; close all
init_env();

%% (Chapter 4-1) Question 2

k = -10:10;
%
figure("Name","ch4-1 Q2_ak")
t_TCL = tiledlayout(2,1,"TileSpacing","tight","Padding","tight");
%
t_Axes_abs = nexttile(t_TCL,1);
t_plot_abs = stem(t_Axes_abs,k,abs(Q2_a(k)),'-o',"DisplayName",'amplitude','Color','#0072BD');
set(t_Axes_abs,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off');
xlabel(t_Axes_abs,"$k$","Interpreter",'latex');
ylabel(t_Axes_abs,"amplitude","Interpreter",'latex');
title(t_Axes_abs,"\bf Amplitude",'Interpreter','latex')
grid on
% legend(t_plot_abs,"Location",'best','Interpreter','latex',"Box","off",'FontSize',10);
%
t_Axes_phase = nexttile(t_TCL,2);
t_plot_phase = stem(t_Axes_phase,k,rad2deg(angle(Q2_a(k))),'-o',"DisplayName",'phase','Color','#D95319');
set(t_Axes_phase,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off');
xlabel(t_Axes_phase,"$k$","Interpreter",'latex');
ylabel(t_Axes_phase,"phase (deg)","Interpreter",'latex');
title(t_Axes_phase,"\bf Phase",'Interpreter','latex')
grid on
% legend(t_plot_phase,"Location",'best','Interpreter','latex',"Box","off",'FontSize',10);
%
[~,t_title_s] = title(t_TCL,"\bf 2022 Spring ICE2301 Homework (Chapter 4-1) Q2","Guorui Wei 516021910080","Interpreter",'latex');
set(t_title_s,'FontSize',8)
%
exportgraphics(t_TCL,"..\\doc\\fig\\hw4_1_Q2_a_k.emf",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')
exportgraphics(t_TCL,"..\\doc\\fig\\hw4_1_Q2_a_k.png",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')

%% (Chapter 4-1) Question 5-1

T = 4;
T1 = 1/4;
t0 = 1;
t = linspace(-1/2*T-t0,3/2*T-t0,1601);
%
figure("Name","ch4-1 Q5_1")
t_TCL = tiledlayout(1,1,"TileSpacing","tight","Padding","tight");
%
t_Axes = nexttile(t_TCL,1);
t_plot_abs = plot(t_Axes,t,Q5_1_x(t),'-',"DisplayName",'amplitude','Color','#0072BD');
set(t_Axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off');
xlabel(t_Axes,"$t$","Interpreter",'latex');
ylabel(t_Axes,"$x(t)$","Interpreter",'latex');
title(t_Axes,"\bf $T_1 = 1/4, \, T = 4$",'Interpreter','latex')
grid on
% legend(t_plot_abs,"Location",'best','Interpreter','latex',"Box","off",'FontSize',10);
[~,t_title_s] = title(t_TCL,"\bf 2022 Spring ICE2301 Homework (Chapter 4-1) Q5(a)","Guorui Wei 516021910080","Interpreter",'latex');
set(t_title_s,'FontSize',8)
%
exportgraphics(t_TCL,"..\\doc\\fig\\hw4_1_Q5a_x.emf",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')
exportgraphics(t_TCL,"..\\doc\\fig\\hw4_1_Q5a_x.png",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')

%% (Chapter 4-3) Question 5(1)

clc; clear;

omega_c = 1;
omega_0 = 2*omega_c;
t_0 = 1;
t = t_0 + linspace(-2*pi/omega_c,2*pi/omega_c,1001);
ch4_3_Q5_h = @(t) 2/pi*cos(omega_0*t).*sin(omega_c*(t-t_0))./(t-t_0);

%
figure("Name","ch4-3 Q5(a)")
t_TCL = tiledlayout(1,1,"TileSpacing","tight","Padding","tight");
%
t_Axes = nexttile(t_TCL,1);
hold on
t_plot_h = plot(t_Axes,t,ch4_3_Q5_h(t),'-',"DisplayName",'$h(t)$');
t_plot_env1 = plot(t_Axes,t,2/pi*cos(omega_0*t),'--',"DisplayName",'$\frac{\pi}{2} \cos{(\omega_0 t)}$');
t_plot_env2 = plot(t_Axes,t,-2/pi*cos(omega_0*t),'--',"DisplayName",'$-\frac{\pi}{2} \cos{(\omega_0 t)}$');
hold off
set(t_Axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off');
xlabel(t_Axes,"$t$","Interpreter",'latex');
ylabel(t_Axes,"","Interpreter",'latex');
title(t_Axes,sprintf("$\\displaystyle h(t) = \\frac{2}{\\pi} \\cos{(\\omega_0 t)} \\frac{\\sin{(\\omega_{\\rm{c}} (t-t_0))}}{t-t_0}, \\quad \\omega_{\\rm{c}} = %d, \\, \\omega_0 = %d, \\, t_0 = %d$",omega_c,omega_0,t_0),'Interpreter','latex')
grid on
legend(t_Axes,"Location",'best','Interpreter','latex',"Box","off",'FontSize',10);
[~,t_title_s] = title(t_TCL,"\bf 2022 Spring ICE2301 Homework (Chapter 4-3) Q5(a)","Guorui Wei 516021910080","Interpreter",'latex');
set(t_title_s,'FontSize',8)
%
exportgraphics(t_TCL,"..\\doc\\fig\\hw4_3_Q5_1.emf",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')
exportgraphics(t_TCL,"..\\doc\\fig\\hw4_3_Q5_1.png",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')

%% (Chapter 4-3) Question 12(a)

ch4_3_Q12a_H = @(omega) 1/2./(omega*1i+1) + 3/2./(omega*1i+3);
omega = linspace(-8,8,1001);

figure("Name","ch4-3 Q12a")
t_TCL = tiledlayout(1,1,"TileSpacing","tight","Padding","tight");
%
t_Axes = nexttile(t_TCL,1);
yyaxis(t_Axes,"left")
t_plot_abs = plot(t_Axes,omega,abs(ch4_3_Q12a_H(omega)),'-',"DisplayName",'amplitude','Color','#0072BD');
set(t_Axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off');
% xlabel(t_Axes,"$\omega$","Interpreter",'latex');
ylabel(t_Axes,"amplitude","Interpreter",'latex');
title(t_Axes,"\bf $\displaystyle H(\omega) = \frac{1/2}{\rm{j} \omega + 1} + \frac{3/2}{\rm{j} \omega + 3}$",'Interpreter','latex')
%
yyaxis(t_Axes,"right")
t_plot_phase = plot(t_Axes,omega,rad2deg(angle(ch4_3_Q12a_H(omega))),'-',"DisplayName",'phase','Color','#D95319');
set(t_Axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off');
xlabel(t_Axes,"$\omega$","Interpreter",'latex');
ylabel(t_Axes,"phase (deg)","Interpreter",'latex');
% title(t_Axes,"\bf Phase",'Interpreter','latex')
grid on
legend(t_Axes,"Location",'best','Interpreter','latex',"Box","off",'FontSize',10);
%
[~,t_title_s] = title(t_TCL,"\bf 2022 Spring ICE2301 Homework (Chapter 4-3) Q12(a)","Guorui Wei 516021910080","Interpreter",'latex');
set(t_title_s,'FontSize',8)
%
exportgraphics(t_TCL,"..\\doc\\fig\\hw4_3_Q12a_H.emf",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')
exportgraphics(t_TCL,"..\\doc\\fig\\hw4_3_Q12a_H.png",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb')

%% local functions

%% Initialize environment
function [] = init_env()
% Initialize environment
%
    % set up project directory
    if ~isfolder("../doc/fig/")
        mkdir ../doc/fig/
    end
    % configure searching path
    mfile_fullpath = mfilename('fullpath'); % the full path and name of the file in which the call occurs, not including the filename extension.
    mfile_fullpath_without_fname = mfile_fullpath(1:end-strlength(mfilename));
    addpath(genpath(mfile_fullpath_without_fname + "../data"), ...
            genpath(mfile_fullpath_without_fname + "../inc")); % adds the specified folders to the top of the search path for the current MATLAB® session.

    return;
end

%% Q2_a
function [Q2_a] = Q2_a(k)
% Q2
%
    arguments
        k
    end

    Q2_a = 1i./(2.*k*pi).*(-1).^k + ((-1).^k - 1)./(2*(k*pi).^2);
    Q2_a(k == 0) = 1/4;
end

%% Q5_1
function [x] = Q5_1_x(t,t0,T1,T,E)
%
%
    arguments
        t
        t0 = 1;
        T1 = 1/4;
        T = 4;
        E = 1;
    end

    x = Q5_1_g(t+t0,T1,T,E)/2 - 1/16;
end

function [g] = Q5_1_g(t,T1,T,E)
%
%
    arguments
        t
        T1
        T
        E
    end

    g = zeros(size(t));
    g(mod(t,T) <= T1 | T - mod(t,T) < T1) = E;
end
