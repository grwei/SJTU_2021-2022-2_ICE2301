%% ch4_1.m
% Description: MATLAB code for Homework (Chapter 4-1) (ICE2301, 2022 Spring)
% Author: Guorui Wei (危国锐) (313017602@qq.com; weiguorui@sjtu.edu.cn)
% Student ID: 516021910080
% Created: 2022-04-17
% Last modified: 2022-04-

%% Initialize project

clc; clear; close all
init_env();

%% Question 2

k = -10:10;
%
figure("Name","Q2_ak")
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

%% Question 5-1

T = 4;
T1 = 1/4;
t0 = 1;
t = linspace(-1/2*T-t0,3/2*T-t0,1601);
%
figure("Name","Q2_ak")
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
