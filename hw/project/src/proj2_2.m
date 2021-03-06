%% proj2_2.m
% Description: MATLAB code for course project (ICE2301, 2022 Spring)
% Author: Guorui Wei (危国锐) (313017602@qq.com; weiguorui@sjtu.edu.cn)
% Student ID: 516021910080
% Created: 2022-05-26
% Last modified: 2022-05-28

%% Initialize project

clc; clear; close all
init_env();

%% Create original signal.

F_ori = 10000./(2.^(0:2));      % frequency of components of original signal [MUST be geometric series!]
phi_ori = [pi/6,pi/2,5*pi/6];   % initial phase of components of origina signal
F_s = 64*max(F_ori);            % sampling frequency
tau_s = 64/min(F_ori);          % sampling duration

t_o = linspace(0,tau_s-1/F_s,F_s*tau_s).'; % time-interval
x_o = x_ori(t_o,F_ori,phi_ori);            % original signal
[omega_ori,CTFT_ori] = CTFT(t_o,x_o);      % CTFT of original signal
cycle_ori = 1/min(F_ori);                  % cycle of original signal

%%% Create figure.

figure('Name',"Fig.0 original signal");
t_TCL = tiledlayout(2,2,"TileSpacing","tight","Padding","tight");
% CTFT of original signal (mag.)
t_axes = nexttile(t_TCL,1);
plot(t_axes,omega_ori/2/pi/1000,abs(CTFT_ori),'-', ...
    "DisplayName","$|X_{\rm o}(f)|$");
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out');
xlim(t_axes,[0,2*max(F_ori)/1000]);
xlabel(t_axes,"frequency (kHz)",FontSize=10,Interpreter="latex");
ylabel(t_axes,"magnitude","FontSize",10,"Interpreter","latex");
title(t_axes,"$X_{\rm o}(\omega) := \mathcal{F}[x_{\rm o}(t)](\omega)$",'Interpreter','latex');
% CTFT of original signal (pha.)
t_axes = nexttile(t_TCL,2);
plot(t_axes,omega_ori/2/pi/1000,angle(CTFT_ori),'-', ...
    "DisplayName","${\rm arg}\{X_{\rm o}(f)\}$");
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLim',[0,2*max(F_ori)/1000]);
xlabel(t_axes,"frequency (kHz)",FontSize=10,Interpreter="latex");
ylabel(t_axes,"phase (rad)","FontSize",10,"Interpreter","latex");
title(t_axes,"$X_{\rm o}(\omega) := \mathcal{F}[x_{\rm o}(t)](\omega)$",'Interpreter','latex');
% original signal
t_axes = nexttile(t_TCL,3,[1 2]);
plot(t_axes,t_o(1:2*F_s*cycle_ori+1)*1e6,x_o(1:2*F_s*cycle_ori+1),'-', ...
    "DisplayName","$x_{\rm o}(t)$");
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
xlabel(t_axes,"$t$ $(\rm \mu s)$",FontSize=10,Interpreter="latex");
ylabel(t_axes,"amplitude","FontSize",10,"Interpreter","latex");
title(t_axes,sprintf("{\\bf original signal $x_{\\rm o}(t)$} (periodic, $T = %g$ ${\\rm \\mu s}$)",1e6/min(F_ori)),'Interpreter','latex');
%
[~,t_title_s] = title(t_TCL,sprintf("\\bf original signal (periodic)"),"Guorui Wei 516021910080",'Interpreter','latex');
set(t_title_s,'FontSize',8);
%
exportgraphics(t_TCL,sprintf("..\\doc\\fig\\proj2_2\\Fig_0_original_signal.emf"),'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');
exportgraphics(t_TCL,sprintf("..\\doc\\fig\\proj2_2\\Fig_0_original_signal.png"),'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');

%% Recover signal by ideal band-limited interpolation.

%%% parameters

F_s = 44e3./(2.^(0:4)); % sampling frequency
t_s0 = 0;               % sampling start time
tau_s = cycle_ori*128;  % sampling duration
F_r = 32*max(F_ori);    % recovering frequency
t_r0 = t_s0;            % recovering start time
tau_r = tau_s;          % recovering duration

for i = 1:length(F_s)
    ibli_fig(F_s(i),t_s0,tau_s,F_r,t_r0,tau_r,F_ori,phi_ori);
end

%% local functions

%% Initialize environment

function [] = init_env()
% Initialize environment
%
    % set up project directory
    if ~isfolder("../doc/fig/proj2_2")
        mkdir ../doc/fig/proj2_2
    end
    % configure searching path
    mfile_fullpath = mfilename('fullpath'); % the full path and name of the file in which the call occurs, not including the filename extension.
    mfile_fullpath_without_fname = mfile_fullpath(1:end-strlength(mfilename));
    addpath(genpath(mfile_fullpath_without_fname + "../data"), ...
            genpath(mfile_fullpath_without_fname + "../inc")); % adds the specified folders to the top of the search path for the current MATLAB® session.
    return;
end

%% Create original signal.

function [x] = x_ori(t,F_ori,phi_ori)
%x_ori - original signal
%
% Syntax: [x] = x_ori(t,F_ori,phi_ori)
%
% create original signal.
    arguments
        t = 0;
        F_ori = 1;
        phi_ori = 0;
    end

    x = zeros(size(t));
    for i = 1:length(F_ori)
        x = x + sin(2*pi*F_ori(i)*t + phi_ori(i));
    end
end

%% ideal band-limited interpolation

function x_r = ibli(t_n,x_n,t_r)
%ibli - ideal band-limited interpolation
%
% Syntax: x_r = ibli(t_n,x_n,t_r)
%
% ideal band-limited interpolation
    arguments
        t_n % sample time
        x_n % sample value
        t_r % recover value
    end

    if(~iscolumn(x_n))
        x_n = x_n.';
    end

    T_s = abs(t_n(2) - t_n(1));
    [grid_t_r,grid_t_n] = ndgrid(t_r,t_n);
    x_r = sinc((grid_t_r - grid_t_n)/T_s)*x_n;
end


%% CTFT by DFT

function [omega,X] = CTFT(t,x)
% Continuous-time Fourier transform, approximated by DFT.
% CAVEAT: MUST be sampled at equal intervals!
    arguments
        t
        x
    end

    if(~iscolumn(t))
        t = t.'; 
    end
    if(~iscolumn(x))
        x = x.'; 
    end
    
    [t,I] = sort(t);
    x = x(I);
    
    T_s = t(2) - t(1);
    N = length(t);
    
    omega = (0:N-1).'/N*2*pi/T_s;
    X = T_s * fft(x) .* exp(-omega*t(1)*1i);
end

%% Recover signal, by ideal band-limited interpolation

function [] = ibli_fig(F_s,t_s0,tau_s,F_r,t_r0,tau_r,F_ori,phi_ori)
    arguments
        F_s = 22000;            % sampling frequency
        t_s0 = 0;               % sampling start time
        tau_s = 1/F_s;          % sampling duration
        F_r = F_s;              % recovering frequency
        t_r0 = t_s0;            % recovering start time
        tau_r = tau_s;          % recovering duration
        F_ori = [10000,5000];   % frequency of components of original signal [MUST be geometric series!]
        phi_ori = [0,0];        % initial phase of components of origina signal
    end

    cycle_ori = 1/min(F_ori);                       % cycle of original signal
    t_n = linspace(t_s0,t_s0+tau_s,F_s*tau_s+1).';  % sampling time
    x_n = x_ori(t_n,F_ori,phi_ori);                 % sampling values
    t_r = linspace(t_r0,t_r0+tau_r,F_r*tau_r+1).';  % recovering time
    x_r = ibli(t_n,x_n,t_r);                        % recovering values
    x_r_ori = x_ori(t_r,F_ori,phi_ori);             % original values at recovering time
    [omega_r,CTFT_r] = CTFT(t_r,x_r);               % CTFT of recovering signal
    % recovering error
    err_r = (x_r - x_r_ori)./abs(x_r_ori);
    err_r(abs(x_r_ori) < 0.05) = NaN;
    
    %%% Create figure.
    
    figure('Name',"Fig.1 ideal band-limited interpolation");
    t_TCL = tiledlayout(2,2,"TileSpacing","tight","Padding","tight");
    % CTFT of recovering signal (mag.)
    t_axes = nexttile(t_TCL,2);
    plot(t_axes,omega_r/2/pi/1000,abs(CTFT_r),'-', ...
        "DisplayName","$|X_{\rm r}(f)|$");
    hold on
    plot(t_axes,F_s/2/1000,0,'o', ...
        "DisplayName",sprintf("Nyquist $f_{\\rm s}/2 = %g$ kHz",F_s/2/1000));
    hold off
    set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out');
    xlim(t_axes,[0,2*max(F_ori)/1000]);
    legend(t_axes,"Location",'best','Interpreter','latex',"Box","off",'FontSize',10);
    xlabel(t_axes,"frequency (kHz)",FontSize=10,Interpreter="latex");
    ylabel(t_axes,"magnitude","FontSize",10,"Interpreter","latex");
    title(t_axes,sprintf("$X_{\\rm r}(\\omega) := \\mathcal{F}[x_{\\rm r}(t)](\\omega),$ $N_{\\rm r} = %g$",length(t_r)),'Interpreter','latex');
    % CTFT of recovering signal (pha.)
    t_axes = nexttile(t_TCL,4);
    plot(t_axes,omega_r/2/pi/1000,angle(CTFT_r),'-', ...
        "DisplayName","${\rm arg}\{X_{\rm r}(f)\}$");
    set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out');
    xlim(t_axes,[0,2*max(F_ori)/1000]);
%     legend(t_axes,"Location",'best','Interpreter','latex',"Box","off",'FontSize',10);
    xlabel(t_axes,"frequency (kHz)",FontSize=10,Interpreter="latex");
    ylabel(t_axes,"phase (rad)","FontSize",10,"Interpreter","latex");
    title(t_axes,"$X_{\rm r}(\omega) := \mathcal{F}[x_{\rm r}(t)](\omega)$",'Interpreter','latex');
    % recovering signal
    TF_t_r_range = t_r >= t_s0-cycle_ori & t_r <= t_s0+tau_s+cycle_ori; % only show recovering signal at sampling time-interval
    t_axes = nexttile(t_TCL,1);
    plot(t_axes,t_r(TF_t_r_range)*1e3,x_r(TF_t_r_range),'-',"DisplayName","$x_{\rm r}(t)$");
    set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
    xlabel(t_axes,"$t$ $\rm (ms)$",FontSize=10,Interpreter="latex");
    ylabel(t_axes,"amplitude","FontSize",10,"Interpreter","latex");
    title(t_axes,"\bf recovering signal $x_{\rm r}(t)$",'Interpreter','latex');
    % recovering errors
    t_axes = nexttile(t_TCL,3);
    plot(t_axes,t_r*1e3,err_r,'-',"DisplayName",'LPF-filtered waveform');
    set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
    xlabel(t_axes,"$t$ $\rm (ms)$",FontSize=10,Interpreter="latex");
    ylabel(t_axes,"relative error","FontSize",10,"Interpreter","latex");
    title(t_axes,"$e_{\rm r}(t) := (x_{\rm r}(t) - x(t))/|x(t)|$",'Interpreter','latex');
    %
    [~,t_title_s] = title(t_TCL,sprintf("{\\bf ideal band-limited interpolation} ($f_{\\rm s} = %g$ kHz, $N_{\\rm s} = %g$)",F_s/1000,length(t_n)),"Guorui Wei 516021910080",'Interpreter','latex');
    set(t_title_s,'FontSize',8);
    %
    % exportgraphics(t_TCL,sprintf("..\\doc\\fig\\proj2_2\\Fig_%s_kHz_ideal_band_limited_interp.emf",num2str(F_s/1000)),'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');
    exportgraphics(t_TCL,sprintf("..\\doc\\fig\\proj2_2\\Fig_%s_kHz_ideal_band_limited_interp.png",num2str(F_s/1000)),'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');
end
