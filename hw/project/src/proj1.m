%% proj1.m
% Description: MATLAB code for course project (ICE2301, 2022 Spring)
% Author: Guorui Wei (危国锐) (313017602@qq.com; weiguorui@sjtu.edu.cn)
% Student ID: 516021910080
% Created: 2022-05-20
% Last modified: 2022-05-27

%% Initialize project

clc; clear; close all
init_env();

%% Load data.

audio_path = '../data/original_simple2.wav';
info = audioinfo(audio_path);
[x,F_s] = audioread(audio_path); % F_s: sample rate (Hz)
% sound(x,F_s)
x_ori = x(:,1); % only one single channel will be analysed
t_ori = 0:1/F_s:info.Duration;
t_ori = t_ori(1:end-1).';
t_1 = -t_ori;
t_2 = t_ori/2;
t_3 = t_ori*2;

%% Continuous-time Fourier transform, approximated by DFT.

%%% Create Fig.2

t_fig_ampl = figure('Name',"Fig.2a CTFT (amplitude) of orginal_wave_and_its_time-domain_transformation");
t_fig_pha = figure('Name',"Fig.2b CTFT (phase) of orginal_wave_and_its_time-domain_transformation");
t_TCL_ampl = tiledlayout(t_fig_ampl,2,2,"TileSpacing","tight","Padding","tight");
t_TCL_pha = tiledlayout(t_fig_pha,2,2,"TileSpacing","tight","Padding","tight");
%
[omega_ori,CTFT_ori] = CTFT_fig(t_ori,x_ori,nexttile(t_TCL_ampl,1),nexttile(t_TCL_pha,1),"$x = f(t)$");
CTFT_fig(t_1,x_ori,nexttile(t_TCL_ampl,2),nexttile(t_TCL_pha,2),"$x = f(-t)$");
CTFT_fig(t_2,x_ori,nexttile(t_TCL_ampl,3),nexttile(t_TCL_pha,3),"$x = f(2t)$");
CTFT_fig(t_3,x_ori,nexttile(t_TCL_ampl,4),nexttile(t_TCL_pha,4),"$x = f(t/2)$");
%
xlabel(t_TCL_ampl,"frequency (kHz)",FontSize=10,Interpreter="latex");
ylabel(t_TCL_ampl,"amplitude","FontSize",10,"Interpreter","latex");
[~,t_title_s] = title(t_TCL_ampl,"\bf CTFT (amplitude)","Guorui Wei 516021910080",'Interpreter','latex');
set(t_title_s,'FontSize',8);
%
xlabel(t_TCL_pha,"frequency (kHz)",FontSize=10,Interpreter="latex");
ylabel(t_TCL_pha,"phase (rad)","FontSize",10,"Interpreter","latex");
[~,t_title_s] = title(t_TCL_pha,"\bf CTFT (phase)","Guorui Wei 516021910080",'Interpreter','latex');
set(t_title_s,'FontSize',8);
%
% exportgraphics(t_TCL_ampl,"..\\doc\\fig\\proj1\\Fig_2a_CTFT_ampl_of_orginal_data_and_its_time-domain_transformation.emf",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');
exportgraphics(t_TCL_ampl,"..\\doc\\fig\\proj1\\Fig_2a_CTFT_ampl_of_orginal_data_and_its_time-domain_transformation.png",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');
% exportgraphics(t_TCL_pha,"..\\doc\\fig\\proj1\\Fig_2b_CTFT_phase_of_orginal_data_and_its_time-domain_transformation.emf",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');
exportgraphics(t_TCL_pha,"..\\doc\\fig\\proj1\\Fig_2b_CTFT_phase_of_orginal_data_and_its_time-domain_transformation.png",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');

%% Inverse continuous-time Fourier transform, approximated by inverse DFT.

mag_spectrum = abs(CTFT_ori);
pha_spectrum = CTFT_ori./mag_spectrum;

[t_ori_IFT,x_ori_IFT] = ICTFT(omega_ori,CTFT_ori,t_ori(1));
[t_mag_IFT,x_mag_IFT] = ICTFT(omega_ori,mag_spectrum,t_ori(1));
[t_pha_IFT,x_pha_IFT] = ICTFT(omega_ori,pha_spectrum,t_ori(1));

%% Ideal filter

omega_c = 2500*2*pi; % [rad/s] cut-off frequency of ideal low-pass filter (LPF) 
LPF_spec = omega_ori < omega_c | omega_ori > (2*pi*F_s - omega_c);
CTFT_filtered = CTFT_ori.*LPF_spec;
[t_filtered_IFT,x_filtered_IFT] = ICTFT(omega_ori,CTFT_filtered,t_ori(1));

%% Create figure.

%% Fig.1 orginal wave and its time-domain transformation

figure('Name',"Fig.1 orginal_wave_and_its_time-domain_transformation");
t_TCL = tiledlayout(2,2,"TileSpacing","tight","Padding","tight");
%
t_axes = nexttile(t_TCL,1);
plot(t_axes,t_ori,x_ori,'-',"DisplayName",'original waves');
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
title(t_axes,"$x = f(t)$",'Interpreter','latex');
% legend(t_axes,"Location",'best','Interpreter','latex',"Box","off",'FontSize',10);
%
t_axes = nexttile(t_TCL,2);
plot(t_axes,t_1,x_ori,'-',"DisplayName",'$x = f(-t)$');
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','YTickLabel',{},'XLimitMethod','tight');
title(t_axes,"$x = f(-t)$",'Interpreter','latex');
%
t_axes = nexttile(t_TCL,3);
plot(t_axes,t_2,x_ori,'-',"DisplayName",'$x = f(2t)$');
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
title(t_axes,"$x = f(2t)$",'Interpreter','latex');
%
t_axes = nexttile(t_TCL,4);
plot(t_axes,t_3,x_ori,'-',"DisplayName",'$x = f(t/2)$');
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','YTickLabel',{},'XLimitMethod','tight');
title(t_axes,"$x = f(t/2)$",'Interpreter','latex');
%
xlabel(t_TCL,"time (seconds)",FontSize=10,Interpreter="latex");
ylabel(t_TCL,"signal","FontSize",10,"Interpreter","latex");
[~,t_title_s] = title(t_TCL,"\bf Orginal waveform and its time-domain transformation","Guorui Wei 516021910080",'Interpreter','latex');
set(t_title_s,'FontSize',8);
%
% exportgraphics(t_TCL,"..\\doc\\fig\\proj1\\Fig_1_orginal_wave_and_its_time-domain_transformation.emf",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');
exportgraphics(t_TCL,"..\\doc\\fig\\proj1\\Fig_1_orginal_wave_and_its_time-domain_transformation.png",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');

%% Fig.3a inverse CTFT

figure('Name',"Fig.3a ICTFT: original signal");
t_TCL = tiledlayout(1,2,"TileSpacing","tight","Padding","tight");
%
t_axes = nexttile(t_TCL,1);
plot(t_axes,t_ori_IFT,x_ori_IFT,'-',"DisplayName",'ICTFT of original waves');
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
ylabel(t_axes,"signal","FontSize",10,"Interpreter","latex");
title(t_axes,"$\widetilde{x}(t) := \mathcal{F}^{-1}[X(\omega)](t)$",'Interpreter','latex');
% errors
x_err = (x_ori_IFT - x_ori)./abs(x_ori);
% x_err(abs(x_ori) < 0.01) = NaN;
t_axes = nexttile(t_TCL,2);
plot(t_axes,t_ori,x_err,'-',"DisplayName",'ICTFT of original waves');
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
ylabel(t_axes,"relative error","FontSize",10,"Interpreter","latex");
title(t_axes,"$(\widetilde{x}(t) - x(t))/|x(t)|$",'Interpreter','latex');
%
xlabel(t_TCL,"time (seconds)",FontSize=10,Interpreter="latex");
[~,t_title_s] = title(t_TCL,"\bf Inverse Fourier Transform","Guorui Wei 516021910080",'Interpreter','latex');
set(t_title_s,'FontSize',8);
%
% exportgraphics(t_TCL,"..\\doc\\fig\\proj1\\Fig_3a_ICTFT_original.emf",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');
exportgraphics(t_TCL,"..\\doc\\fig\\proj1\\Fig_3a_ICTFT_original.png",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');

%% Fig.3b inverse CTFT of magnitude spectrum

figure('Name',"Fig.3b ICTFT: magnitude spectrum");
t_TCL = tiledlayout(2,2,"TileSpacing","tight","Padding","tight");
%
t_axes = nexttile(t_TCL,1);
plot(t_axes,omega_ori/2/pi/1000,mag_spectrum,'-', ...
    "DisplayName",'$|X(\omega)|$');
hold on
plot(t_axes,F_s/2/1000,0,'o', ...
    "DisplayName",sprintf("Nyquist $f_{\\rm s}/2 = %g$ kHz",F_s/2/1000));
hold off
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
legend(t_axes,"Location",'best','Interpreter','latex',"Box","off",'FontSize',10);
xlabel(t_axes,"frequency (kHz)",FontSize=10,Interpreter="latex");
ylabel(t_axes,"magnitude","FontSize",10,"Interpreter","latex");
title(t_axes,"$|X(\omega)| := |\mathcal{F}[x(t)](\omega)|$",'Interpreter','latex');
%
t_axes = nexttile(t_TCL,2);
plot(t_axes,t_mag_IFT,real(x_mag_IFT),'-',"DisplayName",'ICTFT (real) of magnitude spectrum');
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
xlabel(t_axes,"$t$ (seconds)",FontSize=10,Interpreter="latex");
ylabel(t_axes,"real part","FontSize",10,"Interpreter","latex");
title(t_axes,"${\rm Re} \left\{ \widetilde{x}(t) \right\} := {\rm Re} \left\{ \mathcal{F}^{-1}[|X(\omega)|](t) \right\}$",'Interpreter','latex');
%
t_axes = nexttile(t_TCL,3);
plot(t_axes,t_mag_IFT,imag(x_mag_IFT),'-',"DisplayName",'ICTFT (imag.) of magnitude spectrum');
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
xlabel(t_axes,"$t$ (seconds)",FontSize=10,Interpreter="latex");
ylabel(t_axes,"imag. part","FontSize",10,"Interpreter","latex");
title(t_axes,"${\rm Im} \left\{  \widetilde{x}(t) \right\} := {\rm Im} \left\{ \mathcal{F}^{-1}[|X(\omega)|](t) \right\}$",'Interpreter','latex');
% errors
x_err = (x_mag_IFT - x_ori)./abs(x_ori);
x_err(abs(x_ori) < 0.01) = NaN;
t_axes = nexttile(t_TCL,4);
plot(t_axes,t_ori,x_err,'-',"DisplayName",'ICTFT of magnitude spectrum');
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
xlabel(t_axes,"$t$ (seconds)",FontSize=10,Interpreter="latex");
ylabel(t_axes,"relative error","FontSize",10,"Interpreter","latex");
title(t_axes,"$(\widetilde{x}(t) - x(t))/|x(t)|$",'Interpreter','latex');
%
[~,t_title_s] = title(t_TCL,"\bf Inverse Fourier Transform of Magnitude Spectrum","Guorui Wei 516021910080",'Interpreter','latex');
set(t_title_s,'FontSize',8);
%
% exportgraphics(t_TCL,"..\\doc\\fig\\proj1\\Fig_3b_ICTFT_mag_spectrum.emf",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');
exportgraphics(t_TCL,"..\\doc\\fig\\proj1\\Fig_3b_ICTFT_mag_spectrum.png",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');

%% Fig.3c inverse CTFT of phase spectrum

figure('Name',"Fig.3c ICTFT: phase spectrum");
t_TCL = tiledlayout(2,2,"TileSpacing","tight","Padding","tight");
%
t_axes = nexttile(t_TCL,1);
plot(t_axes,omega_ori/2/pi/1000,angle(CTFT_ori),'-',"DisplayName",'phase spectrum');
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
xlabel(t_axes,"frequency (kHz)",FontSize=10,Interpreter="latex");
ylabel(t_axes,"phase","FontSize",10,"Interpreter","latex");
title(t_axes,"${\rm arg} \left\{ X(\omega) \right\} := {\rm arg} \left\{ \mathcal{F}[x(t)](\omega) \right\}$",'Interpreter','latex');
%
t_axes = nexttile(t_TCL,2);
plot(t_axes,t_pha_IFT,real(x_pha_IFT),'-',"DisplayName",'ICTFT (real) of phase spectrum');
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
xlabel(t_axes,"$t$ (seconds)",FontSize=10,Interpreter="latex");
ylabel(t_axes,"real part","FontSize",10,"Interpreter","latex");
title(t_axes,"$\widetilde{x}(t) := \mathcal{F}^{-1} \left[ X(\omega)/|X(\omega)| \right](t)$",'Interpreter','latex');
%
t_axes = nexttile(t_TCL,3);
plot(t_axes,t_pha_IFT,imag(x_pha_IFT),'-',"DisplayName",'ICTFT (imag.) of phase spectrum');
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
xlabel(t_axes,"$t$ (seconds)",FontSize=10,Interpreter="latex");
ylabel(t_axes,"imag. part","FontSize",10,"Interpreter","latex");
title(t_axes,"$\widetilde{x}(t) := \mathcal{F}^{-1} \left[ X(\omega)/|X(\omega)| \right](t)$",'Interpreter','latex');
% errors
x_err = (x_pha_IFT - x_ori)./abs(x_ori);
x_err(abs(x_ori) < 0.01) = NaN;
t_axes = nexttile(t_TCL,4);
plot(t_axes,t_ori,x_err,'-',"DisplayName",'ICTFT of phase spectrum');
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
xlabel(t_axes,"$t$ (seconds)",FontSize=10,Interpreter="latex");
ylabel(t_axes,"relative error","FontSize",10,"Interpreter","latex");
title(t_axes,"$(\widetilde{x}(t) - x(t))/|x(t)|$",'Interpreter','latex');
%
[~,t_title_s] = title(t_TCL,"\bf Inverse Fourier Transform of Phase Spectrum","Guorui Wei 516021910080",'Interpreter','latex');
set(t_title_s,'FontSize',8);
%
% exportgraphics(t_TCL,"..\\doc\\fig\\proj1\\Fig_3c_ICTFT_pha_spectrum.emf",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');
exportgraphics(t_TCL,"..\\doc\\fig\\proj1\\Fig_3c_ICTFT_pha_spectrum.png",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');

%% Fig.4 waveform after passing an ideal LPF

figure('Name',"Fig.4 ideal-LPF-filtered waveform");
t_TCL = tiledlayout(2,2,"TileSpacing","tight","Padding","tight");
%
t_axes = nexttile(t_TCL,1);
plot(t_axes,omega_ori/2/pi/1000,mag_spectrum,'-', ...
    "DisplayName",'$|X(\omega)|$');
hold on
plot(t_axes,[omega_c/2/pi/1000,omega_c/2/pi/1000],[0,max(mag_spectrum)],'--', ...
    "DisplayName",sprintf("cut-off $f_{\\rm c} = %g$ kHz",omega_c/2/pi/1000));
plot(t_axes,F_s/2/1000,0,'o', ...
    "DisplayName",sprintf("Nyquist $f_{\\rm s}/2 = %g$ kHz",F_s/2/1000));
hold off
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
legend(t_axes,"Location",'best','Interpreter','latex',"Box","off",'FontSize',10);
xlabel(t_axes,"frequency (kHz)",FontSize=10,Interpreter="latex");
ylabel(t_axes,"magnitude","FontSize",10,"Interpreter","latex");
title(t_axes,"$|X(\omega)| := |\mathcal{F}[x(t)](\omega)|$",'Interpreter','latex');
%
t_axes = nexttile(t_TCL,2);
plot(t_axes,omega_ori/2/pi/1000,LPF_spec,'-', ...
    "DisplayName",'$H(f)$');
hold on
plot(t_axes,omega_c/2/pi/1000,0,'^', ...
    "DisplayName",sprintf("cut-off $f_{\\rm c} = %g$ kHz",omega_c/2/pi/1000));
plot(t_axes,F_s/2/1000,0,'o', ...
    "DisplayName",sprintf("Nyquist $f_{\\rm s}/2 = %g$ kHz",F_s/2/1000));
hold off
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
legend(t_axes,"Location",'best','Interpreter','latex',"Box","off",'FontSize',10);
xlabel(t_axes,"frequency (kHz)",FontSize=10,Interpreter="latex");
ylabel(t_axes,"magnitude","FontSize",10,"Interpreter","latex");
title(t_axes,sprintf("ideal LPF: $H(\\omega) := \\mathcal{F}[h(t)](\\omega)$"),'Interpreter','latex');
%
t_axes = nexttile(t_TCL,3);
plot(t_axes,t_filtered_IFT,x_filtered_IFT,'-',"DisplayName",'LPF-filtered waveform');
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
xlabel(t_axes,"$t$ (seconds)",FontSize=10,Interpreter="latex");
ylabel(t_axes,"amplitude","FontSize",10,"Interpreter","latex");
title(t_axes,"$\widetilde{x}(t) := \mathcal{F}^{-1}[X(\omega)H(\omega)](t)$",'Interpreter','latex');
% errors
x_err = (x_filtered_IFT - x_ori)./abs(x_ori);
x_err(abs(x_ori) < 0.01) = NaN;
t_axes = nexttile(t_TCL,4);
plot(t_axes,t_ori,x_err,'-',"DisplayName",'LPF-filtered waveform');
set(t_axes,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
xlabel(t_axes,"$t$ (seconds)",FontSize=10,Interpreter="latex");
ylabel(t_axes,"relative error","FontSize",10,"Interpreter","latex");
title(t_axes,"$(\widetilde{x}(t) - x(t))/|x(t)|$",'Interpreter','latex');
%
[~,t_title_s] = title(t_TCL,"\bf ideal-LPF-filtered waveform","Guorui Wei 516021910080",'Interpreter','latex');
set(t_title_s,'FontSize',8);
%
% exportgraphics(t_TCL,"..\\doc\\fig\\proj1\\Fig_4_LPF_filtered_signal.emf",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');
exportgraphics(t_TCL,"..\\doc\\fig\\proj1\\Fig_4_LPF_filtered_signal.png",'Resolution',800,'ContentType','auto','BackgroundColor','none','Colorspace','rgb');

%% local functions

%% Initialize environment

function [] = init_env()
% Initialize environment
%
    % set up project directory
    if ~isfolder("../doc/fig/proj1")
        mkdir ../doc/fig/proj1
    end
    % configure searching path
    mfile_fullpath = mfilename('fullpath'); % the full path and name of the file in which the call occurs, not including the filename extension.
    mfile_fullpath_without_fname = mfile_fullpath(1:end-strlength(mfilename));
    addpath(genpath(mfile_fullpath_without_fname + "../data"), ...
            genpath(mfile_fullpath_without_fname + "../inc")); % adds the specified folders to the top of the search path for the current MATLAB® session.

    return;
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

%% Tool function for Fig.2 (CTFT)

function [omega,X] = CTFT_fig(t,x,axes_amp,axes_pha,str_title)
% H1
% details
    arguments
        t
        x
        axes_amp
        axes_pha
        str_title = "$x = f(t)$"
    end

    [omega,X] = CTFT(t,x);
    %
    N_half = ceil(length(omega)/2);
    plot(axes_amp,omega(1:N_half)/2/pi/1000,abs(X(1:N_half)),'-', ...
        "DisplayName",str_title);
    set(axes_amp,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
    title(axes_amp,str_title,'Interpreter','latex');
    %
    plot(axes_pha,omega(1:N_half)/2/pi/1000,angle(X(1:N_half)),'-', ...
        "DisplayName",str_title);
    set(axes_pha,"YDir",'normal',"TickLabelInterpreter",'latex',"FontSize",10,'Box','off','TickDir','out','XLimitMethod','tight');
    title(axes_pha,str_title,'Interpreter','latex');
end

%% ICTFT by IDFT

function [t,x] = ICTFT(omega,X,t0)
% inverse continuous-time Fourier transform, approximated by IDFT.
% CAVEAT: MUST be sampled at equal intervals!
    arguments
        omega
        X
        t0 = 0;
    end

    if(~iscolumn(omega))
       omega = omega.'; 
    end
    if(~iscolumn(X))
       X = X.'; 
    end
    
    [omega,I] = sort(omega);
    X = X(I);

    N = length(omega);
    omega_s = (omega(2) - omega(1))*N;
    T_s = 2*pi/omega_s;
    
    t = t0 + (0:N-1).'*T_s;
    x = ifft(X.*exp(1i*omega*t0))/T_s;
end
