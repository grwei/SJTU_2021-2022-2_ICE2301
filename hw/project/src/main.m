%% main.m
% Description: MATLAB code for course project (ICE2301, 2022 Spring)
% Author: Guorui Wei (危国锐) (313017602@qq.com; weiguorui@sjtu.edu.cn)
% Student ID: 516021910080
% Created: 2022-05-20
% Last modified: 2022-

%% Initialize project

clc; clear; close all
init_env();

%%

audio_path = '../data/original.wav';
info = audioinfo(audio_path);
[y,Fs] = audioread(audio_path);
% sound(y,Fs)

t = 0:seconds(1/Fs):seconds(info.Duration);
t = t(1:end-1);
plot(t,y)
xlabel('Time')
ylabel('Audio Signal')

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
