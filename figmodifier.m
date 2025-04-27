clc
clear
close all
NumberUncertain = 4
clear InitialGuess IG
% IG =1; ttl = '(a)';%InitialGuess = 0*ones(Ntime,1);
% IG =2; ttl = '(b)';%InitialGuess = H_range(1)*ones(Ntime,1);
% IG =3; ttl = '(c)';%InitialGuess = H_range(2)*ones(Ntime,1);
% IG =4; ttl = '(d)';%InitialGuess = H_range(3)*ones(Ntime,1);
% IG =5; ttl = '(e)';%InitialGuess = H_range(4)*ones(Ntime,1);
% IG =6; ttl = '(f)';%InitialGuess = H_range(5)*ones(Ntime,1);
% IG =7; ttl = '(g)';%indices = [3 5 1 2 4 2 3 5 1 3 2 4 5 2 1 4 5 3 1 5 4 2 3 1 5 2 3 1 5 2]; InitialGuess = H_range(indices);  % 1. Random variation of H bounded by the minimum and maximum values
IG =8; ttl = '(h)'; %[H_var_steps,] = H_variation_between_H_AND_0(Ntime); InitialGuess = 7600*H_var_steps;    % Rob's group H variation Peak amplitude was 7600 [A/m] with 50% duty cycle (60s On and 60s Off)
savepath = 'C:\Users\MNandyala\OneDrive - Inside MD Anderson\Documents\GitHub\fmincon_Entropy_arbitraryFunction\';
savename2 = sprintf('IG%dNU%d_Fe3O4NP_OptH.fig',IG,NumberUncertain)
hFig = openfig([savepath,savename2]);
% Resize the figure window
newWidth = 400; % Width in pixels
newHeight = 275; % Height in pixels
newX = 100; % X-position of the lower-left corner of the figure
newY = 100; % Y-position of the lower-left corner of the figure
set(hFig, 'Position', [newX, newY, newWidth, newHeight]); % Update the figure size
grid off
ylim([-10000,60000])
ylabel([{'Magnetic field'} {'amplitude, H (A/m)'}])
title(ttl,FontWeight="normal")
legendHandle = legend('show','Location','best');
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman', 'FontSize', 16);
% Change the font size of the legend
set(legendHandle, 'FontSize', 14); % Replace 14 with your desired font size
% saveas(gcf,[savepath,savename2],'fig')
% saveas(gcf,[savepath,savename2],'png')
% Save the figure as a .fig file
% figFileName = 'example.fig'; % Replace 'example.fig' with your desired .fig file name
% savefig(figureHandle, figFileName);