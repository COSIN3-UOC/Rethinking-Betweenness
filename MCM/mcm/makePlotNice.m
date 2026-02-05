% Set up some constants for changing font sizes, etc
if(~exist('scale','var'))
    scale = 1;
end

AxisFontsize=20*scale;
TextFontsize=20*scale;
ArrowWidth=10*scale;
titleFontSize=25*scale;
lineWidth = 4*scale;

% set some properties
set(gca, 'FontName', 'Arial');
set(gca,'fontsize',AxisFontsize,'xgrid','off','ygrid','off')

h = get(gca, 'title');
set(h, 'FontName', 'Arial');
set(h,'fontsize',titleFontSize,'fontweight','bold');

%h = findobj(gcf, 'type', 'line');
% set the linewidth for all lines
%set(h,'linewidth',lineWidth);
% set the linewidth for selected line(s)
h=get(gca,'xlabel');
set(h, 'FontSize', TextFontsize) 

h=get(gca,'ylabel');
set(h, 'FontSize', TextFontsize) 