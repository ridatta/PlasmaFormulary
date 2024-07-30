function formatPlots(varargin)
if nargin == 0
    size = 0;
    aspect = 1;
else if nargin == 1
    size = varargin{1};
    aspect=1;
    else
        size = varargin{1};
        aspect = varargin{2};
end
    
% Default plot properties
fnt = 'Arial';
fnt_size = 24;
axis_linewidth = 3; 
linewidth = 2;
marker_size = 20; 

if size ~= 0
set(gcf,'Position',[0 0 size size/aspect]); 
end
%set(gcf, 'Position', [pos(1) pos(2) width*100 height*100]); %<- Set size 
ax = gca;
        set(gca, 'Fontname',fnt);
        set(gca, 'Fontsize',fnt_size);
     set(gca, 'LineWidth',axis_linewidth);
%         set(gca,'XMinorTick','on','YMinorTick','on')
        legend('Fontsize',fnt_size-3);
        legend boxoff
        box on
        % axis square; 
     
end