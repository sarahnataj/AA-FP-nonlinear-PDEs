function plotformat(width,msize)
set(0,'DefaultaxesLineWidth', width) 
set(0,'DefaultaxesFontSize', 12) 
set(0,'DefaultaxesFontWeight', 'bold') 
set(0,'DefaultTextFontSize', 12) 
set(0,'DefaultaxesFontName', 'Times new Roman') 
set(0,'DefaultlegendFontName', 'Times new Roman')
%set(0,'defaultAxesXGrid','on') 
%set(0,'defaultAxesYGrid','on') 
set(0,'DefaultLineMarkerSize',msize);
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
end