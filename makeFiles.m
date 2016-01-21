function [ output_args ] = makeFiles( file )

%.eps black and white
print('-deps',file)

%.eps color
fileColor   = strcat(file,'c');
print('-depsc',fileColor);

%pdf

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
print('-dpdf', file);


end

