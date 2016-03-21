 

function plotMotifs(x ,loc1, loc2, lengs)

 W = 4;

 fh = figure();
    for i = 1:length(lengs)

        if mod(i-1,W) == 0

            figure;

        end

        subplot(W/2,2,mod(i-1,W)+1);

        a = loc1(i)+1;

        b = loc2(i)+1;

        l = lengs(i);

        line(a:a+l,zNorm(x(a:a+l)),'Color','r');
        ax1 = gca;
        set(ax1,'XColor','r','YColor','r')
        xlim(ax1,[a,a+l]);
        ylim(ax1,[-5,5]);
      %  set(ax1,'XTickLabel',num2str(get(ax1,'XTick').'))
        xlabel(strcat('Length : ',num2str(l)),'Color','black');
        
        ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','b','YColor','b');
     
         line(b:b+l,zNorm(x(b:b+l)) , 'Color','b','Parent',ax2);
        xlim(ax2,[b,b+l]);
        ylim(ax2,[-5,5]);

%set(ax2,'XTickLabel',num2str(get(ax2,'XTick').'))
        

    end

end