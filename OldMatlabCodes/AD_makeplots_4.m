function makeplots_4;

% Must have run BT_Analysis and saved the files 

    %Strain Ratios vs Pressure
    figure
    mm=plot(STLP(:,4),D(:,4),'color','b','linewidth',2);
    hold on
    xy=plot(STLP(:,4),D(:,9),'linewidth',2,'color','k');
    [maxs,locs]=max(D(:,1));
    plot(STLP(:,4),D(:,4)+D(:,5),'Color',[.5 .5 1])
    plot(STLP(:,4),D(:,4)-D(:,5),'Color',[.5 .5 1])
    plot(STLP(:,4),D(:,9)+D(:,10),'Color',[.5 .5 .5])
    plot(STLP(:,4),D(:,9)-D(:,10),'Color',[.5 .5 .5])
    plot(STLP(ploc,4),D(ploc,4),'ko','MarkerFaceColor','k','MarkerEdgeColor','k');
    plot(STLP(locs,4),D(locs,4),'rs','MarkerFaceColor','r','MarkerEdgeColor','r');
    plot(STLP(ploc,4),D(ploc,9),'ko','MarkerFaceColor','k','MarkerEdgeColor','k');
    plot(STLP(locs,4),D(locs,9),'rs','MarkerFaceColor','r','MarkerEdgeColor','r');
    title('Minor/Major Logarithmic Strain Ratio','Fontsize',14)
    xlabel ('Pressure (psi)','Fontsize',14)
    l=legend([mm xy],{'e_2/e_1','e_t_r/e_r_o_l_l'});
    set(l,'Location','Southeast')