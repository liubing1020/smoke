function pscatter(ps, w, d1, d2, lu)
	f = figure;

    hold on;
    %colorbar;
    scatter(ps(:,d1),ps(:,d2),15,w,'s');

    xmin = lu(1,d1);
    ymin = lu(1,d2);
    xmax = lu(2,d1);
    ymax = lu(2,d2);
    axis([xmin xmax ymin ymax]);
    xlabel(['p' num2str(d1)],'FontWeight','bold', 'FontSize', 30);
    ylabel(['p' num2str(d2)],'FontWeight','bold', 'FontSize', 30);
end