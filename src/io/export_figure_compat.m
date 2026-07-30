function export_figure_compat(fig, filename, dpi)
%EXPORT_FIGURE_COMPAT Save figure in old and new MATLAB releases.
    [folder,~,~] = fileparts(filename);
    ensure_directory(folder);

    if exist('exportgraphics','file') == 2
        exportgraphics(fig, filename, 'Resolution', dpi);
    else
        set(fig,'PaperPositionMode','auto');
        print(fig, filename, '-dpng', sprintf('-r%d',dpi));
    end
end
