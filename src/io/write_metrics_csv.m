function write_metrics_csv(metrics, cfg, filename)
%WRITE_METRICS_CSV Flatten component peak metrics to CSV.

    mechanism = {};
    order_col = [];
    process = {};
    peak_meV = [];
    fwhm_meV = [];
    hwhm_meV = [];

    mechs = {'optical','piezoelectric'};
    for im = 1:numel(mechs)
        mk = mechs{im};
        for order = cfg.oap.photon_orders
            ok = sprintf('order_%d',order);
            proc = {'emission','absorption'};
            for ip = 1:2
                item = metrics.(mk).(ok).(proc{ip});
                mechanism{end+1,1} = mk; %#ok<AGROW>
                order_col(end+1,1) = order; %#ok<AGROW>
                process{end+1,1} = proc{ip}; %#ok<AGROW>
                peak_meV(end+1,1) = item.peak_x; %#ok<AGROW>
                fwhm_meV(end+1,1) = item.fwhm; %#ok<AGROW>
                hwhm_meV(end+1,1) = item.hwhm; %#ok<AGROW>
            end
        end
    end

    T = table(mechanism,order_col,process,peak_meV,fwhm_meV,hwhm_meV, ...
        'VariableNames',{'Mechanism','PhotonOrder','Process', ...
                         'Peak_meV','FWHM_meV','HWHM_meV'});
    writetable(T,filename);
end
