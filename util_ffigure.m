function util_ffigure(par,out)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function makes a convergence plot comparing several optimization
%methods
%
%by Hans De Sterck, September 2011
% Edited by Asbjørn Nilsen Riseth, April 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nCurves=0;
    % first determine the minimum value of f (because we will plot the
    % convergence of f to the minimum value achieved in any of the
    % tests)
    minval = Inf;
    if par.compareNGMRES_sdls==1
        minval=min(minval,min(out.out_ngmres_sdls.logf));
    end
    if par.compareNGMRES_sd==1
        minval=min(minval,min(out.out_ngmres_sd.logf));
    end
    if par.compareNCG==1
        minval=min(minval,min(out.out_ncg.TraceFunc(2:end)));
    end
    if par.compare_sdls==1
        minval=min(minval,min(out.out_sdls.logf));
    end
    if par.compareLBFGS==1
        minval=min(minval,min(out.out_lbfgs.TraceFunc(2:end)));
    end
    if par.compareNGMRESO_sdls==1
        minval=min(minval,min(out.out_ngmreso_sdls.logf));
    end
    if par.compareNGMRESO_sd==1
        minval=min(minval,min(out.out_ngmreso_sd.logf));
    end
    if par.compareNGMRESO_ALS==1
        minval=min(minval,min(out.out_ngmreso_als.logf));
    end
    if par.compareNGMRES_ALS==1
        minval=min(minval,min(out.out_ngmres_als.logf));
    end
    if par.compareALS==1
        minval=min(minval,min(out.out_als.logf));
    end
    
    if par.compareNGMRES_sdls==1
        var = out.out_ngmres_sdls;
        cnum = 1;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        nCurves=nCurves+1;
        legendstr{nCurves}='N-GMRES-sdls';
        hold on
    end
    if par.compareNGMRES_sd==1
        var = out.out_ngmres_sd;
        cnum = 2;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        nCurves=nCurves+1;
        legendstr{nCurves}='N-GMRES-sd';
        hold on
    end
    if par.compareNCG==1
        var = out.out_ncg;
        cnum = 4;
        mind = 1:par.mispace{cnum}:length(var.TraceFunc(2:end));
        semilogy(var.TraceFunc(2:end)-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        nCurves=nCurves+1;
        legendstr{nCurves}='N-CG';
        hold on
    end
    if par.compare_sdls==1
        var = out.out_sdls;
        cnum = 5;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        nCurves=nCurves+1;
        legendstr{nCurves}='sdls';
        hold on
    end
    if par.compareLBFGS==1
        var = out.out_lbfgs;
        cnum = 7;
        mind = 1:par.mispace{cnum}:length(var.TraceFunc(2:end));
        semilogy(var.TraceFunc(2:end)-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        nCurves=nCurves+1;
        legendstr{nCurves}='L-BFGS';
        hold on
    end
    if par.compareNGMRESO_sdls==1
        var = out.out_ngmreso_sdls;
        cnum = 8;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        nCurves=nCurves+1;
        legendstr{nCurves}='N-GMRES-O-sdls';
        hold on
    end
    if par.compareNGMRESO_sd==1 
        var = out.out_ngmreso_sd;
        cnum = 9;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        nCurves=nCurves+1;
        legendstr{nCurves}='N-GMRES-O-sd';
        hold on
    end
    if par.compareNGMRESO_ALS==1 
        var = out.out_ngmreso_als;
        cnum = 10;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        nCurves=nCurves+1;
        legendstr{nCurves}='N-GMRES-O-ALS';
        hold on
    end
    if par.compareNGMRES_ALS==1 
        var = out.out_ngmres_als;
        cnum = 11;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        nCurves=nCurves+1;
        legendstr{nCurves}='N-GMRES-ALS';
        hold on
    end
    if par.compareALS==1 
        var = out.out_als;
        cnum = 11;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        nCurves=nCurves+1;
        legendstr{nCurves}='ALS';
        hold on
    end

    % optionally indicate the restart points in red for the N-GMRES runs
    if par.par_ngmres.logRestart & par.restartfigure
        if par.compareNGMRES_sdls==1
            logfMod=out.out_ngmres_sdls.logf-minval;
            logfMod(out.out_ngmres_sdls.logRestart==0)=NaN;
            semilogy(logfMod+par.epsi,'or')
        end
        if par.compareNGMRES_sd==1
            logfMod=out.out_ngmres_sd.logf-minval;
            logfMod(out.out_ngmres_sd.logRestart==0)=NaN;
            semilogy(logfMod+par.epsi,'+r')
        end
        % TODO: Add logrestart plots for NGMRES{,O}_ALS
    end
    hold off
    legend(legendstr)
    xlabel('iterations')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
