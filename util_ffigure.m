function util_ffigure(par,out)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this function makes a convergence plot comparing several optimization
%methods
%
%by Hans De Sterck, September 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nCurves=1;
    % first determine the minimum value of f (because we will plot the
    % convergence of f to the minimum value achieved in any of the tests)
    switch par.figRunFirst
      case 1
        minval=min(out.out_ngmres_sdls_precond.logf);
      case 2
        minval=min(out.out_ngmres_descent_precond.logf);
      case 4
        minval=min(out.out_ncg.TraceFunc(2:end));
      case 5
        minval=min(out.out_sdls_precond.logf);
      case 7
        minval=min(out.out_lbfgs.TraceFunc(2:end));
      case 8
        minval=min(out.out_ngmreso_sdls_precond.logf);
      case 9
        minval=min(out.out_ngmreso_descent_precond.logf);
    end
    if par.compareNGMRES_sdls_precond==1 & par.figRunFirst~=1
        minval=min(minval,min(out.out_ngmres_sdls_precond.logf));
    end
    if par.compareNGMRES_descent_precond==1 & par.figRunFirst~=2
        minval=min(minval,min(out.out_ngmres_descent_precond.logf));
    end
    if par.compareNCG==1 & par.figRunFirst~=4
        minval=min(minval,min(out.out_ncg.TraceFunc(2:end)));
    end
    if par.compare_sdls_precond==1 & par.figRunFirst~=5
        minval=min(minval,min(out.out_sdls_precond.logf));
    end
    if par.compareLBFGS==1 & par.figRunFirst~=7
        minval=min(minval,min(out.out_lbfgs.TraceFunc(2:end)));
    end
    if par.compareNGMRESO_sdls_precond==1 & par.figRunFirst~=8
        minval=min(minval,min(out.out_ngmreso_sdls_precond.logf));
    end
    if par.compareNGMRESO_descent_precond==1 & par.figRunFirst~=9
        minval=min(minval,min(out.out_ngmreso_descent_precond.logf));
    end

    % then make the plots
    switch par.figRunFirst
      case 1
        var = out.out_ngmres_sdls_precond
        cnum = 1;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        legendstr{nCurves}='N-GMRES-sdls';
      case 2
        var = out.out_ngmres_descent_precond;
        cnum = 2;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        legendstr{nCurves}='N-GMRES-sd';
      case 4
        var = out.out_ncg;
        cnum = 4;
        mind = 1:par.mispace{cnum}:length(var.TraceFunc(2:end));
        semilogy(var.TraceFunc(2:end)-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        legendstr{nCurves}='N-CG';
      case 5
        var = out.out_sdls_precond;
        cnum = 5;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        legendstr{nCurves}='sdls';
      case 7
        var = out.out_lbfgs;
        cnum = 7;
        mind = 1:par.mispace{cnum}:length(var.TraceFunc(2:end));
        semilogy(var.TraceFunc(2:end)-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        legendstr{nCurves}='L-BFGS';
      case 8
        var = out.out_ngmreso_sdls_precond;
        cnum = 8;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        legendstr{nCurves}='N-GMRES-O-sdls';
      case 9
        var = out.out_ngmreso_descent_precond;
        cnum = 9;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        legendstr{nCurves}='N-GMRES-O-sd';
    end
    if par.compareNGMRES_sdls_precond==1 & par.figRunFirst~=1
        hold on
        var = out.out_ngmres_sdls_precond
        cnum = 1;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        hold off
        nCurves=nCurves+1;
        legendstr{nCurves}='N-GMRES-sdls';
    end
    if par.compareNGMRES_descent_precond==1 & par.figRunFirst~=2
        hold on
        var = out.out_ngmres_descent_precond
        cnum = 2;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        hold off
        nCurves=nCurves+1;
        legendstr{nCurves}='N-GMRES-sd';
    end
    if par.compareNCG==1 & par.figRunFirst~=4
        hold on
        var = out.out_ncg;
        cnum = 4;
        mind = 1:par.mispace{cnum}:length(var.TraceFunc(2:end));
        semilogy(var.TraceFunc(2:end)-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        hold off
        nCurves=nCurves+1;
        legendstr{nCurves}='N-CG';
    end
    if par.compare_sdls_precond==1 & par.figRunFirst~=5
        hold on
        var = out.out_sdls_precond;
        cnum = 5;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        hold off
        nCurves=nCurves+1;
        legendstr{nCurves}='sdls';
    end
    if par.compareLBFGS==1 & par.figRunFirst~=7
        hold on
        var = out.out_lbfgs;
        cnum = 7;
        mind = 1:par.mispace{cnum}:length(var.TraceFunc(2:end));
        semilogy(var.TraceFunc(2:end)-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        hold off
        nCurves=nCurves+1;
        legendstr{nCurves}='L-BFGS';
    end
    if par.compareNGMRESO_sdls_precond==1 & par.figRunFirst~=8
        hold on
        var = out.out_ngmreso_sdls_precond;
        cnum = 8;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        hold off
        nCurves=nCurves+1;
        legendstr{nCurves}='N-GMRES-O-sdls';
    end
    if par.compareNGMRESO_descent_precond==1 & par.figRunFirst~=9
        hold on
        var = out.out_ngmreso_descent_precond;
        cnum = 9;
        mind = 1:par.mispace{cnum}:length(var.logf);
        semilogy(var.logf-minval+par.epsi, ...
                 par.fs{cnum},'LineWidth',1,'MarkerIndices',mind)
        hold off
        nCurves=nCurves+1;
        legendstr{nCurves}='N-GMRES-O-sd';
    end

    % optionally indicate the restart points in red for the N-GMRES runs
    if par.par_ngmres.logRestart & par.restartfigure
        if par.compareNGMRES_sdls_precond==1
            logfMod=out.out_ngmres_sdls_precond.logf-minval;
            logfMod(out.out_ngmres_sdls_precond.logRestart==0)=NaN;
            hold on
            semilogy(logfMod+par.epsi,'or')
            hold off
        end
        if par.compareNGMRES_descent_precond==1
            logfMod=out.out_ngmres_descent_precond.logf-minval;
            logfMod(out.out_ngmres_descent_precond.logRestart==0)=NaN;
            hold on
            semilogy(logfMod+par.epsi,'+r')
            hold off
        end
    end
    legend(legendstr)
    xlabel('iterations')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%