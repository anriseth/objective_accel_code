function qarr = quantlevels(data,qlevels)
    qarr = zeros(6, length(qlevels));

    qarr(1,:) = quantile(data.ngmres_sdls,qlevels);
    qarr(2,:) = quantile(data.ngmres_sd,qlevels);
    qarr(3,:) = quantile(data.ncg,qlevels);
    qarr(4,:) = quantile(data.lbfgs,qlevels);
    qarr(5,:) = quantile(data.ngmreso_sdls,qlevels);
    qarr(6,:) = quantile(data.ngmreso_sd,qlevels);
end
