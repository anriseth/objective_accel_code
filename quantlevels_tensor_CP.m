function qarr = quantlevels_tensor_CP(data,qlevels)
    qarr = zeros(4, length(qlevels));

    qarr(1,:) = quantile(data.ngmres_als,qlevels);
    qarr(2,:) = quantile(data.ncg,qlevels);
    qarr(3,:) = quantile(data.lbfgs,qlevels);
    qarr(4,:) = quantile(data.ngmreso_als,qlevels);
    qarr(5,:) = quantile(data.als,qlevels);
end
