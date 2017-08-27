function retstr = latextable(quantiles, names, order)
    retstr = '';
    for ii = 1:length(order)
        retstr = strcat(retstr, sprintf('%s&%.1f&%.1f&%.1f\\\\\\\\\\n', names{order(ii)}, ...
                                        quantiles(order(ii),1), ...
                                        quantiles(order(ii),2), quantiles(order(ii),3)));
    end
    
end

