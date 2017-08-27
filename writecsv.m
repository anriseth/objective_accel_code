function writecsv(data,fname)
    qlevels = [0.1,0.5,0.9];
    qarr = quantlevels(data,qlevels);
    csvwrite(fname,qarr);
end

