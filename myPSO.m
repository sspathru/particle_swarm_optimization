function[fxbest,xbest] = myPSO(funchandle,xrange,dim,c1,c2,omega)


xrange = xrange*ones(1,dim);

xrange = abs(xrange);
lb = (-1.*xrange)./2;
ub = (xrange./2);
n = 2000;
iter = 1;

phi = c1+c2;
k = 2/(abs(2-phi-sqrt(phi^2 - 4*phi)));

samples = (ub-lb).*rand(n,max(size(xrange))) + lb;
vid_previous = (ub-lb).*rand(n,max(size(xrange)))+ lb;


xpbest_previous = samples;
samples_previous = samples;

fit_samples = funchandle(samples);
fit_samples_previous = fit_samples;
dis = 1;

while dis > 0.0000000001 && iter<2000

    fit_samples = funchandle(samples);

    [fgbest(iter,:),idxgbest] = min(fit_samples);
    xgbest  = samples(idxgbest,:);
    
    if iter > 500
        dis = abs(((abs(fgbest(iter-500,:)) - min(abs(fgbest(iter-500:iter,:))))));
    end
    
    det = fit_samples < fit_samples_previous;
    a = det.*samples;
    b = ~det.*xpbest_previous;
    xpbest = a+b ;

    samples_previous = samples;
    fit_samples_previous = funchandle(samples_previous);

    vid = (omega*vid_previous + c1.*rand(n,max(size(xrange))).*(xpbest - samples) + c2.*rand(n,max(size(xrange))).*(xgbest - samples));
    samples = samples + vid;

    xpbest_previous = xpbest;

    vid_previous = vid;

    iter = iter+1;

end
disp(iter)
xbest = xgbest;
fxbest = funchandle(xbest);
end