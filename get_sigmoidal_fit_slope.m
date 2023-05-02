%% Function to get the slope of the 

function [slope, xFit, yFit] = get_sigmoidal_fit_slope(x,y)

ft = fittype('a/(1+exp(-b*x))','independent','x','dependent','y');
opts = fitoptions('Method','NonlinearLeastSquares');
opts.Display = 'Off';
% opts.StartPoint = [-10 -1];

warning('off','curvefit:fit:noStartPoint')
[fitresult,~] = fit(x',y,ft,opts);

slope = fitresult.b;

F = @(x)(fitresult.a./(1+exp(-fitresult.b*x)));
xFit = linspace(min(x),max(x),1000);
yFit = F(xFit);

end