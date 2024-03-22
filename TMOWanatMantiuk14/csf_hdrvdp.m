function S = csf_hdrvdp( rho, lum )
% Contrast sensitivity function
% S = hdrvdp_csf( rho, lum )

csf_pars = [ ...
   0.0160737   0.991265   3.74038   0.50722   4.46044
   0.383873   0.800889   3.54104   0.682505   4.94958
   0.929301   0.476505   4.37453   0.750315   5.28678
   1.29776   0.405782   4.40602   0.935314   5.61425
   1.49222   0.334278   3.79542   1.07327   6.4635
   1.46213   0.394533   2.7755   1.16577   7.45665 ];

lum_lut = log10( [ 0.002 0.02 0.2 2 20 150] );
log_lum = log10( lum );

metric_par.csf_sa = [30.162 4.0627 1.6596 0.2712];

metric_par.csf_sr_par = [1.1732 1.1478 1.2167 0.5547 2.9899 1.1414]; % rod sensitivity function

par = [0.061466549455263 0.99727370023777070]; % old parametrization of MTF
metric_par.mtf_params_a = [par(2)*0.426 par(2)*0.574 (1-par(2))*par(1) (1-par(2))*(1-par(1))];
metric_par.mtf_params_b = [0.028 0.37 37 360];



par = cell( 4, 1 ); %zeros(length(lum),4);
for k=1:4
    par{k} = interp1( lum_lut, csf_pars(:,k+1), clamp( log_lum, lum_lut(1), lum_lut(end) ) );
end

S =  (par{4} .* 1./((1+(par{1}.*rho).^par{2}) .* 1./(1-exp(-(rho/7).^2)).^par{3}).^0.5) .* ...
    hdrvdp_mtf( rho, metric_par ) .* hdrvdp_joint_rod_cone_sens( lum, metric_par );

end

function MTF = hdrvdp_mtf( rho, metric_par )
% Custom-fit MTF of the eye

MTF = zeros(size(rho));

for kk=1:4
    MTF = MTF + metric_par.mtf_params_a(kk) * exp( -metric_par.mtf_params_b(kk) * rho );
end

end

function S = hdrvdp_joint_rod_cone_sens( la, metric_par )


cvi_sens_drop = metric_par.csf_sa(2); % in the paper - p6
cvi_trans_slope = metric_par.csf_sa(3); % in the paper - p7
cvi_low_slope = metric_par.csf_sa(4); % in the paper - p8


S = metric_par.csf_sa(1) * ( (cvi_sens_drop./la).^cvi_trans_slope+1).^-cvi_low_slope;

end