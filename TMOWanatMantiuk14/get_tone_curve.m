function l_out = get_tone_curve( l_in, l_out )

delta = l_in(2)-l_in(1);

segs = length(l_in);

assert( length(l_in) == length(l_out) );


%E = tc_error( l_out )

% Non-equality constraints
A = cat( 2, eye( segs-1 ), zeros( segs-1,  1 ) ) + cat( 2, zeros( segs-1, 1 ), -eye( segs-1 ) );
b = zeros( segs-1, 1 );

Aeq = zeros( 2, segs );
Aeq(1,1) = 1;
Aeq(end,end) = 1;
beq = [ l_out(1) l_out(end) ]';

l_out = fmincon( @tc_error, l_out, A, b, Aeq, beq );

%E = tc_error( l_out )

    function E = tc_error( l_out )
        
        c = 0.1;
               
        err = c * (diff( l_out )/delta - 1)' + c_thr( l_in(1:(end-1)) ) - c_thr( l_out(1:(end-1)) );
        
        
        E = sum( err .^2 );
        
    end


end

function C = c_thr( l )

C = (1./(8*csf_hdrvdp( 2, 10.^l )))';
%C = michelson2log(1./(8*csf_hdrvdp( 5, 10.^l )))';

end