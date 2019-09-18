

function [z, reps, iter, lmb, kvec, R, info] = ptp(X, maxit, eps, verbose, kvec0, R0)

if verbose >= 0
        printf("\n%s\n\n", "ptp.w,v 1.6 2018/02/03 04:18:31 nurmi Exp");
        printf(" dim = %4d nvec = %5d maxit = %5d eps = %12.4e\n\n", size(X), maxit, eps);
endif


global OPTIMAL;
global OVERSIZED_BASIS;
global LINDEP_BASIS;
global INIT_ERROR;
global NO_WAY_TOI;
global CHOLIN_CRASH;
global BAD_OLD_L;
global MULTI_IN;
global DBLE_VTX;
global EXCESS_ITER;
global RUNNING;
global NEGGRAMM;


OPTIMAL                 =   0;
OVERSIZED_BASIS         = -19;
LINDEP_BASIS            = -20;
INIT_ERROR              = -21;
NO_WAY_TOI              = -22;
CHOLIN_CRASH            = -23;
BAD_OLD_L               = -24;
MULTI_IN                = -25;
DBLE_VTX                = -26;
EXCESS_ITER             = -27;
RUNNING                 = -28;
NEGGRAMM                = -29;
info = RUNNING;

epstol = eps;

[ xdim nvec ] = size(X);

if columns(kvec0) == 0
        
        [ vmin ifac ] = min(sumsq(X));
        curr.kvec = [ ifac ];
        curr.lmb = [ 1 ];
        curr.R = norm(X(:, ifac));
        
else
         
        if rows(R0) != columns(kvec0)
                info = INIT_ERROR;
                printf("XXXX INIT_ERROR: nonmatching sizes of kvec0 (%d %d)\n", size(kvec0));
                printf(" and R0 (%d %d).\n", size(R0));
                printf(" Reverting to the cold start.\n");
                
                [ vmin ifac ] = min(sumsq(X));
                curr.kvec = [ ifac ];
                curr.lmb = [ 1 ];
                curr.R = norm(X(:, ifac));
                
        else
                curr.kvec = kvec0;
                curr.R = R0;
                lmb = R0 \ ( R0' \ ones(size(kvec0')));
                curr.lmb = lmb/sum(lmb);
        endif
        
endif

curr.z = X(:, curr.kvec)*curr.lmb;
 
report.zz = zz2(curr);
report.zx = 0;
report.in = ifac;
report.out = 0;
report.iter = zeros(1, 2);
report.lenbas = columns(curr.kvec);
(verbose >= 0) && report_iter(report);
 
[ vmin ifac ] = get_ifac(X, curr, epstol);
if ifac == 0
        
        z = curr.z;
        reps = eps;
        iter = report.iter;
        lmb = curr.lmb;
        kvec = curr.kvec;
        R = curr.R;
        info = OPTIMAL;
        return
        
        return
endif
[curr.kvec curr.lmb curr.R del_iter] = newbas(curr.kvec, curr.lmb, ifac, curr.R, X);
curr.z = X(:, curr.kvec)*curr.lmb;

report.iter(1)++;
report.iter(2) += del_iter;
report.zz = zz2(curr); report.zx = vmin;
report.in = ifac; report.lenbas = columns(curr.kvec);
(verbose > 0) && (mod(report.iter(1), verbose) == 0) && report_iter(report);

while ifac && ( report.iter <= maxit ) 
        [ vmin ifac ] = get_ifac(X, curr, epstol);
        if ifac == 0
                info = OPTIMAL;
                break
        endif
        [curr.kvec curr.lmb curr.R del_iter] = newbas(curr.kvec, curr.lmb, ifac, curr.R, X);
        curr.z = X(:, curr.kvec)*curr.lmb;
        
        report.iter(1)++;
        report.iter(2) += del_iter;
        report.zz = zz2(curr); report.zx = vmin;
        report.in = ifac; report.lenbas = columns(curr.kvec);
        (verbose > 0) && (mod(report.iter(1), verbose) == 0) && report_iter(report);
        
endwhile
report.in = ifac;
report.zz = zz2(curr); report.zx = vmin;
if (verbose >= 0)
        if info == OPTIMAL
                printf("\nXXXX Solved to optimality\n\n");
        endif
        report_iter(report);
endif


z       = curr.z;
reps    = eps;
lmb     = curr.lmb;
kvec    = curr.kvec;
R       = curr.R;
iter    = report.iter;

endfunction


function [ lambda ] = proplus(kvec, ifac, R, X)

g = X(:, ifac);
r = (g'*X(:, kvec))';
if columns(kvec) == rows(X)
        
        lambda = -R\(R'\r);
        xi = 1/(1 + sum(lambda));
        lambda = xi*[ lambda; 1 ];
        
        return
endif
Re = [ r ones(size(r)) ];
Z = R\(R'\Re);
A = [ [ sumsq(g) 1 ]; [ 1 0 ] ] - Re'*Z;

xit = (inv(A))(:, 2);
lambda = [ -Z*xit; xit(1) ];
endfunction


function [kvec_new lambda_new R_new iter ] = newbas(kvec, lambda_old, ifac, R, X)
iter = iout = 0;
[ lambda_new ] = proplus(kvec, ifac, R, X);
kvec_new = kvec;
if all(lambda_new >= -eps)
        
        kvec_new = [ kvec_new ifac ];
        R_new = lastadd(X(:, kvec_new), R);
        return
        
endif
[ lambda_m izero ] = mid_lambda([ lambda_old; 0 ], lambda_new);
if izero == 0
        printf(" ST-OPT !!!\n");
        exit(-1)
endif
iout = kvec_new(izero);

kvec_new = [ kvec_new(1:izero-1) kvec_new(izero+1: end)];
lambda = [ lambda_m(1:izero-1); lambda_m(izero+1: end) ];
R_new = choldelete(R, izero);
kvec_new = [ kvec_new ifac ];
R_new = lastadd(X(:, kvec_new), R_new);
lambda_new = baric(kvec_new, R_new);

while any(lambda_new < -eps)
        
        [ lambda_m izero ] = mid_lambda( lambda, lambda_new);
        kvec_new = [ kvec_new(1:izero-1) kvec_new(izero+1: end) ];
        lambda = [ lambda_m(1:izero-1); lambda_m(izero+1: end) ];
        R_new = choldelete(R_new, izero);
        lambda_new = baric(kvec_new, R_new);
        z = X(:, kvec_new)*lambda_new;
        
        iter++;
endwhile
endfunction


function [ vmin ifac ] = get_ifac(X, obj, epstol)
[ vmin ifac ] = min(obj.z'*X - sumsq(obj.z));
reps = epstol*norm(X(:, ifac));
( vmin > -reps ) && ( ifac = 0 );
endfunction


function RU = lastadd(X, R)
u = X(:, end)'*X;
q = R'\u(1:end-1)';
zz = sqrt(abs(u(end) - sumsq(q)));
RU = [ R; zeros(1, columns(R)) ];
RU = [ RU [ q; zz ] ];
endfunction


function [ zz ] = zz2(aim)
zz = sumsq(aim.z);
endfunction


function report_iter(report)
printf(" ++");
printf(" iter %4d(+) %4d(-)", report.iter);
printf(" len %4d", report.lenbas);
printf(" zx %12.4e", report.zx);
printf(" in %6d", report.in);
printf(" zz %16.8e", report.zz);
printf("\n");
endfunction




function [ lambda ] = baric(kvec, R)
lambda = R\(R'\ones(rows(R), 1));
lambda = lambda ./ sum(lambda);
endfunction


function [ lambda_m izero ] = mid_lambda(lambda_old, lambda_new)
if all(lambda_new < 0)
        lambda_m = lambda_new;
        izero = 0;
        return;
endif
[ theta imin ] = min(([ lambda_old ] ./ ([ lambda_old ] - lambda_new))( lambda_new < 0 ));
izero = [ 1:rows(lambda_new) ]( lambda_new < 0 )(imin);
lambda_m = theta*lambda_new + (1 - theta)*lambda_old;
endfunction

