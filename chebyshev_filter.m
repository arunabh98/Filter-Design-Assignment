%Chebyshev LPF parameters
D1 = 1/(0.85*0.85)-1;       %since delta is 0.15
epsilon = sqrt(D1);         %epsilon was set to this value to satisfy required inequality
N = 3;

% Open CLHP Poles of the Chebyshev Polynomial of order 4
p1 = -0.12216 - 0.96981i;
p2 = -0.12216 + 0.96981i;
p3 = -0.29492 + 0.40171i;
p4 = -0.29492 - 0.40171i;

%evaluating the Transfer function of Chebyshev Analog LPF
n1 = [1 -p1-p2 p1*p2];
n2 = [1 -p3-p4 p3*p4];
den = conv(n1,n2);          %multiply n1 and n2, which are the two quadratic factors in the denominator
num = [den(5)*sqrt(1/(1+epsilon*epsilon))];         % even order, DC Gain set as 1/(1+ epsilon^2)^0.5

%Band Edge speifications
fs1 = 19.2;
fp1 = 21.2;
fp2 = 27.2;
fs2 = 29.2;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 200;
wp1 = tan(fs1/f_samp*pi);          
ws1 = tan(fp1/f_samp*pi);
ws2 = tan(fp2/f_samp*pi);
wp2 = tan(fs2/f_samp*pi);

%Parameters for Bandstop Transformation
W0 = sqrt(wp1*wp2);
B = wp2-wp1;

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bsf(s) = analog_lpf((B*s)/(s*s + W0*W0));        %bandstop transformation
discrete_bsf(z) = analog_bsf((z-1)/(z+1));              %bilinear transformation

%coeffs of analog bsf
[ns, ds] = numden(analog_bsf(s));                   %numerical simplification to collect coeffs
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete bsf
[nz, dz] = numden(discrete_bsf(z));                     %numerical simplification to collect coeffs                    
nz = sym2poly(expand(nz));
dz = sym2poly(expand(dz));                              %coeffs to matrix form
k = dz(1);                                              %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                           %frequency response

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 200e3);
plot(f,abs(H))
grid
