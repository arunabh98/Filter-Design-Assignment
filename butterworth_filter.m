%Butterworth Analog LPF parameters
Wc = 1.0375;              %cut-off frequency
N = 13;                  %order 

p1 = -1.00735 + 0.24829i;
p2 = -1.00735 - 0.24829i;
p3 = -0.918661 - 0.48215i;
p4 = -0.918661 + 0.48215i;
p5 = -0.77658 + 0.68799i;
p6 = -0.77658 - 0.68799i;
p7 = -0.589367 - 0.853846i;
p8 = -0.589367 + 0.853846i;
p9 = -0.367903 + 0.970079i;
p10 = -0.367903 - 0.970079i;
p11 = -0.125057 - 1.02994i;
p12 = -0.125057 + 1.02994i;
p13 = -1.0375;

%Band Edge speifications
fs1 = 28.4;
fp1 = 30.4;
fp2 = 40.4;
fs2 = 42.4;

%Transformed Band Edge specs using Bilinear Transformation
f_samp = 300;
ws1 = tan(fs1/f_samp*pi);          
wp1 = tan(fp1/f_samp*pi);
wp2 = tan(fp2/f_samp*pi);
ws2 = tan(fs2/f_samp*pi);

W0 = sqrt(wp1*wp2);
B = wp2-wp1; 

[num,den] = zp2tf([],[p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13],Wc^N);

%Evaluating Frequency Response of Final Filter
syms s z;
analog_lpf(s) = poly2sym(num,s)/poly2sym(den,s);    %analog lpf transfer function
analog_bpf(s) = analog_lpf((s*s +W0*W0)/(B*s));     %bandpass transformation
discrete_bpf(z) = analog_bpf((z-1)/(z+1));          %bilinear transformation
discrete_bpf(z) = vpa(simplify(vpa(expand(discrete_bpf(z)), 20)), 20);

%coeffs of analog BPF
[ns, ds] = numden(analog_bpf(s));                   %numerical simplification
ns = sym2poly(expand(ns));                          
ds = sym2poly(expand(ds));                          %collect coeffs into matrix form
k = ds(1);    
ds = ds/k;
ns = ns/k;

%coeffs of discrete BPF
[nz, dz] = numden(discrete_bpf(z));                 %numerical simplification
nz = sym2poly(expand(nz));                          
dz = sym2poly(expand(dz));                          %collect coeffs into matrix form
k = dz(1);                                          %normalisation factor
dz = dz/k;
nz = nz/k;
fvtool(nz,dz)                                       %frequency response in dB

%magnitude plot (not in log scale) 
[H,f] = freqz(nz,dz,1024*1024, 300e3);
plot(f,abs(H))
grid

                                                                                                           
