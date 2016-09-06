function [out] = covoper(vis, st, h, D, M, S)

Nx1=st.Nd(1);
Ny1=st.Nd(2);
Nx2=st.Kd(1);
Ny2=st.Kd(2);

vis2 = S'*vis;
vis2d = reshape(vis2, Ny2, Nx2);
spec = ifft2(full(vis2d));
specd = M*D.* spec(:);
spec_est1=M'*(h*sparse(specd)); % SLOW SLOW SLOW
spec_est2 = D.*spec_est1;
spec_est=reshape(spec_est2,Ny2,Nx2);
im21 = fft2(full(spec_est));
out = S*im21(:);