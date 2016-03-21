function xd = getWaveletDenoisedTrace(x)
deb = x(1);
scal = 'mln'; % one='while noise', sln='unscaled noise', mln='non-white noise'; Use a level-dependent estimation of the level noise
xd = wden(x - deb,'sqtwolog','s',scal,3,'db3') + deb;

end