% derivative of the spherical hankel function
% function dhn=sphdhn(n,x)

function dhn=sphdhn(n,x)

dhn=sphbesldj(n,x)+i*sphbesldy(n,x);
