% stub2.m - double-stub matching
%                           
%  -----------------/-----------/---|
%  main line Z0    /           /    ZL
%  ---------------/---/---l---/---/-|
%                /   d1      /   d2
%               /___/       /___/
%
% Usage: d12 = stub2(zL,l,type)
%        d12 = stub2(zL,l)      (equivalent to type='ss')
%        d12 = stub2(zL)        (equivalent to l=1/8 and type='ss')
%
% zL = normalized load impedance, i.e., zL = ZL/Z0
% l  = fixed separation of stubs in wavelengths, typically, l=1/8
% type = 'ss','so','os','oo' for short/short, short/open, open/short, open/open
%
% d12 = [d1,d2] = 2x2 matrix, where each row is a solution, 
%
% d1 is length of stub-2 located at distance l from load
% d2 is length of stub-1 located at load 
%
% notes: d1,d2 are in wavelengths and are reduced mod lambda/2
%
%        requires that gL <= gmax, where
%        yL = 1/zL = gL + j*bL, 
%        gmax = 1 + cot(kl)^2 = 1/sin(kl)^2
%        if not, use STUB3
%   fprintf('stub separation must be less than lmax = %.4f(lambda)\n\n', lmax);
% Sophocles J. Orfanidis - 1999-2008 - www.ece.rutgers.edu/~orfanidi/ewa

def stub2(zL, l=45, mode='ss'):
    yL = 1 / zL 
    gL, bL = yL.real, bL = yL.imag

    m = (type=='o');                    % selects open stubs    

    c = 1 / np.cot(np.deg2rad(l))
    gmax = 1 + c**2
    if gL > gmax: gL = np.nan

    b = c + np.array([1, -1]) * np.sqrt(gL * (gmax - gL))
    d2 = acot(bL - b) + m(2) * np.pi / 2
    d1 = acot((c-b-gL*c)/gL) + m(1) * np.pi / 2
    d12 = np.mod([d1,d2], np.pi)
    return d12



