#!/usr/bin/python3

import numpy as np
import skrf as rf 
import sys, tempfile, os

def smatch(S):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    if rollet(S)[0] < 1: return np.nan, np.nan
    D = det(S)
    B1 = 1 + np.abs(S11)**2 - np.abs(S22)**2 - np.abs(D)**2;
    B2 = 1 + np.abs(S22)**2 - np.abs(S11)**2 - np.abs(D)**2;
    C1 = S11 - D * np.conj(S22)
    C2 = S22 - D * np.conj(S11)
    GS = (B1 - np.sign(B1) * np.sqrt(B1**2 - 4 * np.abs(C1)**2)) / (2 * C1)
    GL = (B2 - np.sign(B2) * np.sqrt(B2**2 - 4 * np.abs(C2)**2)) / (2 * C2)
    return GS, GL

def gin(S, GL):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    return S11 + S12 * S21 * GL / (1 - S22 * GL)

def gout(S, GS):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    return S22 + S12 * S21 * GS / (1 - S11 * GS)

def z2g(Z, Z0=50):
    return (Z - Z0) / (Z + Z0)

def g2z(G, Z0=50):
    return Z0 * (1 + G) / (1 - G)

def gmag(S):
    K, D = rollet(S)
    return gmsg(S) * (K - np.sqrt(K**2 - 1)) if K < np.inf else gum(S)

def gum(S):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    return np.abs(S21)**2 / ((1 - np.abs(S11)**2) * (1 - np.abs(S22)**2))

def gui(S):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    return 1 / (1 - np.abs(S11)**2)

def guo(S):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    return 1 / (1 - np.abs(S22)**2)

def gmsg(S):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    return np.abs(S21 / S12) if np.abs(S12) > 0 else np.inf

def gu(S):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    U = (S12 * S21 * np.conj(S11 * S22)) * gui(S) * guo(S)
    return 1 / np.abs(1 - U)**2

def det(S):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    return S11 * S22 - S12 * S21

def rollet(S):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    D = det(S)
    K = (1 - np.abs(S11)**2 - np.abs(S22)**2 + np.abs(D)**2) / \
        np.abs(2 * S12 * S21) if np.abs(S12 *S21) > 0 else np.inf
    return K, np.abs(D)

def mu(S):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    D = det(S)
    return (1 - np.abs(S11)**2) / (np.abs(S22 - D * np.conj(S11)) + np.abs(S12 * S21))

def s2abcd(S, Z0=50):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    den = 2 * S21
    A = ((1 + S11) *(1 - S22) + S12 * S21) / den
    B = Z0  * ((1 + S11) * (1 + S22) - S12 * S21) / den
    C = 1 / Z0 * ((1 - S11) * (1 - S22) - S12 * S21) / den
    D = ((1 - S11) * (1 + S22) + S12 * S21) / den
    return np.array([
        [ A, B ],
        [ C, D ]
    ])

def abcd2s(M, Z0=50):
    M = np.array(M)
    A = M[0,0]
    B = M[0,1]
    C = M[1,0]
    D = M[1,1]
    den = A + B / Z0 + C * Z0 + D
    S11 = (A + B / Z0 - C * Z0 - D) / den
    S12 = 2 * (A * D - B * C) / den
    S21 = 2 / den
    S22 = (-A + B / Z0 - C * Z0 + D) / den
    return np.array([
        [S11, S12],
        [S21, S22]
    ])

def tothreeport(S):
    S = np.array(S)
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    D11 = 1 - S11 - S12
    D12 = 1 - S11 - S21
    D21 = 1 - S12 - S22
    D22 = 1 - S21 - S22
    E = S11 + S22 + S12 + S21
    return np.array([
        [ S11+D11*D12/(4-E), S12+D11*D21/(4-E), 2*D11/(4-E) ],
        [ S21+D22*D12/(4-E), S22+D22*D21/(4-E), 2*D22/(4-E) ],
        [       2*D12/(4-E),       2*D21/(4-E),     E/(4-E) ]
    ])

def lift_ground(S, Z, Z0=50):
    S = tothreeport(S)
    G = z2g(Z, Z0)
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    S31 = S[2,0]
    S32 = S[2,1]
    S13 = S[0,2]
    S23 = S[1,2]
    S33 = S[2,2]
    return np.array([
       [ S11+S13*S31/(1/G-S33), S12+S13*S32/(1/G-S33) ],
       [ S21+S23*S31/(1/G-S33), S22+S23*S32/(1/G-S33) ]
    ])

def cbg_transform(S):    
    S = tothreeport(S)
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    S31 = S[2,0]
    S32 = S[2,1]
    S13 = S[0,2]
    S23 = S[1,2]
    S33 = S[2,2]
    return np.array([
       [ S33-S31*S13/(1+S11), S32-S31*S12/(1+S11) ],
       [ S23-S21*S13/(1+S11), S22-S21*S12/(1+S11) ]
    ])

def ccd_transform(S):
    S = tothreeport(S)
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    S31 = S[2,0]
    S32 = S[2,1]
    S13 = S[0,2]
    S23 = S[1,2]
    S33 = S[2,2]
    return np.array([
       [ S11-S12*S21/(1+S22), S13-S12*S23/(1+S22) ],
       [ S31-S32*S21/(1+S22), S33-S32*S23/(1+S22) ]
    ])

############################

def to_halfpi(rin, za):
    """
    RIN <---+---Z---< ZA
            Y
    RIN > ZA.real
    """
    ra, xa = za.real, za.imag
    xd = np.sqrt(ra * (rin - ra))
    if np.iscomplex(xd): raise ValueError
    x2 = np.array([-xa - xd, -xa + xd])
    x1 = -(ra**2 + (x2 + xa)**2) / (x2 + xa)
    return np.transpose([x1 * 1j, x2 * 1j]).tolist()

def to_halftee(rin, za):
    """
    RIN <---Z---+---<  ZA
                Y
    RIN < ZA.real
    """
    ra, xa = za.real, za.imag
    xd = np.sqrt(rin * ra * (ra**2 + xa**2 - rin * ra))
    if np.iscomplex(xd): raise ValueError
    x2 = np.array([(-rin * xa + xd) / (rin - ra),
                   (-rin * xa - xd) / (rin - ra)])
    x1 = -x2 * (ra**2 + xa * (x2 + xa)) / (ra**2 + (x2 + xa)**2)
    return np.transpose([x1 * 1j, x2 * 1j]).tolist()

def component(impedance, fd):
    w = 2 * np.pi * fd
    x = impedance.imag
    return 1 / (w * x) if x < 0 else x / w

##########################

def polar(x): 
    return '{:10.4g} {:7.2f}'.format(np.abs(x), np.angle(x) * 180 / np.pi)

def db(x):
    return 10 * np.log10(x) if x > 0 else np.nan

def read_network(path=None):
    if path is None:
        buf = sys.stdin.read()
        path = tempfile.mktemp() + '.s2p'
        with open(path, 'w') as f:
            f.write(buf)
        nw = rf.Network(path)
        os.unlink(path)
    else:
        nw = rf.Network(path)
    if np.any(nw.z0 != 50):
        print('Only networks referenced to 50 ohms supported', file=sys.stderr)
        sys.exit(1)
    return nw

def write_abcd(nw):
    print('MHZ             A                  B                  C                  D')
    for i in range(len(nw)):
        f = nw.f[i] / 1e6
        S = nw.s[i]
        buf = ' '.join([ polar(x) for x in s2abcd(S).flatten() ])
        print('{:<5g}'.format(f), buf)

def write_summary(nw):
    print('MHZ           ZIN             ZOUT '
          '        GUI    S21    GUO    GUM   GMSG   GMAG     GU'
          '        K       MU')
    for i in range(len(nw)):
        f = nw.f[i] / 1e6
        S = nw.s[i]
        S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
        K, D = rollet(S)
        GMAG = '     -' if K < 1 else '{:6.2f}'.format(db(gmag(S)))
        GU = '     -' if K < 1 else '{:6.2f}'.format(gu(S))
        print('{:<5g} {:16.4g} {:16.4g} {:6.2f} {:6.2f} {:6.2f} {:6.2f} '
            '{:6.2f} {:s} {:s} {:8.4g} {:8.4g}'.format(
            f, g2z(S11), g2z(S22), 
            db(gui(S)), db(np.abs(S21)**2), db(guo(S)), 
            db(gum(S)), db(gmsg(S)), GMAG, GU, K, mu(S)))

def write_sparam(nw):
    print('# MHZ S MA R 50')
    print('! MHZ           S11                S21                S12                S22      '
          '!   GUM        K       MU')
    for i in range(len(nw)):
        f = nw.f[i] / 1e6
        S = nw.s[i]
        K, D = rollet(S)
        buf = ' '.join([ polar(x) for x in S.T.flatten() ])
        print('{:<5g} {:s} ! {:5.1f} {:8.4g} {:8.4g}'.format(f, buf, db(gum(S)), K, mu(S)))

def write_match(nw, match, Z0=50):
    print('MHZ          ZS             ZL')
    for i in range(len(nw)):
        f = nw.f[i] / 1e6
        S = nw.s[i]
        K, D = rollet(S)
        ZS, ZL = (match + ',').split(',')[:2]
        if ZS and ZL:
            ZS, ZL = complex(ZS), complex(ZL)
        elif ZS:
            ZS = complex(ZS)
            ZL = g2z(gout(S, z2g(ZS)))
        elif ZL:
            ZL = complex(ZL)
            ZS = g2z(gin(S, z2g(ZL)))
        else:
            GS, GL = smatch(S)
            ZS, ZL = g2z(GS), g2z(GL)
       
        ZIN, ZOUT = np.conj(ZS), np.conj(ZL)
        if Z0 < np.abs(ZIN):
            res = to_halftee(Z0, ZIN)
            print('ht ', end='') 
        else:
            res = to_halfpi(Z0, ZIN)
            print('hp ', end='') 
        print(component(res[0][0], 2e9), component(res[0][1], 2e9))
        print(component(res[1][0], 2e9), component(res[1][1], 2e9))
        if Z0 < np.abs(ZOUT):
            res = to_halftee(Z0, ZOUT) 
            print('ht ', end='') 
        else:
            res = to_halfpi(Z0, ZOUT)
            print('hp ', end='') 
        print(component(res[0][0], 2e9), component(res[0][1], 2e9))
        print(component(res[1][0], 2e9), component(res[1][1], 2e9))
        print('{:<5g} {:14.4g} {:14.4g}'.format(f, ZS, ZL))

def write_network(nw, mode, match):
    if mode == 'a': 
        write_abcd(nw)
    elif mode == 's':
        write_summary(nw)
    elif mode == 'n':
        print(nw.write_touchstone(form='ma', return_string=True))
    elif mode == 'm':
        write_match(nw, match)
    else:
        write_sparam(nw)


def main(*args):
    args = list(args)
    mode = None
    match = None
    stack = []
    stack.append(read_network())

    while args:
        opt = args.pop(0)
        top = stack[-1]

        if opt == '-n':
            mode = 'n'
        elif opt == '-a':
            mode = 'a'
        elif opt == '-s':
            mode = 's'
        elif opt == '-match':
            match = args.pop(0).strip()
            mode = 'm'

        elif opt == '-swap':
            b = stack.pop()
            a = stack.pop()
            stack.append(b)
            stack.append(a)
        elif opt == '-cascade':
            b = stack.pop()
            a = stack.pop()
            stack.append(a ** b)
        elif opt == '-deembed':
            b = stack.pop()
            a = stack.pop()
            stack.append(a ** b.inv)
        elif opt == '-ideembed':
            b = stack.pop()
            a = stack.pop()
            stack.append(a.inv ** b)

        elif opt == '-p':
            write_output(top, mode=mode)
        elif opt == '-f':
            stack.append(read_network(args.pop(0)))
        elif top and opt == '-series':
            S = abcd2s([[1, float(args.pop(0))], [0, 1]])
            stack.append(rf.Network(frequency=top.frequency, s=[S] * len(top)))
        elif top and opt == '-shunt':
            S = abcd2s([[1, 0], [1/float(args.pop(0)), 1]])
            stack.append(rf.Network(frequency=top.frequency, s=[S] * len(top)))

        elif top and opt == '-cbg':
            top.s = np.array([ cbg_transform(S) for S in top.s ])
        elif top and opt == '-ccd':
            top.s = np.array([ ccd_transform(S) for S in top.s ])
        elif top and opt == '-lift':
            x = args.pop(0)
            top.s = np.array([ lift_ground(
                top.s[i], 
                complex(x) if 'j' in x else 2j * np.pi * top.f[i] * float(x)
            ) for i in range(len(top)) ])
        else:
            print('Unrecognized command line option. Exiting.', file=sys.stderr)
            sys.exit(1)

    if stack: 
        write_network(stack[-1], mode=mode, match=match)


if __name__ == '__main__':
    np.seterr(divide='ignore')
    main(*sys.argv[1:])


