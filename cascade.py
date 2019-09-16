#!/usr/bin/python3

import numpy as np
import skrf as rf 
import sys, tempfile, os

def z2g(Z, Z0=50):
    return (Z - Z0) / (Z + Z0)

def g2z(G, Z0=50):
    return Z0 * (1 + G) / (1 - G)

### matching

def to_stub1(za, zo=50, shorted=True): # match with a stub-series input 
    """
    -----------------/-----------|
    main line zo    /            za
    ---------------/---/----l----|
                  /   d
                 /___/
    """
    GL = z2g(za, zo)
    thL = np.angle(GL)
    bl = thL / 2 + np.array([1, -1]) * np.arccos(-abs(GL)) / 2
    if shorted:
        bd = np.arctan(-np.tan(2 * bl - thL) / 2)
    else:
        bd = np.arctan(1 / (np.tan(2 * bl - thL) / 2))
    d = np.mod([ bd, bl ], np.pi)
    d = np.rad2deg(d)
    return d.T

def to_qwt2(za, zo=50, shorted=True):
    """
    ---------------==========----|--|
    main line zo       z1        |  za
    ---------------==========-|--|--|
                     L1=1/4   |  |
                              |z2| L2=1/8 shorted or
                              |__|    3/8 opened
    """
    ya = 1 / za
    gl, bl = ya.real, ya.imag
    z1 = np.sqrt(zo / gl) * np.array([1, 1])
    z2 = 1 / bl * np.array([1, -1]) * (1 if shorted else -1)
    d = np.array([ 45, 135 ])
    res = np.array([ z1, z2, d ]).T
    return res[1] if z2[0] < 0 else res[0]

def to_qwt3(za, z2, zo=50, shorted=True):
    """
    ---------------==========----|--|
    main line zo       z1        |  za
    ---------------==========-|--|--|
                     L1=1/4   |  |
                              |z2| d=[shorted,opened]
                              |__|
    """
    ya = 1 / za
    gl, bl = ya.real, ya.imag
    z1 = np.nan if zo / gl < 0 else np.sqrt(zo / gl) 
    d = np.arctan([ 1 / (bl * z2), -bl * z2 ])
    d = np.mod(d, np.pi)
    d = np.rad2deg(d)
    return np.array([
        [ z1, d[0] ],
        [ z1, d[1] ]
    ])[0 if shorted else 1]

def lmatch(ZS, ZL, reverse=False):
    """
    ZS <---+---X2--< ZL
           X1   
    """
    if reverse: ZS, ZL = ZL, ZS
    RS, XS = ZS.real, ZS.imag
    RL, XL = ZL.real, ZL.imag
    QS = RS / RL - 1 + XS**2 / (RS * RL)
    if QS < 0: return [[ np.nan, np.nan], [ np.nan, np.nan ]]
    Q = np.sqrt(QS)
    X1 = (XS + np.array([1, -1]) * Q * RS) / (RS / RL - 1)
    X2 = -(XL + np.array([1, -1]) * Q * RL)
    return np.array([X1, X2]).T

### scattering parameters

def smatch(S):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    K = rollet(S)
    if K < 1: return np.nan, np.nan
    D = det(S)
    B1 = 1 + np.abs(S11)**2 - np.abs(S22)**2 - np.abs(D)**2;
    B2 = 1 + np.abs(S22)**2 - np.abs(S11)**2 - np.abs(D)**2;
    C1 = S11 - D * np.conj(S22)
    C2 = S22 - D * np.conj(S11)
    den = B1**2 - 4 * np.abs(C1)**2
    GS = np.nan if den < 0 or C1 == 0 else (B1 - np.sign(B1) * np.sqrt(den)) / (2 * C1)
    den = B2**2 - 4 * np.abs(C2)**2
    GL = np.nan if den < 0 or C2 == 0 else (B2 - np.sign(B2) * np.sqrt(den)) / (2 * C2)
    return GS, GL

def gin(S, GL):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    return S11 + S12 * S21 * GL / (1 - S22 * GL)

def gout(S, GS):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    return S22 + S12 * S21 * GS / (1 - S11 * GS)

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
    return np.nan if S12 == 0 else np.abs(S21) / np.abs(S12)

def gmag(S):
    K = rollet(S)
    if np.isinf(K): return gum(S)
    return np.nan if K < 1 else gmsg(S) * (K - np.sqrt(K**2 - 1))

def gu(S):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    U = (S12 * S21 * np.conj(S11 * S22)) * gui(S) * guo(S)
    K = rollet(S)
    return np.nan if K < 1 else 1 / np.abs(1 - U)**2

def det(S):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    return S11 * S22 - S12 * S21

def rollet(S):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    D = det(S)
    K = (1 - np.abs(S11)**2 - np.abs(S22)**2 + np.abs(D)**2) / np.abs(2 * S12 * S21) 
    return K

def mu(S):
    S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
    D = det(S)
    return (1 - np.abs(S11)**2) / (np.abs(S22 - D * np.conj(S11)) + np.abs(S12 * S21))

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

### abcd

def tline(deg, zo=50):
    theta = np.deg2rad(deg)
    return np.array([
        [ np.cos(theta), 1j * zo * np.sin(theta) ],
        [ 1j * np.sin(theta) / zo, np.cos(theta) ]
    ])

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

### input/output

def open_stub(deg, zo=50):
    theta = np.deg2rad(deg)
    return -1j * zo / np.tan(theta)

def shorted_stub(deg, zo=50):
    theta = np.deg2rad(deg)
    return 1j * zo * np.tan(theta)

def to_complex(s):
    if '/' in s:
        r, theta = s.split('/')
        return float(r) * np.exp(1j * float(theta) * np.pi / 180)
    else:
        return complex(s)

def matching(S, GS, GL):
    if GS is None and GL is None:
        GS, GL = smatch(S)
    elif GL is None:
        GL = np.conj(gout(S, GS))
    elif GS is None:
        GS = np.conj(gin(S, GL))
    GIN, GOUT = gin(S, GL), gout(S, GS)
    return g2z(GS), g2z(GL), g2z(GIN), g2z(GOUT)

def notation(x, precision=5):
    UNITS = 'FH'
    SUFFIX = ["f", "p", "n", "u", "m", "", "k", "M", "G"]
    exp = np.floor(np.log10(np.absolute(x)))
    mant = round(x / 10**exp, precision-1)
    p = int(exp // 3)
    value = (mant * 10**exp) / 10**(3 * p)
    return "%6.4g%s%s" % (np.absolute(value), SUFFIX[p-4], UNITS[0 if x < 0 else 1])

def fm(mode, *d, f=None):
    res = []
    for m, x in zip(list(mode), d):
        if m == 'p':
            res.append('{:>10s}'.format('-') if np.isnan(x) else '{:10.4g}'.format(np.abs(x)))
            res.append('{:>7s}'.format('-') if np.isnan(x) else '{:7.2f}'.format(np.angle(x) * 180 / np.pi))
        if m == 'c':
            res.append('{:>16s}'.format('-') if np.isnan(x) else '{:16.4g}'.format(x))
        if m == 'd':
            x = 10 * np.log10(x) if x > 0 else np.nan
            res.append('{:>6s}'.format('-') if np.isnan(x) else '{:6.2f}'.format(x))
        if m == 'f':
            res.append('{:>6s}'.format('-') if np.isnan(x) else '{:6.2f}'.format(x))
        if m == 'a':
            res.append('{:>7s}'.format('-') if np.isnan(x) else '{:7.2f}'.format(x))
        if m == 'g': 
            res.append('{:>8s}'.format('-') if np.isnan(x) else '{:8.4g}'.format(x))
        if m == 'x': 
            w = 2 * np.pi * f
            x = 1 / (w * x) if x < 0 else x / w
            res.append('{:>8s}'.format('-') if np.isnan(x) or np.isinf(x) else notation(x))
        if m == 'F': 
            res.append('{:<5g}'.format(x))
    return ' '.join(res)

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
        f = nw.f[i]
        S = nw.s[i]
        print(fm('F', f / 1e6), fm('pppp', *s2abcd(S).flatten()))

def write_sparam(nw):
    print('# MHZ S MA R 50')
    print('! MHZ           S11                S21                '
          'S12                S22      !    GUM        K       MU')
    for i in range(len(nw)):
        f = nw.f[i]
        S = nw.s[i]
        K = rollet(S)
        print(fm('F', f / 1e6), fm('pppp', *S.T.flatten()), '!', fm('dgg', gum(S), K, mu(S)))

def write_summary(nw):
    print('MHZ           Z11              Z22         GUI    S21    GUO    '
          'GUM   GMSG   GMAG     GU        K        D       MU')
    for i in range(len(nw)):
        f = nw.f[i]
        S = nw.s[i]
        S11, S12, S21, S22 = S[0,0], S[0,1], S[1,0], S[1,1]
        K = rollet(S)
        print(fm('F', f / 1e6), fm('ccddddddfggg', g2z(S11), g2z(S22), gui(S), 
            np.abs(S21)**2, guo(S), gum(S), gmsg(S), gmag(S), gu(S), K, 
            np.abs(det(S)), mu(S)))

def write_lmatch(nw, data):
    print('MHZ      SHUNT   SERIES !   SERIES    SHUNT          '
          'ZS               ZL         SHUNT   SERIES !   SERIES    SHUNT')
    for i in range(len(nw)):
        f = nw.f[i]
        ZS, ZL, ZIN, ZOUT = matching(nw.s[i], data.get('gs'), data.get('gl'))
        for i in range(2):
            print(fm('F', f / 1e6), 
                fm('xx', *lmatch(50, np.conj(ZS))[i], f=f), '!',
                fm('xx', *lmatch(50, np.conj(ZS), 'r')[i], f=f), 
                fm('cc', ZS, ZL),
                fm('xx', *lmatch(np.conj(ZL), 50)[i], f=f), '!',
                fm('xx', *lmatch(np.conj(ZL), 50, 'r')[i], f=f))

def write_stub1(nw, data):
    print('MHZ    LSHUNT LSERIES          ZS               ZL      LSERIES  LSHUNT')
    for i in range(len(nw)):
        f = nw.f[i]
        ZS, ZL, ZIN, ZOUT = matching(nw.s[i], data.get('gs'), data.get('gl'))
        for i in range(2):
            for shorted in [ False, True ]:
                print(fm('F', f / 1e6),
                    fm('aa', *to_stub1(np.conj(ZS), shorted=shorted)[i]),
                    fm('cc', ZS, ZL),
                    fm('aa', *to_stub1(np.conj(ZL), shorted=shorted)[i][::-1]),
                    'shorted' if shorted else 'open')

def write_qwt2(nw, data):
    print('MHZ       ZQWT   ZSHUNT  LSHUNT          ZS               '
          'ZL       LSHUNT   ZSHUNT     ZQWT')
    for i in range(len(nw)):
        f = nw.f[i]
        ZS, ZL, ZIN, ZOUT = matching(nw.s[i], data.get('gs'), data.get('gl'))
        for shorted in [ False, True ]:
            print(fm('F', f / 1e6),
                fm('gga', *to_qwt2(np.conj(ZS), shorted=shorted)),
                fm('cc', ZS, ZL),
                fm('agg', *to_qwt2(np.conj(ZL), shorted=shorted)[::-1]),
                'shorted' if shorted else 'open')

def write_qwt3(nw, data):
    print('MHZ       ZQWT  LSHUNT          ZS               ZL       LSHUNT     ZQWT')
    for i in range(len(nw)):
        f = nw.f[i]
        z2 = data.get('z2')
        ZS, ZL, ZIN, ZOUT = matching(nw.s[i], data.get('gs'), data.get('gl'))
        for shorted in [ False, True ]:
            print(fm('F', f / 1e6), 
                fm('ga', *to_qwt3(np.conj(ZS), z2, shorted=shorted)),
                fm('cc', ZS, ZIN, ZOUT, ZL),
                fm('ag', *to_qwt3(np.conj(ZL), z2, shorted=shorted)[::-1]),
                'shorted' if shorted else 'open')

def write_match(nw, data):
    print('MHZ            ZS              ZIN             ZOUT               ZL')
    for i in range(len(nw)):
        f = nw.f[i]
        ZS, ZL, ZIN, ZOUT = matching(nw.s[i], data.get('gs'), data.get('gl'))
        print(fm('Fcccc', f / 1e6, ZS, ZIN, ZOUT, ZL))

def write_gamma(nw, data):
    print('MHZ             GS                GIN                GOUT                GL')
    for i in range(len(nw)):
        f = nw.f[i]
        ZS, ZL, ZIN, ZOUT = matching(nw.s[i], data.get('gs'), data.get('gl'))
        print(fm('Fpppp', f / 1e6, z2g(ZS), z2g(ZIN), z2g(ZOUT), z2g(ZL)))

def write_network(nw, data):
    mode = data.get('mode')
    if mode == 'a': 
        write_abcd(nw)
    elif mode == 's':
        write_summary(nw)
    elif mode == 'z':
        write_match(nw, data)
    elif mode == 'g':
        write_gamma(nw, data)
    elif mode == 'lmatch':
        write_lmatch(nw, data)
    elif mode == 'stub':
        write_stub1(nw, data)
    elif mode == 'qwt2':
        write_qwt2(nw, data)
    elif mode == 'qwt3':
        write_qwt3(nw, data)
    else:
        write_sparam(nw)

def main(*args):
    args = list(args)
    data = {}
    stack = []
    a = read_network()

    # make the first network on the stack a zero network

    b = a.copy()
    for S in b.s:
        S[0,0], S[0,1], S[1,0], S[1,1] = 0, 0, 0, 0
    stack.append(b)
    stack.append(a)

    while args:
        opt = args.pop(0)
        top = stack[-1]

        # diplay options

        if opt == '-a':
            data['mode'] = 'a'
        elif opt == '-s':
            data['mode'] = 's'
        elif opt == '-z':
            data['mode'] = 'z'
        elif opt == '-g':
            data['mode'] = 'g'
        elif opt == '-stub':
            data['mode'] = 'stub'
        elif opt == '-lmatch':
            data['mode'] = 'lmatch'
        elif opt == '-qwt':
            data['mode'] = 'qwt2'
        elif opt == '-qwtz':
            data['z2'] = float(args.pop(0))
            data['mode'] = 'qwt3'

        # matching options

        elif opt == '-zs':
            data['gs'] = z2g(to_complex(args.pop(0)))
        elif opt == '-zl':
            data['gl'] = z2g(to_complex(args.pop(0)))
        elif opt == '-gs':
            data['gs'] = to_complex(args.pop(0))
        elif opt == '-gl':
            data['gl'] = to_complex(args.pop(0))
        elif opt == '-zin':
            data['gin'] = z2g(to_complex(args.pop(0)))
        elif opt == '-zout':
            data['gout'] = z2g(to_complex(args.pop(0)))
        elif opt == '-gin':
            data['gin'] = to_complex(args.pop(0))
        elif opt == '-gout':
            data['gout'] = to_complex(args.pop(0))

        # binary operations

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

        # unary operations

        elif opt == '-unilateral':
            GS, GL = data.pop('gs', None), data.pop('gl', None)
            GIN, GOUT = data.pop('gin', None), data.pop('gout', None)
            for S in top.s:
                ZS, ZL, ZIN, ZOUT = matching(S, GS, GL)
                S[0,0] = z2g(ZIN)if GIN is None else GIN 
                S[0,1] = 0
                S[1,1] = z2g(ZOUT)if GOUT is None else GOUT
        elif opt == '-cbg':
            top.s = np.array([ cbg_transform(S) for S in top.s ])
        elif opt == '-ccd':
            top.s = np.array([ ccd_transform(S) for S in top.s ])
        elif opt == '-lift':
            x = args.pop(0)
            top.s = np.array([ lift_ground(
                top.s[i], 
                complex(x) if 'j' in x or '/' in x else 2j * np.pi * top.f[i] * float(x)
            ) for i in range(len(top)) ])

        # push operations

        elif opt == '-f':
            stack.append(read_network(args.pop(0)))
        elif opt == '-tline':
            x = to_complex(args.pop(0))
            S = abcd2s(tline(np.angle(x) * 180 / np.pi, zo=np.abs(x)))
            stack.append(rf.Network(frequency=top.frequency, s=[S] * len(top)))
        elif opt == '-series':
            x = to_complex(args.pop(0))
            S = abcd2s([[1, x], [0, 1]])
            stack.append(rf.Network(frequency=top.frequency, s=[S] * len(top)))
        elif opt == '-shunt':
            x = to_complex(args.pop(0))
            S = abcd2s([[1, 0], [1/x, 1]])
            stack.append(rf.Network(frequency=top.frequency, s=[S] * len(top)))
        elif opt == '-open':
            x = to_complex(args.pop(0))
            x = open_stub(np.angle(x) * 180 / np.pi, zo=np.abs(x))
            S = abcd2s([[1, 0], [1/x, 1]])
            stack.append(rf.Network(frequency=top.frequency, s=[S] * len(top)))
        elif opt == '-shorted':
            x = to_complex(args.pop(0))
            x = shorted_stub(np.angle(x) * 180 / np.pi, zo=np.abs(x))
            S = abcd2s([[1, 0], [1/x, 1]])
            stack.append(rf.Network(frequency=top.frequency, s=[S] * len(top)))

        # unrecognized

        else:
            print('Unrecognized command line option. Exiting.', file=sys.stderr)
            sys.exit(1)

    if stack: 
        write_network(stack[-1], data)


if __name__ == '__main__':
    np.seterr(divide='ignore')
    main(*sys.argv[1:])


