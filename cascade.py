#!/usr/bin/python3

import numpy as np
import skrf as rf 
import sys, tempfile, os


def db(x):
    return 10 * np.log10(x) if x > 0 else np.nan

def z2g(Z, Z0=50):
    return (Z - Z0) / (Z + Z0)

def g2z(G, Z0=50):
    return Z0 * (1 + G) / (1 - G)

def gu(S):
    S11 = S[0,0]
    S22 = S[1,1]
    S21 = S[1,0]
    return np.abs(S21)**2 / ((1 - np.abs(S11)**2) * (1 - np.abs(S22)**2))

def gui(S):
    S11 = S[0,0]
    return 1 / (1 - np.abs(S11)**2)

def guo(S):
    S22 = S[1,1]
    return 1 / (1 - np.abs(S22)**2)

def gmsg(S): # maximum stable gain, use with K<1
    S21 = S[1,0]
    S12 = S[0,1]
    return np.abs(S21 / S12) if np.abs(S12) > 0 else np.inf

def gufm(S):
    S11 = S[0,0]
    S12 = S[0,1]
    S21 = S[1,0]
    S22 = S[1,1]
    U = (S12 * S21 * np.conj(S11 * S22)) * gui(S) * guo(S)
    return 1 / np.abs(1 - U)**2

def s2abcd(S, Z0=50):
    S11 = S[0,0]
    S21 = S[1,0]
    S12 = S[0,1]
    S22 = S[1,1]
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
    S11 = S[0,0]
    S21 = S[1,0]
    S12 = S[0,1]
    S22 = S[1,1]
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
    S11 = S[0,0]
    S21 = S[1,0]
    S31 = S[2,0]
    S12 = S[0,1]
    S22 = S[1,1]
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
    S11 = S[0,0]
    S21 = S[1,0]
    S31 = S[2,0]
    S12 = S[0,1]
    S22 = S[1,1]
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
    S11 = S[0,0]
    S21 = S[1,0]
    S31 = S[2,0]
    S12 = S[0,1]
    S22 = S[1,1]
    S32 = S[2,1]
    S13 = S[0,2]
    S23 = S[1,2]
    S33 = S[2,2]
    return np.array([
       [ S11-S12*S21/(1+S22), S13-S12*S23/(1+S22) ],
       [ S31-S32*S21/(1+S22), S33-S32*S23/(1+S22) ]
    ])


##########################

def read_network(path=None):
    if path is None:
        buf = sys.stdin.read()
        path = tempfile.mktemp() + ".s2p"
        with open(path, "w") as f:
            f.write(buf)
        nw = rf.Network(path)
        os.unlink(path)
    else:
        nw = rf.Network(path)
    if np.any(nw.z0 != 50):
        print('Only networks referenced to 50 ohms supported', file=sys.stderr)
        sys.exit(1)
    return nw

def write_network(nw, mode):
    polar = lambda x: "{:9.4g} {:7.2f}".format(np.abs(x), np.angle(x) * 180 / np.pi)
    if mode == 'a':
        print('MHZ            A                 B                 C                 D')
        for i in range(len(nw)):
            f = nw.f[i] / 1e6
            S = nw.s[i]
            data = ' '.join([ polar(x) for x in s2abcd(S).flatten() ])
            print('{:<5g}'.format(f), data)
    elif mode == 'g':
        print('MHZ      GUM    GUI    GUO     gu   GMSG')
        for i in range(len(nw)):
            f = nw.f[i] / 1e6
            S = nw.s[i]
            K = nw.stability[i]
            # Gmsg = '{:6.2f}'.format(db(gmsg(S))) if K < 1 else '   -'
            print('{:<5g} {:6.2f} {:6.2f} {:6.2f} {:6.2f} {:6.2f}'.format(
                  f, db(gu(S)), db(gui(S)), db(guo(S)), db(gufm(S)), db(gmsg(S))))
    elif mode == 'z':
        print('MHZ           ZIN             ZOUT')
        for i in range(len(nw)):
            f = nw.f[i] / 1e6
            S = nw.s[i]
            print('{:<5g} {:16.4g} {:16.4g}'.format(f, g2z(S[0,0]), g2z(S[1,1])))
    elif mode == 'n':
        print(nw.write_touchstone(form='ma', return_string=True))
    else:
        print('# MHZ S MA R 50')
        print('! MHZ         S11               S21               S12               S22       '
              '!   GUM         K         D')
        for i in range(len(nw)):
            f = nw.f[i] / 1e6
            S = nw.s[i]
            K = nw.stability[i]
            D = np.abs(S[0,0] * S[1,1] - S[0,1] * S[1,0])
            flag = '' if K > 1 and D < 1 else 'pu'
            data = ' '.join([ polar(x) for x in S.T.flatten() ])
            print("{:<5g} {:s} ! {:5.1f} {:9.4g} {:9.4g}".format(f, data, db(gu(S)), K, D), flag)


def main(*args):
    args = list(args)
    mode = None
    stack = []
    stack.append(read_network())

    while args:
        opt = args.pop(0)
        top = stack[-1] if stack else None

        if opt == '-n':
            mode = "n"
        elif opt == '-a':
            mode = "a"
        elif opt == '-z':
            mode = "z"
        elif opt == '-g':
            mode = "g"
        
        elif opt == "-swap":
            b = stack.pop()
            a = stack.pop()
            stack.append(b)
            stack.append(a)
        elif opt == "-cascade":
            b = stack.pop()
            a = stack.pop()
            stack.append(a ** b)
        elif opt == "-deembed":
            b = stack.pop()
            a = stack.pop()
            stack.append(a ** b.inv)
        elif opt == "-ideembed":
            b = stack.pop()
            a = stack.pop()
            stack.append(a.inv ** b)

        elif opt == "-p":
            write_output(top, mode=mode)
        elif opt == "-f":
            stack.append(read_network(args.pop(0)))
        elif top and opt == "-series":
            S = abcd2s([[1, float(args.pop(0))], [0, 1]])
            stack.append(rf.Network(frequency=top.frequency, s=[S] * len(top)))
        elif top and opt == "-shunt":
            S = abcd2s([[1, 0], [1/float(args.pop(0)), 1]])
            stack.append(rf.Network(frequency=top.frequency, s=[S] * len(top)))

        elif top and opt == "-cbg":
            top.s = np.array([ cbg_transform(S) for S in top.s ])
        elif top and opt == "-ccd":
            top.s = np.array([ ccd_transform(S) for S in top.s ])
        elif top and opt == "-lift":
            x = args.pop(0)
            top.s = np.array([ lift_ground(
                top.s[i], 
                complex(x) if 'j' in x else 2j * np.pi * top.f[i] * float(x)
            ) for i in range(len(top)) ])
        else:
            print('Unrecognized command line option. Exiting.', file=sys.stderr)
            sys.exit(1)

    if stack: 
        write_network(stack[-1], mode=mode)


if __name__ == "__main__":
    np.seterr(divide='ignore')
    main(*sys.argv[1:])


