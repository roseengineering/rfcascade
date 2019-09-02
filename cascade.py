#!/usr/bin/python3

import numpy as np
import sys, tempfile, os

def tothreeport(S):
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
    G = (Z - Z0) / (Z + Z0)
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


def read_input(buf=None):
    buf = buf or sys.stdin.read()
    path = tempfile.mktemp() + ".s2p"
    with open(path, "w") as f:
        f.write(buf)
    nw = rf.Network(path)
    os.unlink(path)
    return nw


def write_output(nw, mode):
    polar = lambda x: "{:9.4g} {:7.2f}".format(np.abs(x), np.angle(x) * 180 / np.pi)
    power = lambda x: np.abs(x)**2
    imped = lambda x: 50 * (1 + x) / (1 - x)
    if mode == 'a':
        print('! MHZ          A                 B                 C                 D')
        for i in range(len(nw)):
            print('{:<5g}'.format(nw.f[i] / 1e6), ' '.join([ polar(x) for x in nw.a[i].flatten() ]))
    elif mode == 'z':
        print('! MHZ         Z11               Z22')
        for i in range(len(nw)):
            S = nw.s[i]
            print('{:<5g}'.format(nw.f[i] / 1e6), polar(imped(S[0,0])), polar(imped(S[1,1]))) 
    elif mode == 'n':
        print(nw.write_touchstone(form='ma', return_string=True))
    else:
        print('# MHZ S MA R 50')
        print('! MHZ         S11               S21               S12               S22       '
              '! GUM[dB]       K         D')
        for i in range(len(nw)):
            S = nw.s[i]
            P11, P22, P21 = power(S[0,0]), power(S[1,1]), power(S[1,0])
            GUM = 10 * np.log10(P21 / (1 - P11) / (1 - P22)) if P11 < 1 and P22 < 1 else np.inf
            K = nw.stability[i]
            D = np.abs(S[0,0] * S[1,1] - S[0,1] * S[1,0])
            flag = '' if K > 1 and D < 1 else 'pu'
            data = ' '.join([ polar(x) for x in S.T.flatten() ])
            print("{:<5g} {:s} ! {:5.1f} {:9.4g} {:9.4g}".format(nw.f[i] / 1e6, data, GUM, K, D), flag)


def main(*args):
    args = list(args)
    mode = None
    nw = read_input()
    while args:
        opt = args.pop(0)
        if opt == '-n':
            mode = "n"
        elif opt == '-a':
            mode = "a"
        elif opt == '-z':
            mode = "z"
        elif opt == "-cascade":
            nw **= rf.Network(args.pop(0))
        elif opt == "-cbg":
            nw.s = np.array([ cbg_transform(S) for S in nw.s ])
        elif opt == "-ccd":
            nw.s = np.array([ ccd_transform(S) for S in nw.s ])
        elif opt == "-iseries":
            M = np.matrix([[1, float(args.pop(0))], [0, 1]]) 
            nw.a = np.array([ M * np.matrix(nw.a[i]) for i in range(len(nw)) ])
        elif opt == "-oseries":
            M = np.matrix([[1, float(args.pop(0))], [0, 1]]) 
            nw.a = np.array([ np.matrix(nw.a[i]) * M for i in range(len(nw)) ])
        elif opt == "-oshunt":
            M = np.matrix([[1, 0], [1/float(args.pop(0)), 1]])
            nw.a = np.array([ np.matrix(nw.a[i]) * M for i in range(len(nw)) ])
        elif opt == "-lift":
            z = args.pop(0)
            nw.s = np.array([ lift_ground(nw.s[i], 
                complex(z) if 'j' in z else 2 * np.pi * nw.f[i] * float(z) * 1j) 
                for i in range(len(nw)) ])
        else:
            print('unrecognized command line options, exiting', file=sys.stderr)
            return 1
    write_output(nw, mode=mode)


import skrf as rf 
if __name__ == "__main__":
    np.seterr(divide='ignore')
    sys.exit(main(*sys.argv[1:]))


