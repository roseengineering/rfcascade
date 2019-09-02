
Cascade
-------
Software utility to manipulate touchstone s-parameter files.

Features
--------
The utility takes a touchstone file as standard input
and outputs a new touchstone file (or other formats)
as standard output.  Cascade provides the following
input transformations provided as command line options:

```
-cascade <filename>  : cascades the input with the given touchstone file.
-cbg                 : transforms the input into a common-base arrangement.
-ccd                 : transforms the input into a common-collector arrangement.
-iseries <ohms>      : cascades a series resistor with the input.
-oseries <ohms>      : cascades the input with a series resistor.
-oshunt <ohms>       : cascades the input with a shunt resistor.
-lift <complex>      : lifts ground from input and inserts a complex impedance
-lift <inductance>   : lifts ground from input and inserts an inductor
```

By default the utlitilty outputs a touchstone file with 
GUM and Rollets stability information as comments.  It can
also output the following alternative formats:

```
-n  : outputs the results of write_touchstone() from scikit-rf 
-a  : outputs the results as a ABCD matrix
-z  : outputs the impedance of S11 and S22 normalized to 50 ohms
```

Installation
------------
$ sudo pip install -r requirements.txt
$ sudo cp cascade.py /usr/local/bin/cascade
$ sudo chmod 755 /usr/local/bin/cascade
$ sh test.sh       # to run unit tests

Examples
--------

The pu flag means potentially unstable.

```
$ cascade < 2n5179_5ma.s2p 
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       ! GUM[dB]       K         D
0         0.471    0.00      6.78  180.00         0    0.00     0.844    0.00 !  23.1       inf    0.3975 
100       0.471  -90.00      6.78  122.00     0.023   64.00     0.844  -51.00 !  23.1    0.4623    0.2799 pu
200       0.314 -145.00       4.2  100.00     0.034   58.00      0.78  -93.00 !  17.0     1.109    0.1542 
300        0.23  156.00      2.76   91.00     0.043   65.00     0.768 -134.00 !  12.9     1.819    0.2728 
400       0.171  108.00      2.17   86.00     0.056   63.00     0.756 -177.00 !  10.5     1.874    0.2371 
500       0.168   54.00      1.86   79.00     0.062   62.00     0.741  140.00 !   9.0     1.883    0.1073 
600       0.149   -9.00      1.53   71.00     0.069   66.00      0.74   98.00 !   7.2     2.074   0.08789 
700       0.137  -72.00      1.31   67.00     0.073   71.00     0.739   54.00 !   5.9     2.469    0.1926 
800       0.119 -129.00      1.18   64.00     0.092   74.00     0.744    8.00 !   5.0     2.098    0.1526 
900       0.153 -174.00      1.13   58.00     0.101   68.00     0.742  -38.00 !   4.6     1.875   0.04344 
1000      0.171  122.00     0.979   49.00     0.106   71.00     0.749  -82.00 !   3.5     2.083    0.1502 
```

Display the result as ABCD matrices.

```
$ cascade -a < 2n5179_5ma.s2p 
! MHZ          A                 B                 C                 D
0       0.01692 -180.00        10 -180.00 0.0001217 -180.00   0.07194 -180.00
100     0.05534  -88.08     7.139 -166.65  0.001397  -51.43    0.1243 -120.34
200       0.107  -70.35     6.367 -148.80  0.004076  -59.44      0.17 -131.82
300      0.2218  -60.14     6.235 -130.84  0.007505  -77.69    0.1408 -149.65
400      0.3684  -72.05     4.008  -95.97  0.009134  -94.89   0.03424  -93.31
500      0.4608  -86.97     9.595  -33.46  0.008633 -105.44    0.1543  -28.46
600      0.4624 -105.15     22.05  -37.35  0.008063 -103.97    0.3211  -24.06
700      0.2943 -122.15     32.17  -54.24  0.006801 -106.18    0.5633  -33.08
800     0.06799  -81.22     36.07  -68.63   0.00344  -87.15    0.7668  -53.02
900       0.246    0.28        33  -77.79  0.006156  -18.70    0.8021  -71.03
1000     0.5662    4.77     33.48  -76.44   0.01307  -21.60    0.6974  -89.08
```

Display the result as the impedances of S11 and S22 normalized to 50 ohms.

```
$ cascade -z < 2n5179_5ma.s2p 
! MHZ         Z11               Z22
0           139    0.00       591    0.00
100          50  -50.44     103.3  -77.63
200       30.09  -21.78     47.52  -75.89
300       32.77   11.17     22.18  -69.63
400        45.1   18.52      7.07  -10.46
500       60.73   15.63     19.63   64.67
600       67.24   -2.73     43.73   72.85
700       54.34  -14.87     94.41   69.21
800       43.09  -10.63     307.5   24.88
900        36.8   -1.88     133.6  -63.81
1000      41.85   16.63     57.19  -73.51
```

Shunt a 330 ohm resistor across the two-port output.  This makes the transistor unconditionally stable.

```
$ cascade -oshunt 330 < 2n5179_5ma.s2p 
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       ! GUM[dB]       K         D
0         0.471   -0.00     5.949  180.00 4.877e-18 -180.00     0.618   -0.00 !  18.7 8.289e+15    0.2911 
100      0.4695  -88.72     6.069  124.55   0.02059   66.55    0.6577  -53.05 !  19.2     1.557    0.2057 
200      0.3082 -143.49      3.91  103.15   0.03165   61.15    0.6784  -95.81 !  15.0     1.884    0.1468 
300      0.2213  155.91     2.664   93.31    0.0415   67.31    0.7377 -135.77 !  12.1     2.142    0.2583 
400      0.1643  105.92      2.13   86.17   0.05498   63.17    0.7603 -177.12 !  10.4     1.906    0.2271 
500      0.1675   51.12       1.8   77.00      0.06   60.00    0.7204  141.75 !   8.4     2.151    0.1077 
600       0.155  -10.66     1.431   68.02   0.06452   63.02    0.6532  101.19 !   5.6         3   0.06756 
700      0.1428  -70.78     1.181   64.66   0.06579   68.66    0.5781   57.14 !   3.3     4.307    0.1536 
800      0.1197 -125.52     1.043   63.60    0.0813   73.60    0.5415    8.59 !   1.9     4.169    0.1199 
900      0.1491 -171.48     1.008   59.77   0.09013   69.77    0.5603  -40.49 !   1.8     3.657   0.02895 
1000     0.1638  121.96    0.9022   51.97   0.09769   73.97     0.632  -85.31 !   1.4     3.358    0.1351 
```

Cascade a series a 20 ohm resistor on the two-port output.

```
$ cascade -oseries 20 < 2n5179_5ma.s2p 
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       ! GUM[dB]       K         D
0         0.471   -0.00     6.575  180.00 2.695e-18 -180.00    0.8487   -0.00 !  23.0 6.141e+15    0.3997 
100      0.4714  -93.44     6.155  115.16   0.02088   57.16    0.7407  -46.86 !  20.3     1.194    0.2787 
200      0.3248 -148.73     3.448   92.65   0.02791   50.65    0.5297  -82.32 !  12.7     3.225   0.08236 
300       0.248  155.65     2.105   86.17   0.03279   60.17    0.3788 -122.00 !   7.4     5.894    0.1364 
400       0.185  111.63     1.606   85.66   0.04145   62.66    0.2999 -175.86 !   4.7     6.681    0.1172 
500      0.1686   59.96     1.412   83.15   0.04708   66.15    0.3478  127.85 !   3.7      6.41   0.04192 
600      0.1339   -5.65     1.245   77.85   0.05613   72.85    0.4864   85.50 !   3.2      5.38   0.07833 
700      0.1236  -76.67      1.17   73.13   0.06521   77.13    0.6494   47.27 !   3.8      3.85    0.1566 
800      0.1201 -138.88     1.121   65.13   0.08738   75.13    0.7539    7.12 !   4.7     2.218     0.131 
900       0.163  179.31      1.04   53.18   0.09292   63.18    0.7004  -33.51 !   3.4     2.516   0.05619 
1000     0.1883  121.15    0.8238   41.83   0.08919   63.83    0.5533  -71.64 !   0.1     4.533   0.08788 
```

Lift terminal 3, the port connect to ground, and add an 10nH inductor.  Then cascade the result
with a shunt 100 ohm resistor to stabilize the result.

```
$ cascade -lift 10e-9 -oshunt 100
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       ! GUM[dB]       K         D
0         0.471    0.00     4.641  180.00  1.49e-33   90.00    0.2621   -0.00 !  14.7 5.241e+31    0.1235 
100      0.2894  -31.98      3.67  112.81   0.03295  129.17    0.2602  -55.49 !  12.0     3.527    0.0679 
200      0.2813    0.15     2.299   91.42    0.1193  112.36    0.2716  -97.16 !   7.9     1.653    0.2441 
300      0.3553    4.07     1.669   75.50    0.2257   82.69    0.3018 -146.61 !   5.4     1.188    0.3348 
400      0.3673   -2.21     1.266   61.02    0.2619   55.24    0.2773  162.64 !   3.0     1.297    0.2681 
500      0.3933   -3.40    0.9646   51.78    0.2325   39.45    0.2338  118.67 !   0.7      1.81    0.1452 
600      0.3965  -14.35    0.7704   47.45     0.207   37.88    0.2079   80.31 !  -1.3      2.53   0.08616 
700       0.319  -21.21    0.6623   46.78    0.1807   37.44    0.2023   41.51 !  -2.9     3.631    0.1082 
800      0.2797   -3.83    0.5861   50.23    0.1252   54.88    0.2335    0.19 !  -4.0     5.997    0.1128 
900       0.375    0.32    0.6423   57.44    0.2209  100.88    0.3044  -46.80 !  -2.8     2.922    0.2501 
1000     0.4493  -17.00    0.8127   50.80    0.4842   83.82    0.3872 -103.62 !  -0.1     1.103    0.4691 
```


