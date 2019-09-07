
##########################################

cat > /tmp/cascade_input <<EOF
# mhz s ma r 50
1 .928 -64 10.84 150 .023 70 .529 -27
EOF

cat > /tmp/cascade_output <<EOF
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       !   GUM         K        MU
1         0.928  -64.00     10.84  150.00     0.023   70.00     0.529  -27.00 !  30.7  0.003077    0.2598
EOF
< /tmp/cascade_input python cascade.py | diff - /tmp/cascade_output

cat > /tmp/cascade_output <<EOF
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       !   GUM         K        MU
1        0.9962   -3.97    0.1462   86.68   0.07046   90.57         1   -3.88 !   nan    0.9975    0.5713
EOF
< /tmp/cascade_input python cascade.py -lift 1250j | diff - /tmp/cascade_output

cat > /tmp/cascade_output <<EOF
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       !   GUM         K        MU
1          0.75 -178.23     1.727   -4.26   0.06137    2.79    0.9448   -4.14 !  18.0     0.983    0.9946
EOF
< /tmp/cascade_input python cascade.py -cbg | diff - /tmp/cascade_output

cat > /tmp/cascade_output <<EOF
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       !   GUM         K        MU
1         3.966   18.75     2.995 -163.94     2.043   28.93     1.203 -137.25 !   1.3   -0.9763   -0.9837
EOF
< /tmp/cascade_input python cascade.py -cbg -lift 600j | diff - /tmp/cascade_output

#######################################

cat > /tmp/cascade_input <<EOF
# mhz s ma r 50
1 .495 -158 2.55 75 .132 45 .415 -52
EOF

cat > /tmp/cascade_output <<EOF
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       !   GUM         K        MU
1         0.495 -158.00      2.55   75.00     0.132   45.00     0.415  -52.00 !  10.2    0.9187    0.9417
EOF
< /tmp/cascade_input python cascade.py | diff - /tmp/cascade_output

cat > /tmp/cascade_output <<EOF
MHZ            A                 B                 C                 D
1       0.08687  -26.78     10.33 -116.70  0.004851  -59.83    0.3382  -73.43
EOF
< /tmp/cascade_input python cascade.py -a | diff - /tmp/cascade_output

cat > /tmp/cascade_output <<EOF
MHZ            A                 B                 C                 D
1       0.09409  150.79     11.19   60.87  0.005254  117.74     1.012   18.75
EOF
< /tmp/cascade_input python cascade.py -cbg -a | diff - /tmp/cascade_output

cat > /tmp/cascade_output <<EOF
# MHZ S MA R 50
! MHZ         S11               S21               S12               S22       !   GUM         K        MU
1        0.5361 -133.70     3.851   47.30  0.007624  172.83     1.119  -30.14 !   nan    -3.588    0.8473
EOF
< /tmp/cascade_input python cascade.py -cbg > /tmp/cascade_cbg.s2p
< /tmp/cascade_input python cascade.py -f /tmp/cascade_cbg.s2p -cascade | diff - /tmp/cascade_output
