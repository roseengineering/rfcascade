
import os, subprocess 

def run(command):
    proc = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    buf = proc.stdout.read().decode()
    proc.wait()
    return f"""
```
$ {command}
{buf}\
```
"""

print(f"""
Cascade
-------
Software utility to manipulate touchstone s-parameter files.

Features
--------
The utility takes a touchstone file as standard input
and outputs a new touchstone file (or other formats)
to standard output.  Cascade provides the following
input transformations as command line options:

```
-cascade             : cascade together the two networks on top of stack
-deembed             : de-embed the output of the top two networks on the stack
-ideembed            : de-embed the input of the top two networks on the stack
-swap                : swap the top two networks on the stack
-cbg                 : transform network on top of stack into a common-base arrangement
-ccd                 : transform network on top of stack into a common-collector arrangement
-lift <complex>      : lift network on top of stack from ground and insert an impedance, j required
-lift <henries>      : lift network on top of stack from ground and insert an inductor
-p                   : print network on top of stack
-series <ohms>       : push a series resistor onto the stack
-shunt <ohms>        : push a shunt resistor onto the stack
-f <filename>        : push a touchstone network onto the stack
-gs <complex>        : set the source gamma for matching
-zs <complex>        : set the source impedance for matching
-gl <complex>        : set the load gamma for matching
-zl <complex>        : set the load impedance for matching
```

By default the utility writes out the network on the top of
the stack in touchstone format with GUM and Rollet stability information 
as comments.  It can also output the network in following alternative formats:

```
-a      : display the network as ABCD matrices
-s      : summarize the network in terms of impedances, stability and gain (dB) values
-lmatch : match using l-section networks
-stub   : match using single-stub networks with open shunt stubs
-qwt2   : match using quarter wavelength networks with open shunt stubs
```

Only 50 ohm networks are supported.

Installation
------------
```
$ sudo pip install -r requirements.txt
$ sudo cp cascade.py /usr/local/bin/cascade
$ sudo chmod 755 /usr/local/bin/cascade
$ sh test.sh       # to run unit tests
```

Examples
--------

Add GUM (dB), K and mu as comments to touchstone file.

{ run("< 2n5179_5ma.s2p cascade") }

Display the result as ABCD matrices.

{ run("< 2n5179_5ma.s2p cascade -a") }

Summarize the network.  GU is not in dB.

{ run("< 2n5179_5ma.s2p cascade -s") }

Shunt a 330 ohm resistor across the output of the two-port network.  This makes the transistor unconditionally stable.

{ run("< 2n5179_5ma.s2p cascade -shunt 330 -cascade") }

Cascade a series 20 ohm resistor with the output of the two-port network.

{ run("< 2n5179_5ma.s2p cascade -series 20 -cascade") }

Lift terminal 3, the port connected to ground, and add a 10nH inductor.  Then cascade the result
with a shunt 100 ohm resistor to stabilize the result.

{ run("< 2n5179_5ma.s2p cascade -lift 10e-9 -shunt 100 -cascade") }

Insert a one ohm resistor at the emitter to provide shunt feedback.

{ run("< 2n5179_5ma.s2p cascade -lift 1+0j") }

Transform the common emitter s-parameter input into a common-base.

{ run("< 2n5179_5ma.s2p cascade -cbg | tee cb.s2p") }

Create a cascode amplifier.  S12 is significantly reduced compared to a CE amp.

{ run("< 2n5179_5ma.s2p cascade -f cb.s2p -cascade") }

Stabilize the cascode amp with a 100 ohm resistor across the output.

{ run("< 2n5179_5ma.s2p cascade -f cb.s2p -cascade -shunt 100 -cascade") }

Summarize the stabilized cascode amp.

{ run("< 2n5179_5ma.s2p cascade -f cb.s2p -cascade -shunt 100 -cascade -s") }
""")

os.unlink("cb.s2p")

