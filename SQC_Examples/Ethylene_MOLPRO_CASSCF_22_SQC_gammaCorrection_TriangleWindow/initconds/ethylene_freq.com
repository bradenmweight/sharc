***,ethylene freqency calculation
memory,8,m

file,2,ethylene.wf,new;

basis=6-31++G**
gthresh,energy=1.d-8

geomtyp=xyz
symmetry,nosym
geometry={
   6
       ethylene CARTESIAN COORDINATES
 C    0.0000000031    0.0000000000    0.6768496105
 C   -0.0000000013    0.0000000000   -0.6768496105
 H   -0.0000000017    0.9131426983   -1.2402810929
 H   -0.0000000095    0.9131426983    1.2402810929
 H   -0.0000000159   -0.9131426983    1.2402810929
 H    0.0000000047   -0.9131426983   -1.2402810929
}

hf;accu,14
optg;coord,3n;

{frequencies,analytic
thermo,sym=auto
print,thermo}

mp2
optg;coord,3n
{frequencies
thermo,sym=c2v
print,thermo}
put,molden,ethylene.molden;
