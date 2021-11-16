
import neuron
from neuron import h, rxd
from neuron.units import nM, uM
import numpy
from matplotlib import pyplot
from neuron import crxd as rxd
from neuron.crxd import v
from neuron.crxd.rxdmath import vtrap, exp, log
from math import pi
h.load_file('stdrun.hoc')

#Parameters
h.celsius = 37.0
k = 1.38064852e-23 #Boltzmann constant
T = h.celsius + 273.15 #Temperature in Kelvin
gnabar = 3.39291720 #uS
gkbar = 1.01787516 #uS
gcabar_L = 25e-6 #uS
gcabar_N = 14e-6 #uS
gcabar_T = 8e-6 #uS
gl = 0.001413716 #uS, assuming 30 um soma diameter and transmembrane resistance of 20000 Ohm*cm^2
el = -70
q10 = 3.0**((h.celsius - 37.0)/10.0)
"""
temp1 = 0.0853*T/2
temp2 = v/temp1
temp3 = 0
if (abs(temp2) < 1.10e-4):
    temp3 = 1-temp2/2
else:
    temp3 = temp2/(exp(temp2)-1)
"""
#Na channel
#m-gate
alpha = -0.6*(v + 30)/(exp(-0.1*(v + 30)) - 1)
beta = 20 * exp(-(v+55)/18.0)
mtau = 1/(q10*(alpha + beta))
minf = alpha/(alpha + beta)

"""
Note: "vtrap(x,y) is 1/(exp(x/y)-1) if |x/y|>=1e-6 or y*(1.0 - x/y/2.0) otherwise" 
https://neuron.yale.edu/neuron/docs/hodgkin-huxley-using-rxd
"""

#h-gate
alpha = 0.4*exp(-(v + 50.0)/20.0)
beta = 6.0/(exp(-0.1*(v + 20.0)) + 1.0)
htau = 1/(q10 * (alpha + beta))
hinf = alpha/(alpha + beta)

#K channel
#n-gate
alpha = -0.02 * (v + 40)/(exp(-0.1 * (v + 40)) - 1)
beta = 0.4 * exp(-(v + 50.0)/80.0)
ntau = 1.0/(q10 * (alpha + beta))
ninf = alpha/(alpha + beta)


#Ca channel

#T-channel
#m-gate
alpha = 0.2*vtrap((-v + 19.26), 10)
beta = 9.0e-3*exp(-v/22.03)
minf_T = alpha/(alpha+beta)
tau_Tm = 1/(q10*(alpha + beta))
#h-gate
alpha = 1e-6*exp(-v/16.26)
beta = 1/(exp((-v + 29.79)/10) + 1)
hinf_T = alpha/(alpha + beta)
tau_Th = 1/(q10*(alpha + beta))

#N-channel
#mgate
alpha = 0.19*vtrap((-v + 19.88), 10)
beta = 4.6e-2*exp(-v/20.73)
minf_N = alpha/(alpha + beta)
tau_Nm = 1/(q10*(alpha + beta))
#h-gate
alpha = 1.6e-4*exp(-v/48.4)
beta = 1/(exp((-v + 39)/10) + 1)
hinf_N = alpha/(alpha + beta)
tau_Nh = 1/(q10*(alpha + beta))

#L-channel
alpha = 15.69*vtrap((-v + 81.5), 10)
beta = 0.29*exp(-v/10.86)
minf_L = alpha/(alpha + beta)
tau_Lm = 1/(q10*(alpha + beta))


soma1 = h.Section('somaA')
soma1.pt3dclear()
soma1.pt3dadd(-90, 0, 0, 30)
soma1.pt3dadd(-60, 0, 0, 30)
soma1.nseg = 11

#Where?
cyt = rxd.Region(h.allsec(), name = 'cyt', nrn_region = 'i')
mem = rxd.Region(h.allsec(), name = 'cell_mem', geometry = rxd.membrane())
ecs = rxd.Extracellular(-100, -100, -100, 100, 100, 100, dx = 33)


#Who?
#Intracellular Na, K, Ca
na = rxd.Species([cyt, mem], name = 'na', d = 1, charge = 1, initial = 10)
k = rxd.Species([cyt, mem], name = 'k', d = 1, charge = 1, initial = 54.4)
ca = rxd.Species([cyt, mem], name = 'ca', d = 1, charge = 2, initial = 0.00005)

#Extracellular Na, K, Ca
kecs = rxd.Parameter(ecs, name = 'k', charge = 1, value = 2.5)
naecs = rxd.Parameter(ecs, name = 'na', charge = 1, value = 140.0)
caecs = rxd.Parameter(ecs, name = 'ca', charge = 2, value = 2.0)

#Leak
x = rxd.Species([cyt, mem, ecs], name = 'x', charge = 1)

ki = k[cyt]
ko = kecs[ecs]
nai = na[cyt]
nao = naecs[ecs]
xi = x[cyt]
xo = x[ecs]

cai = ca[cyt]
cao = caecs[ecs]

#gating states -> open probabilities
ngate = rxd.State([cyt, mem], name = 'ngate', initial = 0.2445865)
mgate = rxd.State([cyt, mem], name = 'mgate', initial = 0.0289055)
hgate = rxd.State([cyt, mem], name = 'hgate', initial = 0.7540796)

mgate_T = rxd.State([cyt, mem], name = 'mgate_T', initial = 0.010871215)
hgate_T = rxd.State([cyt, mem], name = 'hgate_T', initial = 0.615047335)
mgate_N = rxd.State([cyt, mem], name = 'mgate_N', initial = 0.001581552)
hgate_N = rxd.State([cyt, mem], name = 'hgate_N', initial = 0.973556944)
mgate_L = rxd.State([cyt, mem], name = 'mgate_L', initial = 3.42574e-6)


#Who?
#Intracellular Na, K, Ca
na = rxd.Species([cyt, mem], name = 'na', d = 1, charge = 1, initial = 10)
k = rxd.Species([cyt, mem], name = 'k', d = 1, charge = 1, initial = 54.4)
ca = rxd.Species([cyt, mem], name = 'ca', d = 1, charge = 2, initial = 0.00005)

#Extracellular Na, K, Ca
kecs = rxd.Parameter(ecs, name = 'k', charge = 1, value = 2.5)
naecs = rxd.Parameter(ecs, name = 'na', charge = 1, value = 140.0)
caecs = rxd.Parameter(ecs, name = 'ca', charge = 2, value = 2.0)

#Leak
x = rxd.Species([cyt, mem, ecs], name = 'x', charge = 1)

ki = k[cyt]
ko = kecs[ecs]
nai = na[cyt]
nao = naecs[ecs]
xi = x[cyt]
xo = x[ecs]

cai = ca[cyt]
cao = caecs[ecs]

#gating states -> open probabilities
ngate = rxd.State([cyt, mem], name = 'ngate', initial = 0.2445865)
mgate = rxd.State([cyt, mem], name = 'mgate', initial = 0.0289055)
hgate = rxd.State([cyt, mem], name = 'hgate', initial = 0.7540796)

mgate_T = rxd.State([cyt, mem], name = 'mgate_T', initial = 0.010871215)
hgate_T = rxd.State([cyt, mem], name = 'hgate_T', initial = 0.615047335)
mgate_N = rxd.State([cyt, mem], name = 'mgate_N', initial = 0.001581552)
hgate_N = rxd.State([cyt, mem], name = 'hgate_N', initial = 0.973556944)
mgate_L = rxd.State([cyt, mem], name = 'mgate_L', initial = 3.42574e-6)


stim = h.IClamp(soma1(0.5))
stim.delay = 5
stim.amp = 1
stim.dur = 50

tvec = h.Vector().record(h._ref_t)

vvec = h.Vector().record(soma1(0.5)._ref_v)

mvec = h.Vector().record(mgate[cyt].nodes(soma1(0.5))._ref_value)
nvec = h.Vector().record(ngate[cyt].nodes(soma1(0.5))._ref_value)
hvec = h.Vector().record(hgate[cyt].nodes(soma1(0.5))._ref_value)
kvec = h.Vector().record(soma1(0.5)._ref_ik)
navec = h.Vector().record(soma1(0.5)._ref_ina)

mT_vec = h.Vector().record(mgate_T[cyt].nodes(soma1(0.5))._ref_value)
hT_vec = h.Vector().record(hgate_T[cyt].nodes(soma1(0.5))._ref_value)

mN_vec = h.Vector().record(mgate_N[cyt].nodes(soma1(0.5))._ref_value)
hN_vec = h.Vector().record(hgate_N[cyt].nodes(soma1(0.5))._ref_value)

mL_vec = h.Vector().record(mgate_L[cyt].nodes(soma1(0.5))._ref_value)


h.finitialize(-70)
h.continuerun(100)
