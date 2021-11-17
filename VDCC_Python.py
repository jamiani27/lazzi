# -*- coding: utf-8 -*-
"""
Created on Tue Nov 16 08:40:30 2021
- Based on Sangwon's file - with Arthur's scale suggestions
- Fixed the V in denominator that halted the simulation when V=0mV (from line 195)

Still remains to be done:
- Low Na/Ca current
- Verify validity of the calculation for the calcium current (from line 195)
- Use Voltage Clamp and plot Ca current with 3 different channels
- Modularize the code. A good example can be found here:
  https://neuron.yale.edu/neuron/static/docs/neuronpython/ballandstick3.html

@author:Sangwon, Arthur, jmbouteiller
"""
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
e = 1.60217662e-19
scale = 1e-14/e #convert S to molecules/um2
gnabar = 0.04*scale # molecules/um2 ms mV (was 3.39291720 #uS)
gkbar = 0.012*scale # molecules/um2 ms mV (was 1.01787516 #uS)
gcabar_L = 25e-12*scale ## molecules/um2 ms mV (was 25e-6 #uS)
gcabar_N = 14e-12*scale # molecules/um2 ms mV (was 14e-6 #uS)
gcabar_T = 8e-12*scale # molecules/um2 ms mV (was 8e-6 #uS)
gl = 0.0003*scale # (was 0.001413716 #uS) molecules/um2 ms mV, assuming 30 um soma diameter and transmembrane resistance of 20000 Ohm*cm^2
el = -70
q10 = 3.0**((h.celsius - 37.0)/10.0)


#Na channel
#m-gate
alpha = 0.1 * vtrap(-(v + 40.0), 10)
beta = 4.0 * exp(-(v + 65)/18.0)
mtau = 1.0/(q10 * (alpha + beta))
minf = alpha/(alpha + beta)

"""
Note: "vtrap(x,y) is 1/(exp(x/y)-1) if |x/y|>=1e-6 or y*(1.0 - x/y/2.0) otherwise"
https://neuron.yale.edu/neuron/docs/hodgkin-huxley-using-rxd
"""

#h-gate
alpha = 0.07 * exp(-(v + 65.0)/20.0)
beta = 1.0/(exp(-(v + 35.0)/10.0) + 1.0)
htau = 1.0/(q10 * (alpha + beta))
hinf = alpha/(alpha + beta)

#K channel
#n-gate
alpha = 0.01 * vtrap(-(v + 55.0), 10.0)
beta = 0.125 * exp(-(v + 65.0)/80.0)
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
ca = rxd.Species([cyt, mem], name = 'ca', d = 1, charge = 2, initial = 0.0001)

#Extracellular Na, K, Ca
kecs = rxd.Parameter(ecs, name = 'k', charge = 1, value = 2.5)
naecs = rxd.Parameter(ecs, name = 'na', charge = 1, value = 140.0)
caecs = rxd.Parameter(ecs, name = 'ca', charge = 2, value = 1.0)

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

#What?
m_gate = rxd.Rate(mgate, (minf-mgate)/mtau)
h_gate = rxd.Rate(hgate, (hinf-hgate)/htau)
n_gate = rxd.Rate(ngate, (ninf-ngate)/ntau)

m_gate_T = rxd.Rate(mgate_T, (minf_T-mgate_T)/tau_Tm)
h_gate_T = rxd.Rate(hgate_T, (hinf_T-hgate_T)/tau_Th)
m_gate_N = rxd.Rate(mgate_N, (minf_N-mgate_N)/tau_Nm)
h_gate_N = rxd.Rate(hgate_N, (hinf_N-hgate_N)/tau_Nh)
m_gate_L = rxd.Rate(mgate_L, (minf_L-mgate_L)/tau_Lm)


"""
  NOTE: 'mgate' and 'm_gate' are different. 'mgate' is an object of class rxd.State, representing
  open probability value. We are telling NEURON 'who' we are interested in, which in this case are the
  various types of gates.
  'm_gate'is an rxd.Rate object which we are telling NEURON to track over time according
  to the equation (minf-mgate)/mtau. Using the rxd.Rate class allows the rxd.State objects
  to change over time.
"""

#Nernst Potentials
Ena = 1e3 * h.R * (h.celsius + 273.15) * log(nao/nai)/h.FARADAY
Ek = 1e3 * h.R * (h.celsius + 273.15) * log(ko/ki)/h.FARADAY
Eca = 1e3 * h.R * (h.celsius + 273.15) * log(cao/cai) / (2* h.FARADAY)

#Conductances
gna = gnabar * mgate**3 * hgate
gk = gkbar * ngate**4


gca_T = gcabar_T * mgate_T**2 * hgate_T
gca_N = gcabar_N * mgate_N**2 * hgate_N
gca_L = gcabar_L * mgate_L**2

#Formulas for Non-modulated Ca currents
# temp1 = 0.0853*T/2
# temp2 = v/temp1
# temp3 = 0
# if abs(temp2) < 1.10e-4:
#     temp3 = 1-temp2/2
# else:
#     temp3 = temp2/(exp(temp2)-1)
#
# dvf = 0.001/(0.001+cai)*temp1*temp3(1-cai/cao*exp(temp2))
# I_non_mod_T = -gca_T*v*((1-cai/cao*exp(2*h.FARADAY*v/k*T))/(1-exp(2*h.FARADAY*v/k*T)))
# I_non_mod_L = -gca_L*v*((1-cai/cao*exp(2*h.FARADAY*v/k*T))/(1-exp(2*h.FARADAY*v/k*T)))
# I_non_mod_N = -gca_N*v*((1-cai/cao*exp(2*h.FARADAY*v/k*T))/(1-exp(2*h.FARADAY*v/k*T)))

I_non_mod_T = -gca_T*v*((1-cai/cao*exp(2*h.FARADAY*(v-Eca)/k*T))/(1-exp(2*h.FARADAY*(v-Eca)/k*T)))
I_non_mod_L = -gca_L*v*((1-cai/cao*exp(2*h.FARADAY*(v-Eca)/k*T))/(1-exp(2*h.FARADAY*(v-Eca)/k*T)))
I_non_mod_N = -gca_N*v*((1-cai/cao*exp(2*h.FARADAY*(v-Eca)/k*T))/(1-exp(2*h.FARADAY*(v-Eca)/k*T)))

# I_non_mod_T = -gca_T*dvf
# I_non_mod_L = -gca_L*dvf
# I_non_mod_N = -gca_N*dvf

#MultiCompartmentReaction Objects that track currents of ions
na_curent = rxd.MultiCompartmentReaction(nai, nao, gna*(v-Ena), mass_action = False,
                                        membrane = mem, membrane_flux = True)

k_current = rxd.MultiCompartmentReaction(ki, ko, gk*(v - Ek), mass_action = False,
                                        membrane = mem, membrane_flux = True)

ca_T_current = rxd.MultiCompartmentReaction(cai, cao, I_non_mod_T, mass_action = False,
                                            membrane = mem, membrane_flux = True)

ca_N_current = rxd.MultiCompartmentReaction(cai, cao, I_non_mod_N, mass_action = False,
                                            membrane = mem, membrane_flux = True)

ca_L_current = rxd.MultiCompartmentReaction(cai, cao, I_non_mod_L, mass_action = False,
                                            membrane = mem, membrane_flux = True)

leak_current = rxd.MultiCompartmentReaction(xi, xo, gl*(v - el), mass_action = False,
                                            membrane = mem, membrane_flux = True)

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

fig = pyplot.figure()
pyplot.plot(tvec, vvec, label = 'rxd')
pyplot.xlabel('t (ms)')
pyplot.ylabel('V$_m$ (mV)')
pyplot.axis([0, 100, -100, 100])

fig = pyplot.figure()
pyplot.plot(tvec, hvec, '-b', label = 'h')
pyplot.plot(tvec, mvec, '-r', label = 'm')
pyplot.plot(tvec, nvec, '-g', label = 'n')
pyplot.xlabel('t (ms)')
pyplot.ylabel('state')
pyplot.axis([0, 100, 0, 1])
pyplot.legend(frameon = False)

fig = pyplot.figure()
pyplot.plot(tvec, kvec.as_numpy(), '-b', label = 'k')
pyplot.plot(tvec, navec.as_numpy(), '-r', label = 'na')
pyplot.xlabel('t (ms)')
pyplot.ylabel('current (mA/cm$^2$)')
pyplot.axis([0, 100, -2, 2])
pyplot.legend(frameon = False)

fig = pyplot.figure()
pyplot.plot(tvec, hT_vec, '-b', label = 'h_T')
pyplot.plot(tvec, mT_vec, '-r', label = 'm_T')
pyplot.xlabel('t (ms)')
pyplot.ylabel('state')
pyplot.axis([0, 100, 0, 1])
pyplot.legend(frameon = False)

fig = pyplot.figure()
pyplot.plot(tvec, hN_vec, '-b', label = 'h_N')
pyplot.plot(tvec, mN_vec, '-r', label = 'm_N')
pyplot.xlabel('t (ms)')
pyplot.ylabel('state')
pyplot.axis([0, 100, 0, 1])
pyplot.legend(frameon = False)

fig = pyplot.figure()
pyplot.plot(tvec, mL_vec, '-r', label = 'm_L')
pyplot.xlabel('t (ms)')
pyplot.ylabel('state')
pyplot.axis([0, 100, 0, 1])
pyplot.legend(frameon = False)