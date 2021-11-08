'''
Created on 11/5/2021
Ver1: AS
1.Load hoc file to make same cell
2.Give current clamp onto the cell
3.Record the voltage
​
Need 'makecells_Axon.hoc' file along with
1. Mod files
2. Morphology files
'''

from neuron import h
from neuron.units import mV, ms, nM, um
import matplotlib.pyplot as plt
import plotly

#Create cell
# h.load_file('init_Axon_TimeCourse.hoc')
h.load_file('makecells_Axons.hoc')
#axon is distal axon, dend_5 is AH, dend_6 is SOCB, dend_7 is NS

axon = h.Cell[0].dend_6[-1] #Access the section defined in hoc file
h.psection(sec = axon) #See biophysics options in the section you want to see
gna = axon.gnabar_spike #Access the density
axon.gnabar_spike = 1 #Change the value of density
axon.gkbar_spike = 1 #Test different values

#Current Clamp
ic = h.IClamp(axon(1))
ic.delay = 2*ms
ic.dur = 0.1*ms
ic.amp = 1

#Record
t = h.Vector().record(h._ref_t) #Time Vector
v = h.Vector().record(axon(1)._ref_v) #Voltage Vector

h.finitialize(-65*mV)
h.continuerun(10*ms)

plt.plot(t,v)
plt.show()