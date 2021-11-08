#python conversion of makecells_Axons.hoc
#started 11/3/21


from neuron import h
h.xopen("makecells_Axons.hoc")


''' 

the following is the file NOT converted into python format:



{load_file("nrngui.hoc")}
{load_file("import3d.hoc")}

//////////////////////////////////////////////////////////////////
// Template for making a cell
begintemplate Cell
   
public soma, axon, dend, apic, nclist, dend_5, dend_6, dend_7  
create soma, axon, dend, apic, dend_5, dend_6, dend_7
public all, somatic, axonal, basal, apical, dend_5, dend_6, dend_7
objref all, somatic, axonal, basal, apical, dendritic_5, dendritic_6, dendritic_7


objectvar nclist

proc init() {
    all = new SectionList()
	dend_5 all.append()
    dend_6 all.append()
    dend_7 all.append()
    axon all.append()

  axonal = new SectionList()
    axon axonal.append()

  dendritic_5 = new SectionList()
    dend_5 dendritic_5.append()

  dendritic_6 = new SectionList()
    dend_6 dendritic_6.append()

  dendritic_7 = new SectionList()
    dend_7 dendritic_7.append()
	
    somatic = new SectionList()
    basal = new SectionList()
	apical = new SectionList()

    nclist = new List()
}
endtemplate Cell

//////////////////////////////////////////////////////////////////
// Function for making a cell
obfunc mkcell() { localobj import, morph, cell
    cell = new Cell()
    morph = new Import3d_SWC_read()
    morph.input($s1)
    import = new Import3d_GUI(morph, 0)
    execute("forall delete_section()", cell)
    import.instantiate(cell)
    return cell
}

//////////////////////////////////////////////////////////////////
// Initiate parameters
objref CBC,GC,AC // list of cells
objref m,cell,strobj
strdef c,stmp,d, celltype

// Ganglion Cell Types
celltype = "D1"

cnt = 0
NUMCELLS = 1
strobj = new StringFunctions()
GC = new List()
CBC = new List()
AC = new List()
objref vc[NUMCELLS]
objref ic[NUMCELLS]

// Open list of cell names and types
ropen("CellTypes_Axon.txt")

// Create 2D array to store cell names and their indices
m = new Matrix(NUMCELLS,2)

// Loop through array and create a cell for each
for(i=0;i<NUMCELLS;i=i+1) {
	// Get morphology filename for cell
	a = fscan()
	sprint(c,"%s%d%s","morphology/",a,".swc")
	
	// Get cell type
	getstr(stmp)

	// Add cell to list and add biophysics depending on type
	// GANGLION CELL
	if (strobj.substr(stmp,"GC") > (-1)) {			
		// Make Cell
		cell = mkcell(c)
		sprint(d,"%s%d%s","Cell[",cnt,"]")
		
		print "GC: ", a, " Cell: ", cnt
		i=0
		forsec d {
			//print secname()
			/// Biophysics for the dendrites
			
			
			if (issection("Cell[0].dend[.*")) {

				insert pas
				e_pas=-60
				g_pas=.00005
				Ra=110
				Cell[0].dend[i].nseg = 5
				
				if (strobj.substr(celltype,"A2i") > (-1)) {	
					//print secname()
				
					insert spike
					gnabar_spike = 0.1
					gkbar_spike = 0.05
					gabar_spike = 3*gkbar_spike
					gcabar_spike = 0.05
					gkcbar_spike = 0.004*gkbar_spike

					insert Ih
					ghbar_Ih = 0
						
					insert IT
					gTbar_IT = 0
				}
				
				if (strobj.substr(celltype,"D1") > (-1)) {	
					print secname()
					insert spike
					gnabar_spike = 0.08
					gkbar_spike = 0.08
					gabar_spike = 3*gkbar_spike
					gcabar_spike = 0.01
					gkcbar_spike = 0.004*gkbar_spike

					insert Ih
					ghbar_Ih = 3e-5
						
					insert IT
					gTbar_IT = 0.001
				}
			i=i+1			
			}


			/// Biophysics for the axon hillock (AH)
			if (issection("Cell[0].dend_5.*")) {
			
				//print secname()
				insert pas
				e_pas=-60
				g_pas=.00005
				Ra=110
				
				Cell[0].axon.nseg = 99
				//print nseg

				insert spike
				gnabar_spike = 0.8
				gkbar_spike = 0.6
				gabar_spike = 3*gkbar_spike
				gcabar_spike = 0
				gkcbar_spike = 0

				insert Ih
				ghbar_Ih = 0
					
				insert IT
				gTbar_IT = 0

			}

			/// Biophysics for the sodium channel band (SOCB)
			if (issection("Cell[0].dend_6.*")) {
			
				//print secname()
				insert pas
				e_pas=-60
				g_pas=.00005
				Ra=110
				
				Cell[0].axon.nseg = 21
				print nseg

				insert spike
				gnabar_spike = 2.4
				gkbar_spike = 0.9
				gabar_spike = 3*gkbar_spike
				gcabar_spike = 0
				gkcbar_spike = 0

				insert Ih
				ghbar_Ih = 0
					
				insert IT
				gTbar_IT = 0
			}

			/// Biophysics for the narrow segment (NS)
			if (issection("Cell[0].dend_7.*")) {
			
				//print secname()
				insert pas
				e_pas=-60
				g_pas=.00005
				Ra=110
				
				Cell[0].axon.nseg = 21
				print nseg

				insert spike
				gnabar_spike = 0.9
				gkbar_spike = 0.7
				gabar_spike = 3*gkbar_spike
				gcabar_spike = 0
				gkcbar_spike = 0

				insert Ih
				ghbar_Ih = 0
					
				insert IT
				gTbar_IT = 0

			}

			/// Biophysics for the distal axon
			if (issection("Cell[0].ax.*")) {
			
				//print secname()
				insert pas
				e_pas=-60
				g_pas=.00005
				Ra=110
				
				Cell[0].axon.nseg = 40
				print nseg

				insert spike
				gnabar_spike = 0.8
				gkbar_spike = 0.6
				gabar_spike = 3*gkbar_spike
				gcabar_spike = 0
				gkcbar_spike = 0

				insert Ih
				ghbar_Ih = 0
					
				insert IT
				gTbar_IT = 0

			}
			/// Biophysics for the soma
			if (issection("Cell[0].s.*")) {
				
				insert pas
				e_pas=-60
				g_pas=.00005
				Ra=110				
				
				if (strobj.substr(celltype,"A2i") > (-1)) {	
					//print secname()
					insert spike
					gnabar_spike = 0.35
					gkbar_spike = 0.12
					gabar_spike = 3*gkbar_spike
					gcabar_spike = 0.137
					gkcbar_spike = 0.004*gkbar_spike

					insert Ih
					ghbar_Ih = 0
						
					insert IT
					gTbar_IT = 0.004
				}
					
				if (strobj.substr(celltype,"D1") > (-1)) {	
					//print secname()
					insert spike
					gnabar_spike = 0.2
					gkbar_spike = 0.211
					gabar_spike = 3*gkbar_spike
					gcabar_spike = 0.013
					gkcbar_spike = 0.004*gkbar_spike

					insert Ih
					ghbar_Ih = 0.0001
						
					insert IT
					gTbar_IT = 0.0024
				}
				
			}
				celsius = 22
				ena=35
				ek=-70
				insert cad
				depth_cad = 3 // (micron)
				taur_cad = 10 // (ms)
		}
		
		cnt = cnt + 1
		
/*				
		// ICLAMP FOR TESTING ION CURRENTS //////////
		access cell.soma
		ic[cnt] = new IClamp(0.5)
		ic[cnt].dur=10
		ic[cnt].amp=0.32
		ic[cnt].del=100
		cnt = cnt + 1
		////////////////////////////////////////////
*/		
		// Add to list
		GC.append(cell)
		
		// Print cell name and index of list array
		m.x[cnt-1][0] = cnt-1
		m.x[cnt-1][1] = a
	} 
}

'''