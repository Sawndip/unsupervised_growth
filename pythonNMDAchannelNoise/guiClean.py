"""

GUI interface for binary chain growth model
"""

import Tkinter as tk
import ttk
import reading as read
import os
import matplotlib
import numpy as np
from math import *
import sys
matplotlib.use("TkAgg")

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure

SIDE = 100
FONT_SIZE = 5
SUPERSYNAPSE_THRESHOLD = 0.015
ACTIVATION_THRESHOLD = 0.0005
TRIAL_DURATION = 1000

filename_weights = "/home/eugene/Output/weights.bin"
filename_spike_times_soma = "/home/eugene/Output/spike_times_soma.bin"
filename_spike_times_dend = "/home/eugene/Output/spike_times_dend.bin"

filename_num_synapses = "/home/eugene/Output/num_synapses.bin"
filename_spike_times_I = "/home/eugene/Output/spike_times_interneuron.bin"
filename_RA_RA_super = "/home/eugene/Output/RA_RA_super_connections.bin"

class My_gui(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        
        tk.Tk.wm_title(self, "Chain growth")        
        
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)        
        container.rowconfigure(0, weight=1)
        container.columnconfigure(0, weight=1)
        
        self.frames = {}
        
        self.pool = Updator()
        self.pIsUpdated = False
        
        for page in (PageOne, PageTwo, PageThree, PageFour, PageFive, PageSix):
            frame = page(container, self, self.pool)
            self.frames[page] = frame
            frame.grid(row=0, column=0, sticky="nswe")
        
        
        self.show_page(PageOne)
        self.update()
    def poolIsUpdated(self):
        return self.pIsUpdated
           
    def update(self):
        #print "Pool is updated!"
        if self.pool.needs_update():
            self.pool.update()
            self.pIsUpdated = True
        else:
            self.pIsUpdated = False
            
        self.after(5000, self.update)
        
    def show_page(self, cont):
        frame = self.frames[cont]        
        frame.tkraise()
    
    def quit(self):
        self.destroy()
        sys.exit()
        

class PageOne(tk.Frame):
    """
    Page with synaptic weight distribution, its mean and standard deviation
    """
    def __init__(self, parent, controller, pool):
        tk.Frame.__init__(self, parent)
         
        button1 = ttk.Button(self, text="Synapses", 
                             command=lambda: controller.show_page(PageTwo))
    
        button2 = ttk.Button(self, text="All spikes HVC(RA)", 
                             command=lambda: controller.show_page(PageThree))
        
        button3 = ttk.Button(self, text="Spikes of strongly connected HVC(RA)", 
                             command=lambda: controller.show_page(PageFour))                
        
        button4 = ttk.Button(self, text="Spikes HVC(I)", 
                             command=lambda: controller.show_page(PageFive))
                             
        button5 = ttk.Button(self, text="Maturation", 
                             command=lambda: controller.show_page(PageSix))
       
        button1.grid(row=0, column=0)                                
        button2.grid(row=1, column=0)
        button3.grid(row=2, column=0)
        button4.grid(row=3, column=0)
        button5.grid(row=4, column=0)
        
        self.pool = pool
        
        f_dist = Figure(figsize=(6,6), dpi=100)
        f_weight_statistics = Figure(figsize=(6,6), dpi=100)
        
        self.ax_dist = f_dist.add_subplot(111)
        
        self.ax1_weight_statistics = f_weight_statistics.add_subplot(211)
        self.ax2_weight_statistics = f_weight_statistics.add_subplot(212)
        
        self.canvas_dist = FigureCanvasTkAgg(f_dist, self)
        self.canvas_dist.show()
        self.canvas_weight_statistics = FigureCanvasTkAgg(f_weight_statistics, self)
        self.canvas_weight_statistics.show()
        
        toolbar_dist = NavigationToolbar2TkAgg(self.canvas_dist, self)
        toolbar_dist.update()
        toolbar_dist.grid(row=0, column=1, sticky="ew")

        toolbar_weight_statistics = NavigationToolbar2TkAgg(self.canvas_weight_statistics, self)
        toolbar_weight_statistics.update()
        toolbar_weight_statistics.grid(row=0, column=21, sticky="ew")

        
        #toolbar.pack(side=tk.TOP)
        #self.canvas_dist.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=False)
        self.canvas_dist.get_tk_widget()
        self.canvas_weight_statistics.get_tk_widget()
        #self.canvas_dist._tkcanvas.pack()
        self.canvas_dist._tkcanvas.grid(row=1, rowspan = 20, column=1, columnspan = 20, sticky="nsew")        
        self.canvas_weight_statistics._tkcanvas.grid(row=1, rowspan = 20, column=21, sticky="nsew")        
        
                  
         
        #self.canvas_dist._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=False)        
        
        self.update(controller)
        
    def update(self, controller):
        """
        Refreshes page of needed
        """
        if controller.poolIsUpdated():
            #print "Distribution page update"            
            self.pool.draw_weight_hist(self.ax_dist)
            self.canvas_dist.draw()
            self.pool.draw_weight_statistics(self.ax1_weight_statistics, self.ax2_weight_statistics)
            self.canvas_weight_statistics.draw()            
        controller.after(5000, self.update, controller)
            
class PageTwo(tk.Frame):
    """
    
    Page with number of active and supersynapses
    """
    def __init__(self, parent, controller, pool):
        tk.Frame.__init__(self, parent)
        
        self.pool = pool        
        
       
        button1 = ttk.Button(self, text="Distribution", 
                             command=lambda: controller.show_page(PageOne))
        
        button2 = ttk.Button(self, text="All spikes HVC(RA)", 
                             command=lambda: controller.show_page(PageThree))
        
        button3 = ttk.Button(self, text="Spikes of strongly connected HVC(RA)", 
                             command=lambda: controller.show_page(PageFour))
        
        button4 = ttk.Button(self, text="Spikes HVC(I)", 
                             command=lambda: controller.show_page(PageFive))
                             
        button5 = ttk.Button(self, text="Maturation", 
                             command=lambda: controller.show_page(PageSix))
                
        
        button1.grid(row=0, column=0)
        button2.grid(row=1, column=0)
        button3.grid(row=2, column=0)
        button4.grid(row=3, column=0)
        button5.grid(row=4, column=0)
        
        f_active = Figure(figsize=(6,6), dpi=100)
        self.ax_active = f_active.add_subplot(111)
        self.canvas_active = FigureCanvasTkAgg(f_active, self)
        self.canvas_active.show()        
        
        f_super = Figure(figsize=(6,6), dpi=100)
        self.ax_super = f_super.add_subplot(111)
        self.canvas_super = FigureCanvasTkAgg(f_super, self)
        self.canvas_super.show()
        
        toolbar_active = NavigationToolbar2TkAgg(self.canvas_active, self)
        toolbar_active.update()
        toolbar_active.grid(row=0, column=1, sticky="ew")        

        toolbar_super = NavigationToolbar2TkAgg(self.canvas_super, self)
        toolbar_super.update()
        toolbar_super.grid(row=0, column=21, sticky="ew")        

        
        self.canvas_active.get_tk_widget()
        self.canvas_active._tkcanvas.grid(row=1, column=1, rowspan = 20,
                                            columnspan=20, sticky="nsew")

        self.canvas_super.get_tk_widget()
        self.canvas_super._tkcanvas.grid(row=1, column=21, rowspan = 20,
                                             sticky="nsew")



        self.update(controller)
        
    def update(self, controller):
        """
        Refreshes page of needed
        """
        if controller.poolIsUpdated():
            self.pool.draw_active_synapse_dynamics(self.ax_active)
            self.canvas_active.draw()
            self.pool.draw_supersynapse_dynamics(self.ax_super)
            self.canvas_super.draw()
        controller.after(5000, self.update, controller)
        

class PageThree(tk.Frame):
    """
    Page with raster plots of HVC(RA) spikes in trial
    """
    def __init__(self, parent, controller, pool):
        tk.Frame.__init__(self, parent)
        
        button1 = ttk.Button(self, text="Distribution", 
                             command=lambda: controller.show_page(PageOne))
        
        button2 = ttk.Button(self, text="Synapses", 
                             command=lambda: controller.show_page(PageTwo))
        
        button3 = ttk.Button(self, text="Spikes of strongly connected HVC(RA)", 
                             command=lambda: controller.show_page(PageFour))
        
        button4 = ttk.Button(self, text="Spikes HVC(I)", 
                             command=lambda: controller.show_page(PageFive))
    
        button5 = ttk.Button(self, text="Maturation", 
                             command=lambda: controller.show_page(PageSix))
                
        
        button1.grid(row=0, column=0)
        button2.grid(row=1, column=0)
        button3.grid(row=2, column=0)
        button4.grid(row=3, column=0)
        button5.grid(row=4, column=0)
        
        self.pool = pool
        
        
        f_soma = Figure(figsize=(6,6), dpi=100)
        f_dend = Figure(figsize=(6,6), dpi=100)
        
        self.ax_soma = f_soma.add_subplot(111)
        self.ax_dend = f_dend.add_subplot(111)        
        
        self.canvas_soma = FigureCanvasTkAgg(f_soma, self)
        self.canvas_soma.show()
        self.canvas_dend = FigureCanvasTkAgg(f_dend, self)
        self.canvas_dend.show()
        
        toolbar_soma = NavigationToolbar2TkAgg(self.canvas_soma, self)
        toolbar_soma.update()
        toolbar_soma.grid(row=0, column=1, sticky="ew")

        toolbar_dend = NavigationToolbar2TkAgg(self.canvas_dend, self)
        toolbar_dend.update()
        toolbar_dend.grid(row=0, column=21, sticky="ew")

        
        #toolbar.pack(side=tk.TOP)
        #self.canvas_dist.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=False)
        self.canvas_soma.get_tk_widget()
        self.canvas_dend.get_tk_widget()
        #self.canvas_dist._tkcanvas.pack()
        self.canvas_soma._tkcanvas.grid(row=1, rowspan = 20, column=1, columnspan = 20, sticky="nsew")        
        self.canvas_dend._tkcanvas.grid(row=1, rowspan = 20, column=21, sticky="nsew")        
        
                  
         
        #self.canvas_dist._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=False)        
        
        self.update(controller)
        
    def update(self, controller):
        """
        Refreshes page if needed
        """
        if controller.poolIsUpdated():
            #print "Distribution page update"            
            self.pool.draw_raster_all_soma(self.ax_soma)
            self.canvas_soma.draw()
            self.pool.draw_raster_all_dend(self.ax_dend)
            self.canvas_dend.draw()            
        controller.after(5000, self.update, controller)

class PageFour(tk.Frame):
    """
    Page with raster plots of strongly connected HVC(RA) spikes in trial
    """
    def __init__(self, parent, controller, pool):
        tk.Frame.__init__(self, parent)
        
        button1 = ttk.Button(self, text="Distribution", 
                             command=lambda: controller.show_page(PageOne))
        
        button2 = ttk.Button(self, text="Synapses", 
                             command=lambda: controller.show_page(PageTwo))
        
        button3 = ttk.Button(self, text="All spikes HVC(RA)", 
                             command=lambda: controller.show_page(PageThree))
        
        button4 = ttk.Button(self, text="Spikes HVC(I)", 
                             command=lambda: controller.show_page(PageFive))
    
        button5 = ttk.Button(self, text="Maturation", 
                             command=lambda: controller.show_page(PageSix))
                
        
        button1.grid(row=0, column=0)
        button2.grid(row=1, column=0)
        button3.grid(row=2, column=0)
        button4.grid(row=3, column=0)
        button5.grid(row=4, column=0)
        
        self.pool = pool
        
        f_strongly_connected_soma = Figure(figsize=(6,6), dpi=100)
        f_strongly_connected_dend = Figure(figsize=(6,6), dpi=100)
        
        self.ax_strongly_connected_soma = f_strongly_connected_soma.add_subplot(111)
        self.ax_strongly_connected_dend = f_strongly_connected_dend.add_subplot(111)        
        
        self.canvas_strongly_connected_soma = FigureCanvasTkAgg(f_strongly_connected_soma, self)
        self.canvas_strongly_connected_soma.show()
        self.canvas_strongly_connected_dend = FigureCanvasTkAgg(f_strongly_connected_dend, self)
        self.canvas_strongly_connected_dend.show()
        
        toolbar_strongly_connected_soma = NavigationToolbar2TkAgg(self.canvas_strongly_connected_soma, self)
        toolbar_strongly_connected_soma.update()
        toolbar_strongly_connected_soma.grid(row=0, column=1, sticky="ew")

        toolbar_strongly_connected_dend = NavigationToolbar2TkAgg(self.canvas_strongly_connected_dend, self)
        toolbar_strongly_connected_dend.update()
        toolbar_strongly_connected_dend.grid(row=0, column=21, sticky="ew")

        
        #toolbar.pack(side=tk.TOP)
        #self.canvas_dist.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=False)
        self.canvas_strongly_connected_soma.get_tk_widget()
        self.canvas_strongly_connected_dend.get_tk_widget()
        #self.canvas_dist._tkcanvas.pack()
        self.canvas_strongly_connected_soma._tkcanvas.grid(row=1, rowspan = 20, column=1, columnspan = 20, sticky="nsew")        
        self.canvas_strongly_connected_dend._tkcanvas.grid(row=1, rowspan = 20, column=21, sticky="nsew")        
        
                  
         
        #self.canvas_dist._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=False)        
        
        self.update(controller)
        
    def update(self, controller):
        """
        Refreshes page if needed
        """
        if controller.poolIsUpdated():
            #print "Distribution page update"            
            self.pool.draw_raster_strongly_connected_soma(self.ax_strongly_connected_soma)
            self.canvas_strongly_connected_soma.draw()
            self.pool.draw_raster_strongly_connected_dend(self.ax_strongly_connected_dend)
            self.canvas_strongly_connected_dend.draw()            
        controller.after(5000, self.update, controller)          

class PageFive(tk.Frame):
    """
    Page with interneuron spikes
    """
    def __init__(self, parent, controller, pool):
        tk.Frame.__init__(self, parent)
        
         
        button1 = ttk.Button(self, text="Distribution", 
                             command=lambda: controller.show_page(PageOne))
        
        button2 = ttk.Button(self, text="Synapses", 
                             command=lambda: controller.show_page(PageTwo))
        
        button3 = ttk.Button(self, text="All spikes HVC(RA)", 
                             command=lambda: controller.show_page(PageThree))

        button4 = ttk.Button(self, text="Spikes of strongly connected HVC(RA)", 
                             command=lambda: controller.show_page(PageFour))
    
        button5 = ttk.Button(self, text="Maturation", 
                             command=lambda: controller.show_page(PageSix))
                
        
        button1.grid(row=0, column=0)
        button2.grid(row=1, column=0)
        button3.grid(row=2, column=0)
        button4.grid(row=3, column=0)
        button5.grid(row=4, column=0)
        
        self.pool = pool
        
        f_interneuron = Figure(figsize=(6,6), dpi=100)
        
        self.ax_interneuron = f_interneuron.add_subplot(111)
        
        self.canvas_interneuron = FigureCanvasTkAgg(f_interneuron, self)
        self.canvas_interneuron.show()
        
        toolbar_interneuron = NavigationToolbar2TkAgg(self.canvas_interneuron, self)
        toolbar_interneuron.update()
        toolbar_interneuron.grid(row=0, column=1, sticky="ew")

        self.canvas_interneuron.get_tk_widget()
        self.canvas_interneuron._tkcanvas.grid(row=1, rowspan = 20, column=1, columnspan = 20, sticky="nsew")        
        
        self.update(controller)
        
    def update(self, controller):
        """
        Refreshes page if needed
        """
        if controller.poolIsUpdated():
            #print "Distribution page update"            
            self.pool.draw_raster_all_I(self.ax_interneuron)
            self.canvas_interneuron.draw()
                       
        controller.after(5000, self.update, controller)

class PageSix(tk.Frame):
    """
    Page with maturation info
    """
    def __init__(self, parent, controller, pool):
        tk.Frame.__init__(self, parent)
        
         
        button1 = ttk.Button(self, text="Distribution", 
                             command=lambda: controller.show_page(PageOne))
        
        button2 = ttk.Button(self, text="Synapses", 
                             command=lambda: controller.show_page(PageTwo))
        
        button3 = ttk.Button(self, text="All spikes HVC(RA)", 
                             command=lambda: controller.show_page(PageThree))
                             
        button4 = ttk.Button(self, text="Spikes of strongly connected HVC(RA)", 
                             command=lambda: controller.show_page(PageFour))
    
        button5 = ttk.Button(self, text="Spikes HVC(I)", 
                             command=lambda: controller.show_page(PageFive))
                
        
        button1.grid(row=0, column=0)
        button2.grid(row=1, column=0)
        button3.grid(row=2, column=0)
        button4.grid(row=3, column=0)
        button5.grid(row=4, column=0)
        
        self.pool = pool
        
        # create text window        
        tx = tk.Text(self,font=('times',12),width=50,height=15,wrap='word')
        
        tx.grid(row=1, rowspan = 20, column=1, columnspan = 20, sticky="nsew")
        tx.insert(1.0, "Mature neurons:\n")
        
        
        self.update(controller)
        
    def update(self, controller):
        """
        Refreshes page if needed
        """
        if controller.poolIsUpdated():
            pass
            #print "Distribution page update"            
                           
        controller.after(5000, self.update, controller)

class Updator:
    """
    class Update updates all graphs
    """
    def __init__(self):
        self.file_weights = filename_weights # file with all synaptic weights
        self.file_spike_times_soma = filename_spike_times_soma # file with somatic spike times of all HVC(RA) neurons
        self.file_spike_times_dend = filename_spike_times_dend # file with dendritic spike times of all HVC(RA) neurons        
        self.file_spike_times_I = filename_spike_times_I # file with spike times of all HVC(I) neurons
        self.file_super = filename_RA_RA_super # file with supersynaptic connections
        self.file_num_synapses = filename_num_synapses # file with number of active and supersynapses                
                
        self.std_weights = [] # standard variation of synaptic weight distribution
        self.mean_weights = [] # mean of synaptic weight distribution
        self.time = [] # simulation time in trials
        self.num_active_synapses = [] # number of active synapses
        self.num_supersynapses = [] # number of supersynapses

        self.mod_time_wgh = 0 # modification time of file with synaptic weights
        self.mod_time_spikes_RA = 0 # modification time of file with RA spikes
        self.mod_time_spikes_I = 0 # modification time of file with I spikes
        self.mod_time_super = 0 # modification time of file with sypersynaptic connections
                
        
        
    def needs_update(self):
        """
        Returns true if all files were modified
        """
        return (self.mod_time_super!=os.path.getmtime(self.file_super) and \
                self.mod_time_wgh!=os.path.getmtime(self.file_weights) and \
                self.mod_time_spikes_RA!=os.path.getmtime(self.file_spike_times_soma) and \
                self.mod_time_spikes_I!=os.path.getmtime(self.file_spike_times_I))
    
    def update(self):
        """
        Update all graphs
        """
        self.update_spike_times()
        self.update_weights_info()
        self.update_num_synapses()
    
    def update_spike_times(self):
        """
        Updates spike times of all neurons in the network
        """
        # update modification time of file with spike times
        self.mod_time_spikes_RA = os.path.getmtime(self.file_spike_times_soma)
        
        # get spike times of all HVC(RA) and HVC(I) neurons in the network
        trial_number, self.spike_times_dend, self.fired_dend_id = self.get_spike_times(self.file_spike_times_dend)
        trial_number, self.spike_times_soma, self.fired_soma_id = self.get_spike_times(self.file_spike_times_soma)
        trial_number, self.spike_times_I, self.fired_I_id = self.get_spike_times(self.file_spike_times_I)
        
        #print self.spike_times_I        
        
        # update time
        self.time.append(trial_number)
        
        # get ordered spike times of strongly connected neurons
        self.get_spike_times_strongly_connected()
    
    def update_num_synapses(self):
        """
        Updates number of active and supersynapses
        """
        # update number of active and super synapses
        self.trial_number_num_synapses, self.num_active_synapses, self.num_supersynapses = read.read_num_synapses(self.file_num_synapses)
          
        
    def update_weights_info(self):
        """
        Updates number of synapses, weight histogram and statistical measures of synaptic weight distribution
        """
        (self.N_RA, weights) = read.read_weights(self.file_weights)
        self.mod_time_wgh = os.path.getmtime(self.file_weights)
        
        self.weights = [item for sublist in weights for item in sublist]
        
        for i in range(len(self.weights)):
            if self.weights[i] < sys.float_info.epsilon:
                self.weights[i] = 0
        
        self.hist, self.bin_edges = np.histogram(self.weights, bins=100)
        
        # update mean and std dependences on time        
        sigma = std(self.weights)
        mu = mean(self.weights)
        
        self.mean_weights.append(mu)
        self.std_weights.append(sigma)        
        
        
        
    def get_spike_times(self, filename):
        """
        Reads spike times of neurons and checks if they are in range
        """
        # read spikes from file
        trial_number, simulation_time, spike_times, fired_id = read.read_time_info(filename)
        
    
        # flatten lists
        spike_times = [item for sublist in spike_times for item in sublist]
        fired_id = [item for sublist in fired_id for item in sublist]
        
        # check that spike times are in range     
        if len(spike_times) > 0:
            if min(spike_times) < 0:
                print "Spike time read from file {0} is negative! Min somatic spike = {1}".format(filename, min(spike_times))
            if max(spike_times) > TRIAL_DURATION:
                print "Spike time read from file {0} exceed trial duration! Max somatic spike = {1}".format(filename, max(spike_times))
        
        return trial_number, spike_times, fired_id
    
    def get_spike_times_strongly_connected(self):
        """
        Returns spike times ordered sequence of strongly connected neurons
        """
        # read supersynapses from file
        (unused, RA_super_targets, unused) = read.read_connections(self.file_super)

        # extract strongly connected neurons
        superTargets = set([element for sublist in RA_super_targets for element in sublist])
        superSources = set([i for i in xrange(len(RA_super_targets)) if len(RA_super_targets[i]) > 0])
        stronglyConnected = superTargets | superSources

        self.fired_dend_id_strongly_connected = [i for i in self.fired_dend_id if i in stronglyConnected]
        self.spike_times_dend_strongly_connected = [t for ind, t in enumerate(self.spike_times_dend) if self.fired_dend_id[ind] in stronglyConnected]
        
        self.fired_soma_id_strongly_connected = [i for i in self.fired_soma_id if i in stronglyConnected]
        self.spike_times_soma_strongly_connected = [t for ind, t in enumerate(self.spike_times_soma) if self.fired_soma_id[ind] in stronglyConnected]
        
        self.order_strongly_connected_fired()
        
   
                         
    def order_strongly_connected_fired(self):
        """
        Order strongly connected fired HVC(RA) neurons
        """
        if len(self.spike_times_dend_strongly_connected) > 0:
            # start spikes from the first fired strongly connected neuron
            first_spike = min(self.spike_times_dend_strongly_connected)
            self.spike_times_strongly_dend_connected = map(lambda x: x - first_spike, self.spike_times_dend_strongly_connected)
            # sort spikes in ascending order
            self.spike_times_dend_strongly_connected, self.fired_dend_id_strongly_connected = \
                            zip(*sorted(zip(self.spike_times_dend_strongly_connected, self.fired_dend_id_strongly_connected)))
            # get neuron ids in ascending order
            self.fired_dend_id_strongly_connected_ascending = remap_ID(self.fired_dend_id_strongly_connected)
        
        if len(self.spike_times_soma_strongly_connected) > 0:
            # start spikes from the first fired strongly connected neuron
            first_spike = min(self.spike_times_soma_strongly_connected)
            self.spike_times_strongly_soma_connected = map(lambda x: x - first_spike, self.spike_times_soma_strongly_connected)
            # sort spikes in ascending order
            self.spike_times_soma_strongly_connected, self.fired_soma_id_strongly_connected = \
                            zip(*sorted(zip(self.spike_times_soma_strongly_connected, self.fired_soma_id_strongly_connected)))
            # get neuron ids in ascending order
            self.fired_soma_id_strongly_connected_ascending = remap_ID(self.fired_soma_id_strongly_connected)
    
    def draw_raster_all_dend(self, axes):
        """
        Draw raster plot of dendritic spike of all HVC(RA) neurons
        """
        axes.clear()
        
        for i in range(len(self.spike_times_dend)):
            axes.vlines(self.spike_times_dend[i], self.fired_dend_id[i]-0.5, self.fired_dend_id[i]+0.5)
        
        axes.set_ylabel("Neuron ID")
        axes.set_xlabel("spike time dendrite (ms)")
       
        if len(self.fired_dend_id) > 1:
            axes.set_ylim(-0.5, max(self.fired_dend_id) + 0.5)
            axes.set_xlim(0, TRIAL_DURATION)
        
        axes.set_title("Dendritic spikes of HVC(RA) neurons")
    
    def draw_raster_all_soma(self, axes):
        """
        Draw raster plot of somatic spike of all HVC(RA) neurons
        """
        axes.clear()
        
        for i in range(len(self.spike_times_soma)):
            axes.vlines(self.spike_times_soma[i], self.fired_soma_id[i]-0.5, self.fired_soma_id[i]+0.5)
        
        axes.set_ylabel("Neuron ID")
        axes.set_xlabel("spike time soma (ms)")
       
        if len(self.fired_soma_id) > 1:
            axes.set_ylim(-0.5, max(self.fired_soma_id) + 0.5)
            axes.set_xlim(0, TRIAL_DURATION)
        
        axes.set_title("Somatic spikes of HVC(RA) neurons")
        
    
    def draw_raster_all_I(self, axes):
        """
        Draw all spikes of I neurons
        """
        #print "Refresh spike times!"
        if len(self.spike_times_I) > 0:        
            axes.clear()
            
            for i in range(len(self.spike_times_I)):
                axes.vlines(self.spike_times_I[i], self.fired_I_id[i]-0.5, self.fired_I_id[i]+0.5)
            
            axes.set_ylabel("neuron ID")
            axes.set_xlabel("spike time (ms)")
           
            axes.set_ylim(-0.5, max(self.fired_I_id)+0.5)
            axes.set_xlim(0, TRIAL_DURATION)
            axes.set_title("Spikes of HVC(I) neurons")
    

    def draw_raster_strongly_connected_dend(self, axes):
        """
        Draw ordered dendritic spikes of strongly connected neurons
        """
        #print "Refresh spike times!"
        # if there are some spikes to show
        if len(self.spike_times_dend_strongly_connected) > 0:
            axes.clear()
            
            for i in range(len(self.spike_times_dend_strongly_connected)):
                axes.vlines(self.spike_times_dend_strongly_connected[i], self.fired_dend_id_strongly_connected_ascending[i]-0.5,\
                            self.fired_dend_id_strongly_connected_ascending[i]+0.5)
        
        #ax.scatter(spike_times, random_ID)
            axes.set_yticks(self.fired_dend_id_strongly_connected_ascending)
            
            axes.set_yticklabels(self.fired_dend_id_strongly_connected)
            axes.set_xlim([-5, max(self.spike_times_dend_strongly_connected)+5])
            axes.set_ylabel("neuron ID")
            axes.set_xlabel("relative spike time (ms)")
            axes.set_title("Ordered dendritic spikes of strongly connected neurons")
    
    def draw_raster_strongly_connected_soma(self, axes):
        """
        Draw ordered somatic spikes of strongly connected neurons
        """
        #print "Refresh spike times!"
        # if there are some spikes to show
        if len(self.spike_times_soma_strongly_connected) > 0:
            axes.clear()
            
            for i in range(len(self.spike_times_soma_strongly_connected)):
                axes.vlines(self.spike_times_soma_strongly_connected[i], self.fired_soma_id_strongly_connected_ascending[i]-0.5,\
                            self.fired_soma_id_strongly_connected_ascending[i]+0.5)
        
        #ax.scatter(spike_times, random_ID)
            axes.set_yticks(self.fired_soma_id_strongly_connected_ascending)
            
            axes.set_yticklabels(self.fired_soma_id_strongly_connected)
            axes.set_xlim([-5, max(self.spike_times_soma_strongly_connected)+5])
            axes.set_ylabel("neuron ID")
            axes.set_xlabel("relative spike time (ms)")
            axes.set_title("Ordered somatic spikes of strongly connected neurons")
    
    def draw_active_synapse_dynamics(self, axes):
        """
        Draw number of active synapses
        """
        #print "Resfresh synapses!"
        #print self.active_synapses
        #print self.supersynapses        
        axes.clear()
        axes.plot(self.trial_number_num_synapses, self.num_active_synapses, 'g', label="Active")
          
        axes.set_xlabel("time (s)")
        axes.set_ylabel("number of active synapses")
    
    def draw_supersynapse_dynamics(self, axes):
        """
        Draw number of supersynapses
        """
        axes.clear()
        axes.plot(self.trial_number_num_synapses, self.num_supersynapses, 'b', label="Supersynapses")
        
        axes.set_xlabel("time (s)")
        axes.set_ylabel("number of super synapses")
        
    
    def draw_weight_hist(self, axes):
        """
        Draw synaptic weight histogram
        """
        #print "Refresh weights!"
        #print type(axes)
        width = 0.7 * (self.bin_edges[1] - self.bin_edges[0])
        center = (self.bin_edges[:-1] + self.bin_edges[1:]) / 2
        axes.clear()
        axes.bar(center, self.hist, align="center", width=width, log=1)
        axes.set_ylabel("Number of synapses")
        #ymin, ymax = axes.get_ylim()        
        #axes.set_ylim([0, ymax + 1])
        axes.set_xlabel("Synaptic weight")
        axes.ticklabel_format(style="sci", axis="x", scilimits=(0,0))      
        
    def draw_weight_statistics(self, ax1, ax2):
        """
        Draw mean and standard deviation of synaptic weight distribution
        """
        ax1.clear()
        ax1.plot(self.time, self.mean_weights)
        ax1.set_ylabel("mean synaptic weight")
        ax1.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
        ax1.get_xaxis().set_visible(False)        
        
        ax2.clear()
        ax2.plot(self.time, self.std_weights)
        ax2.set_ylabel("std of synaptic weights")
        ax2.set_xlabel("time")
        ax2.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
        ax2.ticklabel_format(style="sci", axis="x", scilimits=(0,0))
        

def remap_ID(order):
    """
    Assigns id labels ti neurons according to their spiking order
    """
    already_counted = []
    remapped = []
    ind = 0
    
    for i in order:
        if i in already_counted:
            ind_in_counted = already_counted.index(i)
            remapped.append(ind_in_counted)
        else:
            remapped.append(ind)
            already_counted.append(i)
            ind += 1
            
    return remapped

def mean(array):
    """
    Calculates mean of the numbers in array
    """
    n = len(array)
    if n < 1:
        raise ValueError("Mean requires at least one element in array")
    
    return sum(array)/n
 
def std(array):
    """
    Calculates standard deviation of the number in array
    """
    #copy = array.flatten()
    n = len(array)
    if n < 1:
        raise ValueError("Std requires at least one element in array")
    mu = mean(array)
    return sqrt(sum([(x - mu) ** 2 for x in array])/n)

app = My_gui()
app.mainloop()
