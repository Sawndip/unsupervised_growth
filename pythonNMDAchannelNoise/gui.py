import Tkinter as tk
import ttk
import reading as read
import networkx as nx
import os
import matplotlib
import numpy as np
from math import *
import sys
matplotlib.use("TkAgg")

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib.pyplot import figure

SIDE = 100
FONT_SIZE = 5
SUPERSYNAPSE_THRESHOLD = 0.010
TRIAL_DURATION = 1000

filenameRA_RA = "/home/eugene/Output/RA_RA_connections.bin"
filenameRA_I = "/home/eugene//Output/RA_I_connections.bin"
filenameI = "/home/eugene/Output/I_RA_connections.bin"
filename_xy_RA = "/home/eugene/Output/RA_xy.bin"
filename_xy_I = "/home/eugene/Output/I_xy.bin"
filename_weights = "/home/eugene/Output/weights.bin"
filename_time_dend = "/home/eugene/Output/time_info_dend.bin"
filename_time_soma = "/home/eugene/Output/time_info_soma.bin"
filenameRA_RA_super = "/home/eugene/Output/RA_RA_super_connections.bin"

class My_gui(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        
        tk.Tk.wm_title(self, "Chain growth")        
        
        container = tk.Frame(self)
        container.pack(side="top", fill="both", expand=True)        
        container.rowconfigure(0, weight=1)
        container.columnconfigure(0, weight=1)
        
        self.frames = {}
        
        self.pool = Pool()
        self.pIsUpdated = False
        
        for page in (PageOne, PageTwo, PageThree, PageFour, PageFive, PageSix, PageSeven):
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
    def __init__(self, parent, controller, pool):
    
        tk.Frame.__init__(self, parent)
        
        button1 = ttk.Button(self, text="Ra2I", 
                             command=lambda: controller.show_page(PageFour))
        button2 = ttk.Button(self, text="I2RA", 
                             command=lambda: controller.show_page(PageFive))
        button3 = ttk.Button(self, text="Distribution", command=lambda: controller.show_page(PageTwo))        
        button4 = ttk.Button(self, text="Synapses", command=lambda: controller.show_page(PageThree))
        button5 = ttk.Button(self, text="Raster", 
                             command=lambda: controller.show_page(PageSix))
                             
        button6 = ttk.Button(self, text="Ordered spikes", 
                             command=lambda: controller.show_page(PageSeven)) 
                             
        button7 = ttk.Button(self, text="Quit", command = lambda: controller.quit())
        
        button1.grid(row=0, column=0)        
        button2.grid(row=1, column=0)  	    	
        button3.grid(row=2, column=0)
        button4.grid(row=3, column=0)
        button5.grid(row=4, column=0)
        button6.grid(row=5, column=0)
        button7.grid(row=6, column=0)
        
        self.pool = pool        
        
        f = figure(figsize=(6,3), dpi=200)
        self.ax = f.add_subplot(111)
        
        
        self.canvas = FigureCanvasTkAgg(f, self)
        self.canvas.show()
        
        #self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        toolbar.update()
        toolbar.grid(row=0, column=1, sticky="ew")
        self.canvas.get_tk_widget()
        self.canvas._tkcanvas.grid(row=1, column=1, rowspan=20, columnspan=20, sticky="nsew")        
        #self.canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.update(controller)
        
    def update(self, controller):
        if controller.poolIsUpdated():
            self.pool.draw_graph(self.ax)
            self.canvas.draw()
        controller.after(5000, self.update, controller)
	
        
class PageTwo(tk.Frame):
    def __init__(self, parent, controller, pool):
        tk.Frame.__init__(self, parent)
        
        button1 = ttk.Button(self, text="Graph", 
                            command=lambda: controller.show_page(PageOne))        
       
        button2 = ttk.Button(self, text="RA2I", 
                             command=lambda: controller.show_page(PageFour))
        
        button3 = ttk.Button(self, text="I2RA", 
                             command=lambda: controller.show_page(PageFive))
         
        button4 = ttk.Button(self, text="Synapses", 
                             command=lambda: controller.show_page(PageThree))
    
        button5 = ttk.Button(self, text="Raster", 
                             command=lambda: controller.show_page(PageSix))
        
        button6 = ttk.Button(self, text="Ordered spikes", 
                             command=lambda: controller.show_page(PageSeven)) 
         
        button1.grid(row=0, column=0)                                
        button2.grid(row=1, column=0)
        button3.grid(row=2, column=0)
        button4.grid(row=3, column=0)
        button5.grid(row=4, column=0)
        button6.grid(row=5, column=0)
         
        self.pool = pool
        
        f_dist = Figure(figsize=(6,6), dpi=100)
        f_time_statistics = Figure(figsize=(6,6), dpi=100)
        
        self.ax_dist = f_dist.add_subplot(111)
        
        self.ax1_time_statistics = f_time_statistics.add_subplot(211)
        self.ax2_time_statistics = f_time_statistics.add_subplot(212)
        
        self.canvas_dist = FigureCanvasTkAgg(f_dist, self)
        self.canvas_dist.show()
        self.canvas_time_statistics = FigureCanvasTkAgg(f_time_statistics, self)
        self.canvas_time_statistics.show()
        
        toolbar_dist = NavigationToolbar2TkAgg(self.canvas_dist, self)
        toolbar_dist.update()
        toolbar_dist.grid(row=0, column=1, sticky="ew")

        toolbar_time_statistics = NavigationToolbar2TkAgg(self.canvas_time_statistics, self)
        toolbar_time_statistics.update()
        toolbar_time_statistics.grid(row=0, column=21, sticky="ew")

        
        #toolbar.pack(side=tk.TOP)
        #self.canvas_dist.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=False)
        self.canvas_dist.get_tk_widget()
        self.canvas_time_statistics.get_tk_widget()
        #self.canvas_dist._tkcanvas.pack()
        self.canvas_dist._tkcanvas.grid(row=1, rowspan = 20, column=1, columnspan = 20, sticky="nsew")        
        self.canvas_time_statistics._tkcanvas.grid(row=1, rowspan = 20, column=21, sticky="nsew")        
        
                  
         
        #self.canvas_dist._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=False)        
        
        self.update(controller)
        
    def update(self, controller):
        
        if controller.poolIsUpdated():
            #print "Distribution page update"            
            self.pool.draw_hist(self.ax_dist)
            self.canvas_dist.draw()
            self.pool.draw_time_statistics(self.ax1_time_statistics, self.ax2_time_statistics)
            self.canvas_time_statistics.draw()            
        controller.after(5000, self.update, controller)
            
class PageThree(tk.Frame):
    def __init__(self, parent, controller, pool):
        tk.Frame.__init__(self, parent)
        
        self.pool = pool        
        
        button1 = ttk.Button(self, text="Graph",
                             command=lambda: controller.show_page(PageOne))        
        
        button2 = ttk.Button(self, text="Ra2I",
                             command=lambda: controller.show_page(PageFour))        
        
        button3 = ttk.Button(self, text="I2RA", 
                             command=lambda: controller.show_page(PageFive))
       
        button4 = ttk.Button(self, text="Distribution", 
                             command=lambda: controller.show_page(PageTwo))
        
        button5 = ttk.Button(self, text="Raster", 
                             command=lambda: controller.show_page(PageSix))
        
        button6 = ttk.Button(self, text="Ordered spikes", 
                             command=lambda: controller.show_page(PageSeven))        
        
        button1.grid(row=0, column=0)
        button2.grid(row=1, column=0)
        button3.grid(row=2, column=0)
        button4.grid(row=3, column=0)
        button5.grid(row=4, column=0)
        button6.grid(row=5, column=0)
        
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
        if controller.poolIsUpdated():
            self.pool.draw_active_synapses(self.ax_active)
            self.canvas_active.draw()
            self.pool.draw_supersynapses(self.ax_super)
            self.canvas_super.draw()
        controller.after(5000, self.update, controller)
        

class PageFour(tk.Frame):
    def __init__(self, parent, controller, pool):
        tk.Frame.__init__(self, parent)
        
        self.pool = pool        
        
        button1 = ttk.Button(self, text="Graph",
                             command=lambda: controller.show_page(PageOne))        
        
        button2 = ttk.Button(self, text="I2RA", 
                             command=lambda: controller.show_page(PageFive))
        
        button3 = ttk.Button(self, text="Distribution", 
                             command=lambda: controller.show_page(PageTwo))
                             
        button4 = ttk.Button(self, text="Synapses", 
                             command=lambda: controller.show_page(PageThree))
        
        button5 = ttk.Button(self, text="Raster", 
                             command=lambda: controller.show_page(PageSix))
        
        button6 = ttk.Button(self, text="Ordered spikes", 
                             command=lambda: controller.show_page(PageSeven))        
        
        button1.grid(row=0, column=0)
        button2.grid(row=1, column=0)
        button3.grid(row=2, column=0)
        button4.grid(row=3, column=0)
        button5.grid(row=4, column=0)
        button6.grid(row=5, column=0)
        
        f = figure(figsize=(6,3), dpi=200)
        self.ax = f.add_subplot(111)
        
        
        self.canvas = FigureCanvasTkAgg(f, self)
        self.canvas.show()
        
        #self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        toolbar.update()
        toolbar.grid(row=0, column=1, sticky="ew")
        self.canvas.get_tk_widget()
        self.canvas._tkcanvas.grid(row=1, column=1, rowspan=20, columnspan=20, sticky="nsew")        
        self.pool.draw_RA2I(self.ax)
        self.canvas.draw()

class PageFive(tk.Frame):
    def __init__(self, parent, controller, pool):
        tk.Frame.__init__(self, parent)
        
        self.pool = pool        
        
        button1 = ttk.Button(self, text="Graph",
                             command=lambda: controller.show_page(PageOne))        
        
        button2 = ttk.Button(self, text="RA2I", 
                             command=lambda: controller.show_page(PageFour))
        
        button3 = ttk.Button(self, text="Distribution", 
                             command=lambda: controller.show_page(PageTwo))
                             
        button4 = ttk.Button(self, text="Synapses", 
                             command=lambda: controller.show_page(PageThree))
        
        button5 = ttk.Button(self, text="Raster", 
                             command=lambda: controller.show_page(PageSix))
        
        button6 = ttk.Button(self, text="Ordered spikes", 
                             command=lambda: controller.show_page(PageSeven))
                             
        button1.grid(row=0, column=0)
        button2.grid(row=1, column=0)
        button3.grid(row=2, column=0)
        button4.grid(row=3, column=0)
        button5.grid(row=4, column=0)
        button6.grid(row=5, column=0)
        
        f = figure(figsize=(6,3), dpi=200)
        self.ax = f.add_subplot(111)
        
        
        self.canvas = FigureCanvasTkAgg(f, self)
        self.canvas.show()
        
        #self.canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2TkAgg(self.canvas, self)
        toolbar.update()
        toolbar.grid(row=0, column=1, sticky="ew")
        self.canvas.get_tk_widget()
        self.canvas._tkcanvas.grid(row=1, column=1, rowspan=20, columnspan=20, sticky="nsew")        
        self.pool.draw_I2RA(self.ax)
        self.canvas.draw()

class PageSix(tk.Frame):
    def __init__(self, parent, controller, pool):
        tk.Frame.__init__(self, parent)
        
        button1 = ttk.Button(self, text="Graph", 
                            command=lambda: controller.show_page(PageOne))        
       
        button2 = ttk.Button(self, text="RA2I", 
                             command=lambda: controller.show_page(PageFour))
        
        button3 = ttk.Button(self, text="I2RA", 
                             command=lambda: controller.show_page(PageFive))
         
        button4 = ttk.Button(self, text="Distribution", 
                             command=lambda: controller.show_page(PageTwo))
        
        button5 = ttk.Button(self, text="Synapses", 
                             command=lambda: controller.show_page(PageThree))
        
        button6 = ttk.Button(self, text="Ordered spikes", 
                             command=lambda: controller.show_page(PageSeven))
    
        button1.grid(row=0, column=0)                                
        button2.grid(row=1, column=0)
        button3.grid(row=2, column=0)
        button4.grid(row=3, column=0)
        button5.grid(row=4, column=0)
        button6.grid(row=5, column=0)
        
        self.pool = pool
        
        f_dend = Figure(figsize=(6,6), dpi=100)
        f_soma = Figure(figsize=(6,6), dpi=100)
        
        self.ax_dend = f_dend.add_subplot(111)
        self.ax_soma = f_soma.add_subplot(111)
        
        self.canvas_dend = FigureCanvasTkAgg(f_dend, self)
        self.canvas_dend.show()
        self.canvas_soma = FigureCanvasTkAgg(f_soma, self)
        self.canvas_soma.show()
        
        toolbar_dend = NavigationToolbar2TkAgg(self.canvas_dend, self)
        toolbar_dend.update()
        toolbar_dend.grid(row=0, column=1, sticky="ew")

        toolbar_soma = NavigationToolbar2TkAgg(self.canvas_soma, self)
        toolbar_soma.update()
        toolbar_soma.grid(row=0, column=21, sticky="ew")

        
        #toolbar.pack(side=tk.TOP)
        #self.canvas_dist.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=False)
        self.canvas_dend.get_tk_widget()
        self.canvas_soma.get_tk_widget()
        #self.canvas_dist._tkcanvas.pack()
        self.canvas_dend._tkcanvas.grid(row=1, rowspan = 20, column=1, columnspan = 20, sticky="nsew")        
        self.canvas_soma._tkcanvas.grid(row=1, rowspan = 20, column=21, sticky="nsew")        
        
                  
         
        #self.canvas_dist._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=False)        
        
        self.update(controller)
        
    def update(self, controller):
        
        if controller.poolIsUpdated():
            #print "Distribution page update"            
            self.pool.draw_raster_dend(self.ax_dend)
            self.canvas_dend.draw()
            self.pool.draw_raster_soma(self.ax_soma)
            self.canvas_soma.draw()            
        controller.after(5000, self.update, controller)
            
class PageSeven(tk.Frame):
    def __init__(self, parent, controller, pool):
        tk.Frame.__init__(self, parent)
        
        button1 = ttk.Button(self, text="Graph", 
                            command=lambda: controller.show_page(PageOne))        
       
        button2 = ttk.Button(self, text="RA2I", 
                             command=lambda: controller.show_page(PageFour))
        
        button3 = ttk.Button(self, text="I2RA", 
                             command=lambda: controller.show_page(PageFive))
         
        button4 = ttk.Button(self, text="Distribution", 
                             command=lambda: controller.show_page(PageTwo))
        
        button5 = ttk.Button(self, text="Synapses", 
                             command=lambda: controller.show_page(PageThree))
        
        button6 = ttk.Button(self, text="Raster", 
                             command=lambda: controller.show_page(PageSix))
    
        button1.grid(row=0, column=0)                                
        button2.grid(row=1, column=0)
        button3.grid(row=2, column=0)
        button4.grid(row=3, column=0)
        button5.grid(row=4, column=0)
        button6.grid(row=5, column=0)
        
        self.pool = pool
        
        f_dend = Figure(figsize=(6,6), dpi=100)
        
        self.ax_dend = f_dend.add_subplot(111)
        
        self.canvas_dend = FigureCanvasTkAgg(f_dend, self)
        self.canvas_dend.show()
        
        f_soma = Figure(figsize=(6,6), dpi=100)
        
        self.ax_soma = f_soma.add_subplot(111)
        
        self.canvas_soma = FigureCanvasTkAgg(f_soma, self)
        self.canvas_soma.show()
        
        toolbar_dend = NavigationToolbar2TkAgg(self.canvas_dend, self)
        toolbar_dend.update()
        toolbar_dend.grid(row=0, column=1, sticky="ew")

        toolbar_soma = NavigationToolbar2TkAgg(self.canvas_soma, self)
        toolbar_soma.update()
        toolbar_soma.grid(row=0, column=21, sticky="ew")

        
        #toolbar.pack(side=tk.TOP)
        #self.canvas_dist.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=False)
        self.canvas_dend.get_tk_widget()
        self.canvas_soma.get_tk_widget()
        #self.canvas_dist._tkcanvas.pack()
        self.canvas_dend._tkcanvas.grid(row=1, rowspan = 20, column=1, columnspan = 20, sticky="nsew")        
        self.canvas_soma._tkcanvas.grid(row=1, rowspan = 20, column=21, columnspan = 20, sticky="nsew")        
        
                  
         
        #self.canvas_dist._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=False)        
        
        self.update(controller)
        
    def update(self, controller):
        
        if controller.poolIsUpdated():
            #print "Distribution page update"            
            self.pool.draw_ordered_dend(self.ax_dend)
            self.canvas_dend.draw()
            self.pool.draw_ordered_soma(self.ax_soma)
            self.canvas_soma.draw()
                       
        controller.after(5000, self.update, controller)


class Pool:
    def __init__(self):
        self.file_weights = filename_weights
        self.file_time_dend = filename_time_dend
        self.file_time_soma = filename_time_soma
                
        self.std = []
        self.mean = []
        self.time = []
        self.active_synapses = []
        self.supersynapses = []

        self.mod_time_RA_RA = 0
        self.mod_time_wgh = 0
        self.mod_time_time = 0
                
        
        
        self.pos = {}
        self.node_color = []
        self.edgewidth = []
        
                 
        (self.N_RA, self.RA2I_edges, self.RA2I_weights) = read.get_RA2I_graph(filenameRA_I)
        (self.N_I, self.I2RA_edges, self.I2RA_weights) = read.get_I2RA_graph(self.N_RA, filenameI)
        (unused, self.edges, self.weights) = read.get_RA2RA_graph(filenameRA_RA)
        
        #print self.I2RA_edges[0][0]
        #print self.edges        
        #print self.RA2I_edges
        # read coordinates
        (xx_RA, yy_RA) = read.read_coordinates(filename_xy_RA)
        (xx_I, yy_I) = read.read_coordinates(filename_xy_I)
    
        # fill dictionary for positions of neurons
        
        for i in range(self.N_RA):
            self.pos[i] = (xx_RA[i], yy_RA[i])
        for i in range(self.N_I):   
            self.pos[i + self.N_RA] = (xx_I[i], yy_I[i])
        
        # assign colors to nodes
        for i in range(self.N_RA):
            self.node_color.append("g")
        for i in range(self.N_I):
            self.node_color.append("r")
        
        #self.weights = map(lambda x: x*10, self.weights)        
        
        # create bidirectional graph
        
        self.G = nx.DiGraph()
        self.G.add_nodes_from(self.pos.keys())
       
        # biderectional graph for fixed synapses       
        self.G_I2RA = nx.DiGraph()
        self.G_I2RA.add_nodes_from(self.pos.keys())
        
        self.G_RA2I = nx.DiGraph()
        self.G_RA2I.add_nodes_from(self.pos.keys())
                
        
        self.edge_color_I2RA = []        
        self.edge_color_RA2I = []        
        
        for i in range(len(self.I2RA_edges)):
            self.G_I2RA.add_edge(self.I2RA_edges[i][0], self.I2RA_edges[i][1], weight=self.I2RA_weights[i])
            self.edge_color_I2RA.append('r')
            
        for i in range(len(self.RA2I_edges)):
            self.G_RA2I.add_edge(self.RA2I_edges[i][0], self.RA2I_edges[i][1], weight=self.RA2I_weights[i])
            self.edge_color_RA2I.append('g')
            
        self.edgewidth_I2RA = [ d['weight'] for (u,v,d) in self.G_I2RA.edges(data=True)]
        self.edgewidth_RA2I = [ d['weight'] for (u,v,d) in self.G_RA2I.edges(data=True)]

    def needs_update(self):
        return (self.mod_time_RA_RA!=os.path.getmtime(filenameRA_RA) and self.mod_time_wgh!=os.path.getmtime(self.file_weights)
            and self.mod_time_time!=os.path.getmtime(self.file_time_soma))
            
    def update(self):
        (self.N_RA, w) = read.read_weights(self.file_weights)
        self.mod_time_RA_RA = os.path.getmtime(filenameRA_RA)
        self.mod_time_wgh = os.path.getmtime(self.file_weights)
        self.mod_time_time = os.path.getmtime(self.file_time_soma)
        
        self.weights = [item for sublist in w for item in sublist]
        
        for i in range(len(self.weights)):
            if self.weights[i] < sys.float_info.epsilon:
                self.weights[i] = 0
        
        self.hist, self.bin_edges = np.histogram(self.weights, bins=100)
        
        # update mean and std dependences on time        
        
        sigma = std(self.weights)
        mu = mean(self.weights)
        
        self.mean.append(mu)
        self.std.append(sigma)        

        # update time pattern of spikes
         
        self.trial_number, simulation_time, spike_times_dend, neuronID_dend = read.read_time_info(self.file_time_dend)
        
        unused, unused, spike_times_soma, neuronID_soma = read.read_time_info(self.file_time_soma)
        (unused, RA_super_targets, unused) = read.read_connections(filenameRA_RA_super)

        # extract strongly connected neurons
        superTargets = set([element for sublist in RA_super_targets for element in sublist])
        superSources = set([i for i in xrange(len(RA_super_targets)) if len(RA_super_targets[i]) > 0])
        stronglyConnected = superTargets | superSources


        self.time.append(simulation_time / 1000)
        
        self.spike_times_soma = [item for sublist in spike_times_soma for item in sublist]
        self.ID_soma = [item for sublist in neuronID_soma for item in sublist]
        
        self.spike_times_dend = [item for sublist in spike_times_dend for item in sublist]
        self.ID_dend = [item for sublist in neuronID_dend for item in sublist]
        
        # check that spike times are in range     
        if len(self.spike_times_soma) > 1:
            if min(self.spike_times_soma) < 0:
                print "Somatic spike time is negative! Min somatic spike = ", min(self.spike_times_soma)
            if max(self.spike_times_soma) > TRIAL_DURATION:
                print "Somatic spike time exceeds trial duration!"
                
        if len(self.spike_times_dend) > 1:
            if min(self.spike_times_dend) < 0:
                print "Dendritic spike time is negative! Min dendritic spike = ", min(self.spike_times_dend)
            if max(self.spike_times_dend) > TRIAL_DURATION:
                print "Dendritic spike time exceeds trial duration!"
        
        # extract only spikes of strongly connected neurons
        self.ID_soma_strong = [i for i in self.ID_soma if i in stronglyConnected]
        self.spike_times_soma_strong = [t for ind, t in enumerate(self.spike_times_soma) if self.ID_soma[ind] in stronglyConnected]
        
        self.ID_dend_strong = [i for i in self.ID_dend if i in stronglyConnected]
        self.spike_times_dend_strong = [t for ind, t in enumerate(self.spike_times_dend) if self.ID_dend[ind] in stronglyConnected]
        
        self.order_dend()
        self.order_soma()
        
        #ID = range(self.N_RA)
        #self.spike_times = sorted(spike_times)
        #temp = sorted(zip(ID, spike_times), key=lambda par:par[1])
        
        #self.ID, self.spike_times = zip(*temp)
                
        (unused, self.edges, self.active_weights) = read.get_RA2RA_graph(filenameRA_RA)
        
        active_synapses = len(self.active_weights)
        #print active_synapses
        #print self.active_weights
        supersynapses = sum(1 for w in self.active_weights if w > SUPERSYNAPSE_THRESHOLD)
        #print supersynapses
        self.active_synapses.append(active_synapses)
        self.supersynapses.append(supersynapses)
        #self.weights = map(lambda x: x*10, self.weights)                
        
        # create new empty graph
        self.G_temp = nx.DiGraph()
        self.G_temp.add_nodes_from(self.G.nodes(data=True)) # read nodes from the previous graph
        self.G = self.G_temp        
        for i in range(len(self.edges)):
            self.G.add_edge(self.edges[i][0], self.edges[i][1], 
                            weight=self.active_weights[i])
        
        
        self.edgewidth = [ d['weight'] for (u,v,d) in self.G.edges(data=True)]

        self.update_edge_color()
        
    def update_edge_color(self):
        # assign color to edges
        self.edge_color = []        
        
        for i in range(len(self.edges)):
            if (self.edgewidth[i] >= 1 * SUPERSYNAPSE_THRESHOLD):
                self.edge_color.append("b")
            else:
                self.edge_color.append("g")
                
    def draw_raster_dend(self, axes):
        #print "Refresh spike times!"
        axes.clear()
        
        for i in range(len(self.spike_times_dend)):
            axes.vlines(self.spike_times_dend[i], self.ID_dend[i]-0.5, self.ID_dend[i]+0.5)
        
        #axes.scatter(self.spike_times, self.ID)
        axes.set_ylabel("Neuron ID")
        axes.set_xlabel("spike time dendrite (ms)")
       
        if len(self.ID_dend)>1:
            axes.set_ylim(0, max(self.ID_dend))
            axes.set_xlim(0, TRIAL_DURATION)

    def draw_raster_soma(self, axes):
        #print "Refresh spike times!"
        axes.clear()
        
        for i in range(len(self.spike_times_soma)):
            axes.vlines(self.spike_times_soma[i], self.ID_soma[i]-0.5, self.ID_soma[i]+0.5)
        
        #axes.scatter(self.spike_times, self.ID)
        axes.set_ylabel("Neuron ID")
        axes.set_xlabel("spike time soma (ms)")
       
        #print "somatic spikes: ", self.spike_times_soma

        if len(self.ID_soma)>1:
            axes.set_ylim(0, max(self.ID_soma))
            axes.set_xlim(0, TRIAL_DURATION)
    
    def order_dend(self):
        if len(self.spike_times_dend_strong) > 1:
            first_spike = min(self.spike_times_dend_strong)
    
            spike_times = map(lambda x: x - first_spike, self.spike_times_dend_strong)
            neuron_fired = list(self.ID_dend_strong)
            
            
            self.spike_times_ordered_dend, self.neuron_fired_ordered_dend = zip(*sorted(zip(spike_times, neuron_fired)))
           
            self.random_ID_dend = self.remap_ID(self.neuron_fired_ordered_dend)
        
    def order_soma(self):
        if len(self.spike_times_soma_strong) > 1 :
            first_spike = min(self.spike_times_soma_strong)
    
            spike_times = map(lambda x: x - first_spike, self.spike_times_soma_strong)
            neuron_fired = list(self.ID_soma_strong)
            
            
            self.spike_times_ordered_soma, self.neuron_fired_ordered_soma = zip(*sorted(zip(spike_times, neuron_fired)))
           
            self.random_ID_soma = self.remap_ID(self.neuron_fired_ordered_soma)
            
    
    
    def remap_ID(self, order):
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
     
            
    def draw_ordered_dend(self, axes):
        #print "Refresh spike times!"
        if len(self.spike_times_dend_strong) > 1:
            axes.clear()
            
            for i in range(len(self.spike_times_ordered_dend)):
                axes.vlines(self.spike_times_ordered_dend[i], self.random_ID_dend[i]-0.5, self.random_ID_dend[i]+0.5)
        
        #ax.scatter(spike_times, random_ID)
            axes.set_yticks(self.random_ID_dend)
            
            axes.set_yticklabels(self.neuron_fired_ordered_dend)
            axes.set_xlim([-5, max(self.spike_times_ordered_dend)+5])
            axes.set_ylabel("real neuron ID")
            axes.set_xlabel("relative spike time (ms)")
            axes.set_title("Ordered dendritic spikes")
     
    def draw_ordered_soma(self, axes):
        #print "Refresh spike times!"
        if len(self.spike_times_soma_strong) > 1:
            axes.clear()
            
            for i in range(len(self.spike_times_ordered_soma)):
                axes.vlines(self.spike_times_ordered_soma[i], self.random_ID_soma[i]-0.5, self.random_ID_soma[i]+0.5)
        
        #ax.scatter(spike_times, random_ID)
            axes.set_yticks(self.random_ID_soma)
            
            axes.set_yticklabels(self.neuron_fired_ordered_soma)
            axes.set_xlim([-5, max(self.spike_times_ordered_soma)+5])
            axes.set_ylabel("real neuron ID")
            axes.set_xlabel("relative spike time (ms)")
            axes.set_title("Ordered somatic spikes")
       
    
    def draw_active_synapses(self, axes):
        #print "Resfresh synapses!"
        #print self.active_synapses
        #print self.supersynapses        
        axes.clear()
        axes.plot(self.time, self.active_synapses, 'g', label="Active")
        
        
        #axes.legend()
        
        #axes.legend([line_active, line_super], ["Active", "Supersynapses"])
        axes.set_xlabel("time (s)")
        axes.set_ylabel("Number of active synapses")
    
    def draw_supersynapses(self, axes):
        axes.clear()
        axes.plot(self.time, self.supersynapses, 'b', label="Supersynapses")
        
        axes.set_xlabel("time (s)")
        axes.set_ylabel("Number of super synapses")
        
    
    def draw_hist(self, axes):
        #print "Refresh weights!"
        #print type(axes)
        width = 0.7 * (self.bin_edges[1] - self.bin_edges[0])
        center = (self.bin_edges[:-1] + self.bin_edges[1:]) / 2
        axes.clear()
        axes.bar(center, self.hist, align="center", width=width)
        axes.set_ylabel("Number of synapses")
        axes.set_xlabel("Synaptic weight")
        axes.ticklabel_format(style="sci", axis="x", scilimits=(0,0))
        #axes.figure.canvas.draw()
        #controller.after(5000, self.draw, axes, controller)            
        
    def draw_time_statistics(self, ax1, ax2):
        #print "Refresh time statistics"
        
        ax1.clear()
        ax1.plot(self.time, self.mean)
        ax1.set_ylabel("mean")
        ax1.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
        ax1.get_xaxis().set_visible(False)        
        
        ax2.clear()
        ax2.plot(self.time, self.std)
        ax2.set_ylabel("std")
        ax2.set_xlabel("time (s)")
        ax2.ticklabel_format(style="sci", axis="y", scilimits=(0,0))
        ax2.ticklabel_format(style="sci", axis="x", scilimits=(0,0))
        
    def draw_graph(self, axes):
        axes.clear()
        #print "Refresh graph!"
        nx.draw_networkx(self.G, ax=axes, pos=self.pos, node_size = 130, width = self.edgewidth, 
                    edge_color = self.edge_color, node_color=self.node_color, with_labels=True, font_size=FONT_SIZE)
        axes.set_xlim((-5,SIDE+5))
        axes.set_ylim((-5,SIDE+5))

    def draw_RA2I(self, axes):
        axes.clear()
        #print "Refresh graph!"
        nx.draw_networkx(self.G_RA2I, ax=axes, pos=self.pos, node_size = 130, width = self.edgewidth_RA2I, 
                    edge_color = self.edge_color_RA2I, node_color=self.node_color, with_labels=True, font_size=FONT_SIZE)
        axes.set_xlim((-5,SIDE+5))
        axes.set_ylim((-5,SIDE+5))
        
    def draw_I2RA(self, axes):
        axes.clear()
        #print "Refresh graph!"
        nx.draw_networkx(self.G_I2RA, ax=axes, pos=self.pos, node_size = 130, width = self.edgewidth_I2RA, 
                    edge_color = self.edge_color_I2RA, node_color=self.node_color, with_labels=True, font_size=FONT_SIZE)
        axes.set_xlim((-5,SIDE+5))
        axes.set_ylim((-5,SIDE+5))

def mean(array):
    n = len(array)
    if n < 1:
        raise ValueError("Mean requires at least one element in array")
    
    return sum(array)/n
 
def std(array):
    #copy = array.flatten()
    n = len(array)
    if n < 1:
        raise ValueError("Std requires at least one element in array")
    mu = mean(array)
    return sqrt(sum([(x - mu) ** 2 for x in array])/n)

app = My_gui()
app.mainloop()
c
