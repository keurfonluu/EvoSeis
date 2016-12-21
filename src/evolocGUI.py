# -*- coding: utf-8 -*-
#==============================================================================
# evoloc_gui.py
#------------------------------------------------------------------------------
# evoloc_gui displays a graphical user interface for EvoLoc.
#
# Created by
#       Keurfon Luu <keurfon.luu@mines-paristech.fr>
#       MINES ParisTech - Centre de Géosciences
#       PSL - Research University
#
# Last updated
#       2016-12-22 00:32
#==============================================================================

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from mpl_toolkits.mplot3d import axes3d
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from datetime import datetime
from HoverInfo import HoverInfo

import multiprocessing, os, sys
if sys.version_info[0] < 3:
    import Tkinter as tk
    import tkFileDialog as tkfile
    import tkMessageBox as tkmessage
    from ttk import *
    import tkFont as font
    FileNotFoundError = IOError
else:
    import tkinter as tk
    import tkinter.filedialog as tkfile
    import tkinter.messagebox as tkmessage
    from tkinter.ttk import *
    from tkinter import font

# Evoloc GUI
#============
class evolocGUI():

    master = None

    # Initialize master window
    #==========================
    def __init__(self, master):
        self.master = master
        master.title("EvoLoc")
        master.protocol("WM_DELETE_WINDOW", self.close_window)
        master.geometry("1000x620")
        master.minsize(1000, 620)

        default_font = font.nametofont("TkDefaultFont")
        default_font.configure(family = "Helvetica", size = 9)
        master.option_add("*Font", default_font)

        self.menubar()
        self.define_variables()
        self.trace_variables()
        self.init_variables()
        self.init_figure()
        self.frame1()
        self.frame2()
        self.frame3()
        self.frame4()
        self.frame5()
        self.frame6()
        self.frame7()
        self.footer()

    # Menu bar
    #==========
    def menubar(self):
        menubar = tk.Menu(self.master)

        # File
        filemenu = tk.Menu(menubar, tearoff = 0)
        filemenu.add_command(label = "Import parameters", command = self.import_parameters)
        filemenu.add_command(label = "Export parameters", command = self.export_parameters)
        filemenu.add_separator()
        filemenu.add_command(label = "Exit", command = self.close_window)

        # Help
        helpmenu = tk.Menu(menubar, tearoff = 0)
        helpmenu.add_command(label = "About", command = self.about_window)

        # Display menu bar
        menubar.add_cascade(label = "File", menu = filemenu)
        menubar.add_cascade(label = "Help", menu = helpmenu)
        self.master.config(menu = menubar)

    # Frame 1: Import velocity models
    #=================================
    def frame1(self):
        self.frame1 = tk.Frame(self.master, borderwidth = 2, relief = "groove")
        self.frame1.place(bordermode = "outside", relwidth = 0.49, relheight = 0.3, relx = 0, x = 5, y = 5, anchor = "nw")

        # P-velocity label
        velp_label = tk.Label(self.frame1, text = "P-velocity model (m/s)")

        # P-velocity entry
        velp_entry = tk.Entry(self.frame1, textvariable = self.velp_filename)

        # P-velocity import button
        velp_import_button = tk.Button(self.frame1, text = "Import", command = self.import_velp_file)

        # S-velocity label
        vels_label = tk.Label(self.frame1, text = "S-velocity model (m/s)")

        # S-velocity entry
        vels_entry = tk.Entry(self.frame1, textvariable = self.vels_filename)

        # S-velocity import button
        vels_import_button = tk.Button(self.frame1, text = "Import", command = self.import_vels_file)

        # Grid size (nz, nx, ny)
        nz_label = tk.Label(self.frame1, text = "nz")
        nx_label = tk.Label(self.frame1, text = "nx")
        ny_label = tk.Label(self.frame1, text = "ny")

        nz_entry = tk.Entry(self.frame1, textvariable = self.nz, justify = "right")
        nx_entry = tk.Entry(self.frame1, textvariable = self.nx, justify = "right")
        ny_entry = tk.Entry(self.frame1, textvariable = self.ny, justify = "right")

        # Mesh size (dz, dx, dy)
        dz_label = tk.Label(self.frame1, text = "dz (m)")
        dx_label = tk.Label(self.frame1, text = "dx (m)")
        dy_label = tk.Label(self.frame1, text = "dy (m)")

        dz_entry = tk.Entry(self.frame1, textvariable = self.dz, justify = "right")
        dx_entry = tk.Entry(self.frame1, textvariable = self.dx, justify = "right")
        dy_entry = tk.Entry(self.frame1, textvariable = self.dy, justify = "right")

        # View velocity models buttons
        view_vp_button = tk.Button(self.frame1, text = "View Vp >>", command = self.view_vp)
        view_vs_button = tk.Button(self.frame1, text = "View Vs >>", command = self.view_vs)

        # Layout
        velp_label.place(relx = 0, x = 5, y = 5, anchor = "nw")
        velp_entry.place(relwidth = 0.78, relx = 0, x = 5, y = 22, anchor = "nw")
        velp_import_button.place(relwidth = 0.19, relx = 0.8, y = 18, anchor = "nw")

        vels_label.place(relx = 0, x = 5, y = 50, anchor = "nw")
        vels_entry.place(relwidth = 0.78, relx = 0, x = 5, y = 67, anchor = "nw")
        vels_import_button.place(relwidth = 0.19, relx = 0.8, y = 63, anchor = "nw")

        nz_label.place(relx = 0, x = 5, y = 95, anchor = "nw")
        nx_label.place(relx = 0.15, x = 5, y = 95, anchor = "nw")
        ny_label.place(relx = 0.3, x = 5, y = 95, anchor = "nw")

        nz_entry.place(relwidth = 0.15, relx = 0, x = 5, y = 112, anchor = "nw")
        nx_entry.place(relwidth = 0.15, relx = 0.15, x = 5, y = 112, anchor = "nw")
        ny_entry.place(relwidth = 0.15, relx = 0.3, x = 5, y = 112, anchor = "nw")

        dz_label.place(relx = 0.55, x = -5, y = 95, anchor = "nw")
        dx_label.place(relx = 0.7, x = -5, y = 95, anchor = "nw")
        dy_label.place(relx = 0.85, x = -5, y = 95, anchor = "nw")

        dz_entry.place(relwidth = 0.15, relx = 0.7, x = -5, y = 112, anchor = "ne")
        dx_entry.place(relwidth = 0.15, relx = 0.85, x = -5, y = 112, anchor = "ne")
        dy_entry.place(relwidth = 0.15, relx = 1, x = -5, y = 112, anchor = "ne")

        view_vp_button.place(relwidth = 0.3, relx = 0.1, y = 143, anchor = "nw")
        view_vs_button.place(relwidth = 0.3, relx = 0.9, y = 143, anchor = "ne")

    # Frame 2: View velocity models
    #===============================
    def frame2(self):
        self.frame2 = tk.Frame(self.master, borderwidth = 2, relief = "groove", bg = "white")
        self.frame2.place(bordermode = "outside", relwidth = 0.49, relheight = 0.52, relx = 1, x = -5, y = 5, anchor = "ne")
        self.frame2.first_view = True



    # Frame 3: Import data files
    #============================
    def frame3(self):
        self.frame3 = tk.Frame(self.master, borderwidth = 2, relief = "groove")
        self.frame3.place(bordermode = "outside", relwidth = 0.49, relheight = 0.2, relx = 0, rely = 0.32, x = 5, y = 5, anchor = "nw")

        # Stations file label
        stat_label = tk.Label(self.frame3, text = "Stations file (depth in m)")

        # Stations file entry
        stat_entry = tk.Entry(self.frame3, textvariable = self.stat_filename)

        # Stations file import button
        stat_import_button = tk.Button(self.frame3, text = "Import", command = self.import_stat_file)

        # Arrival times label
        tobs_label = tk.Label(self.frame3, text = "Arrival times file (s)")

        # Arrival times entry
        tobs_entry = tk.Entry(self.frame3, textvariable = self.tobs_filename)

        # Arrival times import button
        tobs_import_button = tk.Button(self.frame3, text = "Import", command = self.import_tobs_file)

        # Reference date label
        refdate_label = tk.Label(self.frame3, text = "Reference date (YYYY-MM-DD)")

        # Reference date entry
        refdate_entry = tk.Entry(self.frame3, textvariable = self.refdate)

        # Station Z coordinates label
        zcoord_label = tk.Label(self.frame3, text = "Z coordinates")

        # Station Z coordinates options menu
        zcoord_option = tk.OptionMenu(self.frame3, self.zcoord, "Depth", "Elevation")

        # Layout
        stat_label.place(relx = 0, x = 5, y = 5, anchor = "nw")
        stat_entry.place(relwidth = 0.78, relx = 0, x = 5, y = 22, anchor = "nw")
        stat_import_button.place(relwidth = 0.19, relx = 0.8, y = 18, anchor = "nw")

        tobs_label.place(relx = 0, x = 5, y = 50, anchor = "nw")
        tobs_entry.place(relwidth = 0.78, relx = 0, x = 5, y = 67, anchor = "nw")
        tobs_import_button.place(relwidth = 0.19, relx = 0.8, y = 63, anchor = "nw")

        zcoord_label.place(relx = 0, x = 5, y = 95, anchor = "nw")
        zcoord_option.place(relwidth = 0.2, height = 25, x = 100, y = 92, anchor = "nw")

        refdate_label.place(relx = 0.45, y = 95, anchor = "nw")
        refdate_entry.place(relwidth = 0.15, relx = 0.45, x = 185, y = 95, anchor = "nw")

        # Hover
        HoverInfo(refdate_label, "If no reference date have been used, leave it empty.")

    # Frame 4: Choose formulation
    #=============================
    def frame4(self):
        self.frame4 = tk.Frame(self.master, borderwidth = 2, relief = "groove")
        self.frame4.place(bordermode = "outside", relwidth = 0.49, relheight = 0.17, relx = 0, rely = 0.54, x = 5, y = 5, anchor = "nw")
        self.frame4.first_view = True

        # Probabilistic formulation label
        formulation_label = tk.Label(self.frame4, text = "Probabilistic formulation")

        # Probabilistic formulation selection
        self.bayes = tk.Radiobutton(self.frame4, text = "Bayesian", variable = self.do_bayesh, value = False, command = self.bayes_widget)
        self.bayesh = tk.Radiobutton(self.frame4, text = "Hierarchical bayesian", variable = self.do_bayesh, value = True, command = self.bayesh_widget)

        # Layout
        formulation_label.place(relx = 0, x = 5, y = 5, anchor = "nw")
        self.bayes.place(relx = 0, x = 5, y = 22, anchor = "nw")
        self.bayesh.place(relx = 0.5, x = 5, y = 22, anchor = "nw")

        if not self.do_bayesh.get():
            self.bayes.invoke()
        else:
            self.bayesh.invoke()

    # Frame 5: Configure coordinates
    #================================
    def frame5(self):
        self.frame5 = tk.Frame(self.master, borderwidth = 2, relief = "groove")
        self.frame5.place(bordermode = "outside", relwidth = 0.49, relheight = 0.2, relx = 0, rely = 0.73, x = 5, y = 5, anchor = "nw")
        self.frame5.first_view = True

        # Coordinates unit label
        coord_label = tk.Label(self.frame5, text = "Coordinates unit")

        # Coordinates unit selection
        self.meter = tk.Radiobutton(self.frame5, text = "Meter", variable = self.unit, value = 0, command = self.meter_widget)
        self.degree = tk.Radiobutton(self.frame5, text = "Degree", variable = self.unit, value = 1, command = self.degree_widget)

        # Layout
        coord_label.place(relx = 0, x = 5, y = 5, anchor = "nw")
        self.meter.place(relx = 0, x = 5, y = 22, anchor = "nw")
        self.degree.place(relx = 0.5, x = 5, y = 22, anchor = "nw")

        if self.unit.get() == 0:
            self.meter.invoke()
        else:
            self.degree.invoke()

    # Frame 6: Inversion parameters
    #===============================
    def frame6(self):
        self.frame6 = tk.Frame(self.master, borderwidth = 2, relief = "groove")
        self.frame6.place(bordermode = "outside", relwidth = 0.49, relheight = 0.21, relx = 1, rely = 0.54, x = -5, y = 5, anchor = "ne")
        self.frame6.first_view = True

        # Optimization method label
        method_label = tk.Label(self.frame6, text = "Optimization method")

        # Optimization method selection
        self.cpso = tk.Radiobutton(self.frame6, text = "Competitive PSO", variable = self.method, value = 1, command = self.cpso_widget)
        self.de = tk.Radiobutton(self.frame6, text = "Differential Evolution", variable = self.method, value = 2, command = self.de_widget)

        # Population size label
        population_label = tk.Label(self.frame6, text = "Population size")

        # Population size spinbox
        population_spinbox = tk.Spinbox(self.frame6, from_ = 2, to_ = 9999, increment = 1, textvariable = self.popsize)

        # Iterations label
        itermax_label = tk.Label(self.frame6, text = "Maximum number of iterations")

        # Iterations spinbox
        itermax_spinbox = tk.Spinbox(self.frame6, from_ = 1, to_ = 9999, increment = 1, textvariable = self.itermax)

        # Layout
        method_label.place(relx = 0, x = 5, y = 5, anchor = "nw")
        self.cpso.place(relx = 0, x = 5, y = 22, anchor = "nw")
        self.de.place(relx = 0.5, x = 5, y = 22, anchor = "nw")

        if self.method.get() == 1:
            self.cpso.invoke()
        else:
            self.de.invoke()

        population_label.place(relx = 0, x = 5, y = 100, anchor = "nw")
        population_spinbox.place(relwidth = 0.1, relx = 0.21, y = 100, anchor = "nw")

        itermax_label.place(relx = 0.5, x = 3, y = 100, anchor = "nw")
        itermax_spinbox.place(relwidth = 0.1, relx = 0.87, y = 100, anchor = "nw")

    # Frame 7: Output
    #=================
    def frame7(self):
        self.frame7 = tk.Frame(self.master, borderwidth = 2, relief = "groove")
        self.frame7.place(bordermode = "outside", relwidth = 0.49, relheight = 0.16, relx = 1, rely = 0.77, x = -5, y = 5, anchor = "ne")

        # Output directory label
        output_label = tk.Label(self.frame7, text = "Output directory")

        # Output directory entry
        output_entry = tk.Entry(self.frame7, textvariable = self.outdir)

        # Output directory open button
        output_open_button = tk.Button(self.frame7, text = "Open", command = self.open_outdir)

        # PDF check button
        pdf_button = tk.Checkbutton(self.frame7, text = "Compute Probability Density Function", variable = self.pdf)

        # Layout
        output_label.place(relx = 0, x = 5, y = 5, anchor = "nw")
        output_entry.place(relwidth = 0.78, relx = 0, x = 5, y = 22, anchor = "nw")
        output_open_button.place(relwidth = 0.19, relx = 0.8, y = 18, anchor = "nw")

        pdf_button.place(relx = 0, x = 5, y = 60, anchor = "nw")


    # Footer: Locate and exit
    #=========================
    def footer(self):

        # Threads label
        num_threads_label = tk.Label(self.master, text = "Number of threads")

        # Threads spinbox
        num_threads_spinbox = tk.Spinbox(self.master, from_ = 1, to_ = multiprocessing.cpu_count(), increment = 1, textvariable = self.num_threads)

        # Locate button
        locate_button = tk.Button(self.master, text = "Locate", command = self.locate)

        # Exit button
        exit_button = tk.Button(self.master, text = "Exit", command = self.close_window)

        # Layout
        num_threads_label.place(relx = 0, rely = 1, x = 5, y = -8, anchor = "sw")
        num_threads_spinbox.place(width = 30, rely = 1, x = 120, y = -8, anchor = "sw")
        locate_button.place(relwidth = 0.1, relx = 0.9, rely = 1, x = -5, y = -5, anchor = "se")
        exit_button.place(relwidth = 0.1, relx = 1, rely = 1, x = -5, y = -5, anchor = "se")

    # Widget: Bayesian formulation
    #==============================
    def bayes_widget(self):
        if not self.frame4.first_view:
            self.frame4_in.forget()
            self.frame4_in = tk.Frame(self.frame4, borderwidth = 0)
        else:
            self.frame4_in = tk.Frame(self.frame4, borderwidth = 0)
            self.frame4.first_view = False
        self.frame4_in.place(relwidth = 0.99, height = 40, relx = 0, y = 50, x = 5, anchor = "nw")

        # Uncertainties file label
        sigma_label = tk.Label(self.frame4_in, text = "Uncertainties file (s)")

        # Uncertainties file entry
        sigma_entry = tk.Entry(self.frame4_in, textvariable = self.sigma_filename)

        # Uncertainties file import button
        sigma_import_button = tk.Button(self.frame4_in, text = "Import", command = self.import_sigma_file)

        # Layout
        sigma_label.place(relx = 0, x = 0, y = 0, anchor = "nw")
        sigma_entry.place(relwidth = 0.79, relx = 0, x = 0, y = 17, anchor = "nw")
        sigma_import_button.place(relwidth = 0.1925, relx = 0.799, y = 13, anchor = "nw")

    # Widget: Hierarchical bayesian formulation
    #===========================================
    def bayesh_widget(self):
        if not self.frame4.first_view:
            self.frame4_in.forget()
            self.frame4_in = tk.Frame(self.frame4, borderwidth = 0)
        else:
            self.frame4_in = tk.Frame(self.frame4, borderwidth = 0)
            self.frame4.first_view = False
        self.frame4_in.place(relwidth = 0.99, height = 40, relx = 0, y = 50, x = 5, anchor = "nw")

        # P-uncertainties label
        sigmap_label = tk.Label(self.frame4_in, text = "P picking uncertainty (ms)")
        spmin_label = tk.Label(self.frame4_in, text = "Min.")
        spmax_label = tk.Label(self.frame4_in, text = "Max.")

        # P-uncertainty entries
        spmin_entry = tk.Entry(self.frame4_in, textvariable = self.spmin, justify = "right")
        spmax_entry = tk.Entry(self.frame4_in, textvariable = self.spmax, justify = "right")

        # S-uncertainties label
        sigmas_label = tk.Label(self.frame4_in, text = "S picking uncertainty (ms)")
        ssmin_label = tk.Label(self.frame4_in, text = "Min.")
        ssmax_label = tk.Label(self.frame4_in, text = "Max.")

        # S-uncertainty entries
        ssmin_entry = tk.Entry(self.frame4_in, textvariable = self.ssmin, justify = "right")
        ssmax_entry = tk.Entry(self.frame4_in, textvariable = self.ssmax, justify = "right")

        # Layout
        sigmap_label.place(relx = 0, x = 0, y = 0, anchor = "nw")
        spmin_label.place(relx = 0.1, x = -5, y = 22, anchor = "ne")
        spmax_label.place(relx = 0.32, x = -5, y = 22, anchor = "ne")
        spmin_entry.place(relwidth = 0.1, relx = 0.1, x = 0, y = 20, anchor = "nw")
        spmax_entry.place(relwidth = 0.1, relx = 0.32, x = 0, y = 20, anchor = "nw")
        ssmin_entry.place(relwidth = 0.1, relx = 0.6, x = 0, y = 20, anchor = "nw")
        ssmax_entry.place(relwidth = 0.1, relx = 0.82, x = 0, y = 20, anchor = "nw")

        sigmas_label.place(relx = 0.5, x = 0, y = 0, anchor = "nw")
        ssmin_label.place(relx = 0.6, x = -5, y = 22, anchor = "ne")
        ssmax_label.place(relx = 0.82, x = -5, y = 22, anchor = "ne")

    # Widget: Meter
    #===============
    def meter_widget(self):
        if not self.frame5.first_view:
            self.frame5_in.forget()
            self.frame5_in = tk.Frame(self.frame5, borderwidth = 0)
        else:
            self.frame5_in = tk.Frame(self.frame5, borderwidth = 0)
            self.frame5.first_view = False
        self.frame5_in.place(relwidth = 0.99, height = 65, relx = 0, y = 50, x = 5, anchor = "nw")

    # Widget: Degree
    #================
    def degree_widget(self):
        if not self.frame5.first_view:
            self.frame5_in.forget()
            self.frame5_in = tk.Frame(self.frame5, borderwidth = 0)
        else:
            self.frame5_in = tk.Frame(self.frame5, borderwidth = 0)
            self.frame5.first_view = False
        self.frame5_in.place(relwidth = 0.99, height = 65, relx = 0, y = 50, x = 5, anchor = "nw")

        # West
        west_entry = tk.Entry(self.frame5_in, textvariable = self.west, justify = "right")

        # South
        south_entry = tk.Entry(self.frame5_in, textvariable = self.south, justify = "right")

        # East
        east_entry = tk.Entry(self.frame5_in, textvariable = self.east, justify = "right")

        # North
        north_entry = tk.Entry(self.frame5_in, textvariable = self.north, justify = "right")

        # Layout
        west_entry.place(relwidth = 0.15, relx = 0.275, x = 0, y = 20, anchor = "nw")
        south_entry.place(relwidth = 0.15, relx = 0.425, x = 0, y = 40, anchor = "nw")
        east_entry.place(relwidth = 0.15, relx = 0.575, x = 0, y = 20, anchor = "nw")
        north_entry.place(relwidth = 0.15, relx = 0.425, x = 0, anchor = "nw")

    # Widget: Competitive Particle Swarm Optimization
    #=================================================
    def cpso_widget(self):
        if not self.frame6.first_view:
            self.frame6_in.forget()
            self.frame6_in = tk.Frame(self.frame6, borderwidth = 0)
        else:
            self.frame6_in = tk.Frame(self.frame6, borderwidth = 0)
            self.frame6.first_view = False
        self.frame6_in.place(relwidth = 0.99, height = 50, relx = 0, y = 50, x = 5, anchor = "nw")

        # Omega label
        omega_label = tk.Label(self.frame6_in, text = "Inertia")

        # Omega scale
        omega_scale = tk.Scale( self.frame6_in, from_ = 0., to_ = 1., resolution = 0.01, \
                                variable = self.omega, showvalue = 0, \
                                orient = "horizontal", borderwidth = 1, \
                                width = 15, sliderlength = 20, sliderrelief = "ridge" \
                                )

        # Omega entry
        omega_entry = tk.Entry(self.frame6_in, textvariable = self.omega, justify = "right")

        # Phi label
        phi_label = tk.Label(self.frame6_in, text = "Sociability and cognition")

        # Phi scale
        phi_scale = tk.Scale(   self.frame6_in, from_ = 0.5, to_ = 2., resolution = 0.01, \
                                variable = self.phi, showvalue = 0, \
                                orient = "horizontal", borderwidth = 1, \
                                width = 15, sliderlength = 20, sliderrelief = "ridge" \
                                )

        # Phi entry
        phi_entry = tk.Entry(self.frame6_in, textvariable = self.phi, justify = "right")

        # Layout
        omega_label.place(relx = 0, x = 0, y = 0, anchor = "nw")
        omega_scale.place(relwidth = 0.35, relx = 0, x = 0, y = 20, anchor = "nw")
        omega_entry.place(relwidth = 0.1, relx = 0.35, x = -3, y = 21, anchor = "nw")

        phi_label.place(relx = 0.5, x = 0, y = 0, anchor = "nw")
        phi_scale.place(relwidth = 0.35, relx = 0.5, x = 0, y = 20, anchor = "nw")
        phi_entry.place(relwidth = 0.1, relx = 0.85, x = -3, y = 21, anchor = "nw")

    # Widget: Differential Evolution
    #================================
    def de_widget(self):
        if not self.frame6.first_view:
            self.frame6_in.forget()
            self.frame6_in = tk.Frame(self.frame6, borderwidth = 0)
        else:
            self.frame6_in = tk.Frame(self.frame6, borderwidth = 0)
            self.frame6.first_view = False
        self.frame6_in.place(relwidth = 0.99, height = 50, relx = 0, y = 50, x = 5, anchor = "nw")

        # CR label
        CR_label = tk.Label(self.frame6_in, text = "Crossover probability")

        # CR scale
        CR_scale = tk.Scale(    self.frame6_in, from_ = 0., to_ = 1., resolution = 0.01, \
                                variable = self.CR, showvalue = 0, \
                                orient = "horizontal", borderwidth = 1, \
                                width = 15, sliderlength = 20, sliderrelief = "ridge" \
                                )

        # CR entry
        CR_entry = tk.Entry(self.frame6_in, textvariable = self.CR, justify = "right")

        # F label
        F_label = tk.Label(self.frame6_in, text = "Differential weight")

        # F scale
        F_scale = tk.Scale(     self.frame6_in, from_ = 0., to_ = 2., resolution = 0.01, \
                                variable = self.F, showvalue = 0, \
                                orient = "horizontal", borderwidth = 1, \
                                width = 15, sliderlength = 20, sliderrelief = "ridge" \
                                )

        # F entry
        F_entry = tk.Entry(self.frame6_in, textvariable = self.F, justify = "right")

        # Layout
        CR_label.place(relx = 0, x = 0, y = 0, anchor = "nw")
        CR_scale.place(relwidth = 0.35, relx = 0, x = 0, y = 20, anchor = "nw")
        CR_entry.place(relwidth = 0.1, relx = 0.35, x = -3, y = 21, anchor = "nw")

        F_label.place(relx = 0.5, x = 0, y = 0, anchor = "nw")
        F_scale.place(relwidth = 0.35, relx = 0.5, x = 0, y = 20, anchor = "nw")
        F_entry.place(relwidth = 0.1, relx = 0.85, x = -3, y = 21, anchor = "nw")

    # Window: About
    #===============
    def about_window(self):
        about = "EvoLoc 1.0" + "\n" \
                + "An EVOlutionary earthquake LOCation program" + "\n\n" \
                + "Created by Keurfon Luu"
        tkmessage.showinfo("About", about)

    # Window: Close
    #===============
    def close_window(self):
        yes = tkmessage.askyesno("Exit", "Do you really want to quit?")
        if yes: self.close()

    # Function import_parameters
    #============================
    def import_parameters(self):
        filename = tkfile.askopenfilename(   title = "Import parameters", \
                                             filetypes = [ ("Parameters files", ".par") ], \
                                             )
        if len(filename) > 0:
            self.load_parameters(filename)
            if not self.do_bayesh.get():
                self.bayes.invoke()
            else:
                self.bayesh.invoke()
            if self.unit.get() == 0:
                self.meter.invoke()
            else:
                self.degree.invoke()
            if self.method.get() == 1:
                self.cpso.invoke()
            else:
                self.de.invoke()

    # Function export_parameters
    #============================
    def export_parameters(self):
        filename = tkfile.asksaveasfilename( title = "Export parameters", \
                                             filetypes = [ ("Parameters files", ".par") ], \
                                             defaultextension = ".par" \
                                             )
        if len(filename) > 0:
            if self.check_parameters():
                self.save_parameters(filename)

    # Function import_velp_file
    #===========================
    def import_velp_file(self):
        filename = tkfile.askopenfilename(   title = "Import P-velocity model", \
                                             filetypes = [ ("All files", ".*") ] \
                                             )
        if len(filename) > 0: self.velp_filename.set(filename)

    # Function import_vels_file
    #===========================
    def import_vels_file(self):
        filename = tkfile.askopenfilename(   title = "Import S-velocity model", \
                                             filetypes = [ ("All files", ".*") ] \
                                             )
        if len(filename) > 0: self.vels_filename.set(filename)

    # Function view_vp
    #==================
    def view_vp(self):
        if not self.frame2.first_view:
            self.fig.clear()
            self.init_figure()
            self.frame2_in.forget()
            self.frame2_in = tk.Frame(self.frame2, borderwidth = 0)
        else:
            self.frame2_in = tk.Frame(self.frame2, borderwidth = 0)
            self.frame2.first_view = False
        self.frame2_in.place(relwidth = 1, relheight = 1, relx = 0, anchor = "nw")

        self.canvas = FigureCanvasTkAgg(self.fig, master = self.frame2_in)
        self.canvas.show()
        self.canvas.get_tk_widget().pack()
        self.ax1.mouse_init()
        self.velp = self.import_vel(self.velp_filename.get())
        if not self.VelocityError:
            self.view_vel(self.velp)
            self.figure_parameters("P-velocity model")

    # Function view_vs
    #==================
    def view_vs(self):
        if not self.frame2.first_view:
            self.fig.clear()
            self.init_figure()
            self.frame2_in.forget()
            self.frame2_in = tk.Frame(self.frame2, borderwidth = 0)
        else:
            self.frame2_in = tk.Frame(self.frame2, borderwidth = 0)
            self.frame2.first_view = False
        self.frame2_in.place(relwidth = 1, relheight = 1, relx = 0, anchor = "nw")

        self.canvas = FigureCanvasTkAgg(self.fig, master = self.frame2_in)
        self.canvas.show()
        self.canvas.get_tk_widget().pack()
        self.ax1.mouse_init()
        self.vels = self.import_vel(self.vels_filename.get())
        if not self.VelocityError:
            self.view_vel(self.vels)
            self.figure_parameters("S-velocity model")

    # Function import_tobs_file
    #===========================
    def import_tobs_file(self):
        filename = tkfile.askopenfilename(   title = "Import arrival times file", \
                                             filetypes = [ ("All files", ".*") ] \
                                             )
        if len(filename) > 0: self.tobs_filename.set(filename)

    # Function import_stat_file
    #===========================
    def import_stat_file(self):
        filename = tkfile.askopenfilename(   title = "Import stations file", \
                                             filetypes = [ ("All files", ".*") ] \
                                             )
        if len(filename) > 0: self.stat_filename.set(filename)

    # Function import_sigma_file
    #============================
    def import_sigma_file(self):
        filename = tkfile.askopenfilename(   title = "Import uncertainties file", \
                                             filetypes = [ ("All files", ".*") ] \
                                             )
        if len(filename) > 0: self.sigma_filename.set(filename)

    # Function open_outdir
    #======================
    def open_outdir(self):
        dirname = tkfile.askdirectory(title = "Open output directory")
        if len(dirname) > 0: self.outdir.set(dirname)

    # Function define_variables
    #===========================
    def define_variables(self):
        self.velp_filename = tk.StringVar(self.master)
        self.vels_filename = tk.StringVar(self.master)

        self.nz = tk.Variable(self.master)
        self.nx = tk.Variable(self.master)
        self.ny = tk.Variable(self.master)
        self.dz = tk.Variable(self.master)
        self.dx = tk.Variable(self.master)
        self.dy = tk.Variable(self.master)

        self.stat_filename = tk.StringVar(self.master)
        self.tobs_filename = tk.StringVar(self.master)
        self.refdate = tk.StringVar(self.master)
        self.zcoord = tk.StringVar(self.master)

        self.do_bayesh = tk.BooleanVar(self.master)
        self.sigma_filename = tk.StringVar(self.master)
        self.spmin = tk.Variable(self.master)
        self.spmax = tk.Variable(self.master)
        self.ssmin = tk.Variable(self.master)
        self.ssmax = tk.Variable(self.master)

        self.unit = tk.IntVar(self.master)
        self.west = tk.Variable(self.master)
        self.south = tk.Variable(self.master)
        self.east = tk.Variable(self.master)
        self.north = tk.Variable(self.master)

        self.method = tk.IntVar(self.master)
        self.popsize = tk.IntVar(self.master)
        self.itermax = tk.IntVar(self.master)
        self.omega = tk.Variable(self.master)
        self.phi = tk.Variable(self.master)
        self.CR = tk.Variable(self.master)
        self.F = tk.Variable(self.master)

        self.outdir = tk.StringVar(self.master)
        self.pdf = tk.BooleanVar(self.master)

        self.num_threads = tk.IntVar(self.master)

    # Function trace_variables
    #==========================
    def trace_variables(self):
        self.velp_filename.trace("w", self.callback)
        self.vels_filename.trace("w", self.callback)

        self.nz.trace("w", self.callback)
        self.nx.trace("w", self.callback)
        self.ny.trace("w", self.callback)
        self.dz.trace("w", self.callback)
        self.dx.trace("w", self.callback)
        self.dy.trace("w", self.callback)

        self.stat_filename.trace("w", self.callback)
        self.tobs_filename.trace("w", self.callback)
        self.refdate.trace("w", self.callback)
        self.zcoord.trace("w", self.callback)

        self.do_bayesh.trace("w", self.callback)
        self.sigma_filename.trace("w", self.callback)
        self.spmin.trace("w", self.callback)
        self.spmax.trace("w", self.callback)
        self.ssmin.trace("w", self.callback)
        self.ssmax.trace("w", self.callback)

        self.unit.trace("w", self.callback)
        self.west.trace("w", self.callback)
        self.south.trace("w", self.callback)
        self.east.trace("w", self.callback)
        self.north.trace("w", self.callback)

        self.method.trace("w", self.callback)
        self.popsize.trace("w", self.callback)
        self.itermax.trace("w", self.callback)
        self.omega.trace("w", self.callback)
        self.phi.trace("w", self.callback)
        self.CR.trace("w", self.callback)
        self.F.trace("w", self.callback)

        self.outdir.trace("w", self.callback)
        self.pdf.trace("w", self.callback)

        self.num_threads.trace("w", self.callback)

    # Function init_variables
    #=========================
    def init_variables(self):
        self.nz.set(0)
        self.nx.set(0)
        self.ny.set(0)
        self.dz.set(0.)
        self.dx.set(0.)
        self.dy.set(0.)
        self.refdate.set("")
        self.zcoord.set("Depth")
        self.do_bayesh.set(False)
        self.spmin.set(0.)
        self.spmax.set(0.)
        self.ssmin.set(0.)
        self.ssmax.set(0.)
        self.unit.set(0)
        self.west.set(0.)
        self.south.set(0.)
        self.east.set(0.)
        self.north.set(0.)
        self.method.set(1)
        self.popsize.set(20)
        self.itermax.set(100)
        self.omega.set(0.7)
        self.phi.set(1.5)
        self.CR.set(0.5)
        self.F.set(1.)
        self.outdir.set("./output/")
        self.pdf.set(False)
        self.num_threads.set(1)

    # Function check_parameters
    #===========================
    def check_parameters(self):
        # Check nz
        if int(self.nz.get()) < 1:
            tkmessage.showerror("Error", "nz < 1.")
            return False

        # Check nx
        if int(self.nx.get()) < 1:
            tkmessage.showerror("Error", "nx < 1.")
            return False

        # Check ny
        if int(self.ny.get()) < 1:
            tkmessage.showerror("Error", "ny < 1.")
            return False

        # Check dz
        if float(self.dz.get()) <= 0.:
            tkmessage.showerror("Error", "dz <= 0.0.")
            return False

        # Check dx
        if float(self.dx.get()) <= 0.:
            tkmessage.showerror("Error", "dx <= 0.0.")
            return False

        # Check dy
        if float(self.dy.get()) <= 0.:
            tkmessage.showerror("Error", "dy <= 0.0.")
            return False

        # Check P-velocity model
        try:
            vel = np.fromfile(self.velp_filename.get(), dtype = "float32")
            if (int(self.nz.get())-1)*(int(self.nx.get())-1)*(int(self.ny.get())-1) != len(vel):
                tkmessage.showerror("Error", "An error occurred when checking " + self.velp_filename.get() \
                                    + "\nPlease check the grid size ((nz-1)*(nx-1)*(ny-1) should be equal to " \
                                    + str(len(vel)) + ").")
                return False
        except FileNotFoundError:
            tkmessage.showerror("Error", "P-velocity model file not found.")
            return False

        # Check S-velocity model
        try:
            vel = np.fromfile(self.vels_filename.get(), dtype = "float32")
            if (int(self.nz.get())-1)*(int(self.nx.get())-1)*(int(self.ny.get())-1) != len(vel):
                tkmessage.showerror("Error", "An error occurred when checking " + self.vels_filename.get() \
                                    + "\nPlease check the grid size ((nz-1)*(nx-1)*(ny-1) should be equal to " \
                                    + str(len(vel)) + ").")
                return False
        except FileNotFoundError:
            tkmessage.showerror("Error", "S-velocity model file not found.")
            return False

        # Check stations file
        try:
            rcv = np.loadtxt(self.stat_filename.get())
            if rcv.shape[1] != 3:
                tkmessage.showerror("Error", "Stations file has " + str(rcv.shape[1]) + " columns instead of 3.")
                return False
        except FileNotFoundError:
            tkmessage.showerror("Error", "Stations file not found.")
            return False
        except ValueError:
            tkmessage.showerror("Error", "Stations file is inadequate.")
            return False
        nstat = rcv.shape[0]

        # Check arrival times file
        try:
            tobs = np.loadtxt(self.tobs_filename.get())
            if tobs.shape[1] != 2:
                tkmessage.showerror("Error", "Arrival times file has " + str(tobs.shape[1]) + " columns instead of 2.")
                return False
        except FileNotFoundError:
            tkmessage.showerror("Error", "Arrival times file not found.")
            return False
        except ValueError:
            tkmessage.showerror("Error", "Arrival times file is inadequate.")
            return False
        nobs = tobs.shape[0]

        # Check stations file and arrival times file
        if nobs % nstat > 0:
            tkmessage.showerror("Error", "Stations and arrival times files are unconsistent. " \
                                + "The number of observations does not match the number of stations.")
            return False

        # Check the reference date format
        if self.refdate.get() != "":
            try:
                tmp = datetime.strptime(self.refdate.get(), "%Y-%m-%d")
            except ValueError:
                tkmessage.showerror("Error", "Time data " + self.refdate.get() + " does not match format YYYY-MM-DD.")
                return False

        # Check uncertainties file
        if not self.do_bayesh.get():
            try:
                sigma = np.loadtxt(self.sigma_filename.get())
                if sigma.shape[1] != 2:
                    tkmessage.showerror("Error", "Uncertainties file has " + str(sigma.shape[1]) + " columns instead of 2.")
                    return False
            except FileNotFoundError:
                tkmessage.showerror("Error", "Uncertainties file not found.")
                return False
            except ValueError:
                tkmessage.showerror("Error", "Uncertainties file is inadequate.")
                return False

            if sigma.shape[0] != nobs:
                tkmessage.showerror("Error", "Arrival times and uncertainties files are unconsistent. " \
                                + "The number of observations does not match the number of available uncertainties.")
                return False

        # Check uncertainty boundaries
        else:
            # Check spmin
            if float(self.spmin.get()) <= 0.:
                tkmessage.showerror("Error", "spmin <= 0.0.")
                return False

            # Check spmax
            if float(self.spmax.get()) <= 0.:
                tkmessage.showerror("Error", "spmax <= 0.0.")
                return False
            elif float(self.spmax.get()) < float(self.spmin.get()):
                tkmessage.showerror("Error", "spmax < spmin.")
                return False

            # Check ssmin
            if float(self.ssmin.get()) <= 0.:
                tkmessage.showerror("Error", "ssmin <= 0.0.")
                return False

            # Check ssmax
            if float(self.ssmax.get()) <= 0.:
                tkmessage.showerror("Error", "ssmax <= 0.0.")
                return False
            elif float(self.ssmax.get()) < float(self.ssmin.get()):
                tkmessage.showerror("Error", "ssmax < ssmin.")
                return False

        # Check coordinates
        if int(self.unit.get()) == 1:
            # Check west
            if float(self.west.get()) < -180.:
                tkmessage.showerror("Error", "west < -180.0.")
                return False

            # Check south
            if float(self.south.get()) < -90.:
                tkmessage.showerror("Error", "south < -90.0.")
                return False

            # Check east
            if float(self.east.get()) > 180.:
                tkmessage.showerror("Error", "east > 180.0.")
                return False
            elif float(self.east.get()) <= float(self.west.get()):
                tkmessage.showerror("Error", "east <= west.")
                return False

            # Check north
            if float(self.north.get()) > 90.:
                tkmessage.showerror("Error", "north > 90.0.")
                return False
            elif float(self.north.get()) <= float(self.south.get()):
                tkmessage.showerror("Error", "north <= south.")
                return False

        return True

    # Function save_parameters
    #==========================
    def save_parameters(self, filename):
        fid = open(filename, "w")

        fid.write("# Grid size \n")
        fid.write(str(self.nz.get()) + " " + str(self.nx.get()) + " " + str(self.ny.get()) + "\n")
        fid.write("# Mesh size \n")
        fid.write(str(self.dz.get()) + " " + str(self.dx.get()) + " " + str(self.dy.get()) + "\n")
        fid.write("# P-velocity model \n")
        fid.write('"' + self.velp_filename.get() + '"' + "\n")
        fid.write("# S-velocity model \n")
        fid.write('"' + self.vels_filename.get() + '"' + "\n")

        fid.write("# Stations file \n")
        fid.write('"' + self.stat_filename.get() + '"' + "\n")
        fid.write("# Arrival times file \n")
        fid.write('"' + self.tobs_filename.get() + '"' + "\n")
        fid.write("# Reference date used to create data (year / month / day) \n")
        if self.refdate.get() == "":
            fid.write("0 0 0" + "\n")
        else:
            tmp = datetime.strptime(self.refdate.get(), "%Y-%m-%d")
            fid.write(str(tmp.year) + " " + str(tmp.month) + " " + str(tmp.day) + "\n")
        fid.write("# Z coordinates of stations (0: depth, 1: elevation) \n")
        if self.zcoord.get() == "Depth":
            fid.write(str(0) + "\n")
        elif self.zcoord.get() == "Elevation":
            fid.write(str(1) + "\n")

        fid.write("# Probabilistic formulation \n")
        fid.write(str(int(self.do_bayesh.get())) + "\n")
        fid.write("# Uncertainties file \n")
        fid.write('"' + self.sigma_filename.get() + '"' + "\n")
        fid.write("# P picking uncertainty boundaries (in ms, spmin / spmax) \n")
        fid.write(str(self.spmin.get()) + " " + str(self.spmax.get()) + "\n")
        fid.write("# S picking uncertainty boundaries (in ms, ssmin / ssmax) \n")
        fid.write(str(self.ssmin.get()) + " " + str(self.ssmax.get()) + "\n")

        fid.write("# Coordinates unit (0: meter, 1: degree) \n")
        fid.write(str(self.unit.get()) + "\n")
        fid.write("# Coordinate boundaries if degree (west / south / east / north) \n")
        fid.write(str(self.west.get()) + " " + str(self.south.get()) + " " + str(self.east.get()) + " " + str(self.north.get()) + "\n")

        fid.write("# Optimization method (1: CPSO, 2: DE) \n")
        fid.write(str(self.method.get()) + "\n")
        fid.write("# Population size \n")
        fid.write(str(self.popsize.get()) + "\n")
        fid.write("# Maximum number of iterations \n")
        fid.write(str(self.itermax.get()) + "\n")
        fid.write("# CPSO parameters (omega / phi) \n")
        fid.write(str(self.omega.get()) + " " + str(self.phi.get()) + "\n")
        fid.write("# DE parameters (CR / F) \n")
        fid.write(str(self.CR.get()) + " " + str(self.F.get()) + "\n")

        fid.write("# Output directory \n")
        fid.write('"' + self.outdir.get() + '"' + "\n")
        fid.write("# Compute PDF\n")
        fid.write(str(int(self.pdf.get())) + "\n")

        fid.write("# Number of threads \n")
        fid.write(str(self.num_threads.get()) + "\n")

        fid.close()

    # Function load_parameters
    #==========================
    def load_parameters(self, filename):
        with open(filename) as fid:
            par = fid.readlines()

        tmp = [ int(each) for each in par[1].split() ]
        self.nz.set(tmp[0])
        self.nx.set(tmp[1])
        self.ny.set(tmp[2])
        tmp = [ float(each) for each in par[3].split() ]
        self.dz.set(tmp[0])
        self.dx.set(tmp[1])
        self.dy.set(tmp[2])
        tmp = par[5].split("\n")[0][1:-1].replace('"', "")
        self.velp_filename.set(tmp)
        tmp = par[7].split("\n")[0][1:-1].replace('"', "")
        self.vels_filename.set(tmp)

        tmp = par[9].split("\n")[0][1:-1].replace('"', "")
        self.stat_filename.set(tmp)
        tmp = par[11].split("\n")[0][1:-1].replace('"', "")
        self.tobs_filename.set(tmp)
        tmp = [ int(each) for each in par[13].split() ]
        if tmp == [ 0, 0, 0 ]:
            self.refdate.set("")
        else:
            self.refdate.set(str(tmp[0]).zfill(4) + "-" + str(tmp[1]).zfill(2) + "-" + str(tmp[2]).zfill(2))
        tmp = int(par[15].split()[0])
        if tmp == 0:
            self.zcoord.set("Depth")
        elif tmp == 1:
            self.zcoord.set("Elevation")

        tmp = bool(int(par[17].split()[0]))
        self.do_bayesh.set(tmp)
        tmp = par[19].split("\n")[0][1:-1].replace('"', "")
        self.sigma_filename.set(tmp)
        tmp = [ float(each) for each in par[21].split() ]
        self.spmin.set(tmp[0])
        self.spmax.set(tmp[1])
        tmp = [ float(each) for each in par[23].split() ]
        self.ssmin.set(tmp[0])
        self.ssmax.set(tmp[1])

        tmp = int(par[25].split()[0])
        self.unit.set(tmp)
        tmp = [ float(each) for each in par[27].split() ]
        self.west.set(tmp[0])
        self.south.set(tmp[1])
        self.east.set(tmp[2])
        self.north.set(tmp[3])

        tmp = int(par[29].split()[0])
        self.method.set(tmp)
        tmp = int(par[31].split()[0])
        self.popsize.set(tmp)
        tmp = int(par[33].split()[0])
        self.itermax.set(tmp)
        tmp = [ float(each) for each in par[35].split() ]
        self.omega.set(tmp[0])
        self.phi.set(tmp[1])
        tmp = [ float(each) for each in par[37].split() ]
        self.CR.set(tmp[0])
        self.F.set(tmp[1])

        tmp = par[39].split("\n")[0][1:-1].replace('"', "")
        self.outdir.set(tmp)
        tmp = int(par[41].split()[0])
        self.pdf.set(tmp)

        tmp = int(par[43].split()[0])
        self.num_threads.set(min(tmp, multiprocessing.cpu_count()))

    # Function init_figure
    #======================
    def init_figure(self):
        self.fig = Figure(figsize = (12, 8), facecolor = "white", dpi = 100)
        self.ax1 = self.fig.add_subplot(1, 1, 1, projection = "3d")

    # Function figure_axis
    #======================
    def figure_axis(self):
        az = np.linspace(0, self.nz.get()-1, self.nz.get()-1)
        ax = np.linspace(0, self.nx.get()-1, self.nx.get()-1)
        ay = np.linspace(0, self.ny.get()-1, self.ny.get()-1)
        return az, ax, ay

    # Function figure_parameters
    #============================
    def figure_parameters(self, title):
        self.ax1.set_title(title)
        self.ax1.set_xlabel("X (km)")
        self.ax1.set_ylabel("Y (km)")
        self.ax1.set_zlabel("Z (km)")
        xmin, xmax = 0., float(self.nx.get()) * float(self.dx.get()) * 0.001
        ymin, ymax = 0., float(self.ny.get()) * float(self.dy.get()) * 0.001
        zmin, zmax = 0., float(self.nz.get()) * float(self.dz.get()) * 0.001
        self.ax1.set_xticks([ xmin, xmax ])
        self.ax1.set_yticks([ ymin, ymax ])
        self.ax1.set_zticks([ zmin, zmax ])
        self.ax1.set_xlim(xmin, xmax)
        self.ax1.set_ylim(ymin, ymax)
        self.ax1.set_zlim(zmin, zmax)
        self.ax1.invert_zaxis()
        self.ax1.view_init(elev = 20, azim = 240)

    # Function import_vel
    #=====================
    def import_vel(self, filename):
        nz = int(self.nz.get()) - 1
        nx = int(self.nx.get()) - 1
        ny = int(self.ny.get()) - 1

        try:
            vel = np.fromfile(filename, dtype = "float32")
            try:
                vel = np.reshape(vel, (nz, nx, ny), order = "F")
                self.VelocityError = False
                return vel
            except ValueError:
                tkmessage.showerror("Error", "An error occurred when loading " + filename \
                                    + "\nPlease check the grid size ((nz-1)*(nx-1)*(ny-1) should be equal to " \
                                    + str(len(vel)) + ").")
                self.VelocityError = True
        except FileNotFoundError:
            tkmessage.showerror("Error", "Velocity model file not found.")
            self.VelocityError = True

    # Function view_vel
    #===================
    def view_vel(self, vel):
        nz = int(self.nz.get()) - 1
        nx = int(self.nx.get()) - 1
        ny = int(self.ny.get()) - 1
        dz = float(self.dz.get()) * 0.001
        dx = float(self.dx.get()) * 0.001
        dy = float(self.dy.get()) * 0.001
        az = np.linspace(0, nz*dz, nz)
        ax = np.linspace(0, nx*dx, nx)
        ay = np.linspace(0, ny*dy, ny)
        vel *= 0.001
        vmin = np.min(vel)
        vmax = np.max(vel)
        X, Y = np.meshgrid(ax, ay)
        self.ax1.contourf(X, Y, vel[1,:,:].T, 50, vmin = vmin, vmax = vmax, zdir = "z", offset = 0, cmap = "jet")
        X, Z = np.meshgrid(ax, az)
        self.ax1.contourf(X, vel[:,:,0], Z, 50, vmin = vmin, vmax = vmax, zdir = "y", offset = 0, cmap = "jet")
        Y, Z = np.meshgrid(ay, az)
        im = self.ax1.contourf(vel[:,-1,:], Y, Z, 50, vmin = vmin, vmax = vmax, zdir = "x", offset = 0, cmap = "jet")
        self.cb = self.fig.colorbar(im)
        self.cb.set_clim(vmin, vmax)
        self.cb.set_label("Velocity (km/s)")
        self.cb.ax.tick_params(labelsize = 11)

    # Function locate
    #=================
    def locate(self):
        if self.check_parameters():
            if not os.path.exists(".evotmp"):
                os.makedirs(".evotmp")
            fid = open(".evotmp/signal_emit", "wb")
            fid.write(np.array([ 1. ], dtype = "float32").tostring())
            fid.close()
            self.save_parameters(".evotmp/parameters")
            self.master.quit()
            self.master.destroy()
        else:
            return

    # Function close
    #================
    def close(self):
        if not os.path.exists(".evotmp"):
            os.makedirs(".evotmp")
        fid = open(".evotmp/signal_emit", "wb")
        fid.write(np.array([ 0. ], dtype = "float32").tostring())
        fid.close()
        self.master.quit()
        self.master.destroy()

    # Function callback
    #===================
    def callback(self, *args):
        return

root = tk.Tk()
gui = evolocGUI(root)
root.mainloop()
