__author__ = 'szieger'
__project__ = 'in silico sensor response'

# coding: utf-8
import basics_sensorSignal as bs
import sys
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
import matplotlib.pyplot as plt
import seaborn as sns
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QMainWindow, QPushButton, QFileDialog,
                              QAction, qApp, QGridLayout, QLabel, QInputDialog, QLineEdit, QCheckBox, QTextEdit,
                              QGroupBox, QMessageBox, QTableWidget, QTableWidgetItem, QFrame)
from PyQt5.QtGui import QIcon, QValidator, QDoubleValidator, QIntValidator, QColor
from PyQt5.QtCore import Qt, QSize
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT, FigureCanvasQTAgg
import pandas as pd
import numpy as np
import seaborn as sns
import os
from matplotlib.patches import Ellipse
import matplotlib.patches as mpatches
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

fs = 10.


class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.initUI()
        self.xcoords = []
        self.setWindowIcon(QIcon('icon.png'))

    def initUI(self):
        # creating main window
        w = QWidget()
        self.setCentralWidget(w)
        self.setWindowTitle('Sensor response')

        # Menu bar - Load data, Save report, Save all, Exit
        loadAction = QAction('&Load data', self)
        loadAction.setStatusTip('Load data')
        saveAction_report = QAction('&Save report', self)
        saveAction_report.setStatusTip('Save results')
        saveAction_all = QAction('&Save all', self)
        saveAction_all.setStatusTip('Save all')
        exitAction = QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(qApp.quit)
        self.statusBar()

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(loadAction)
        fileMenu.addAction(saveAction_report)
        fileMenu.addAction(saveAction_all)
        fileMenu.addAction(exitAction)

        # (Invisible) structure of main window (=grid)
        mlayout = QVBoxLayout(w)

        # 1st layer: box with vertical alignment to split in top and bottom
        vbox_top, vbox_middle, vbox_bottom = QHBoxLayout(), QHBoxLayout(), QHBoxLayout()
        mlayout.addLayout(vbox_top), mlayout.addLayout(vbox_middle), mlayout.addLayout(vbox_bottom)

        # 2nd layer: box with vertical alignment to split in top and bottom
        hbox_left, hbox_m1, hbox_middle, hbox_right = QVBoxLayout(), QVBoxLayout(), QVBoxLayout(), QVBoxLayout()
        vbox_top.addLayout(hbox_left), vbox_top.addLayout(hbox_m1), vbox_top.addLayout(hbox_middle)
        vbox_top.addLayout(hbox_right)

        # 3rd layer: left box with parameter settings
        hbox_ltop, hbox_lmiddle, hbox_lbottom = QHBoxLayout(), QHBoxLayout(), QHBoxLayout()
        hbox_ltop.setContentsMargins(5, 10, 50, 5)
        hbox_left.addLayout(hbox_ltop), hbox_left.addLayout(hbox_lmiddle), hbox_left.addLayout(hbox_lbottom)

        # 4th layer: middle and right box with plotting sections
        hbox_mtop, hbox_mbottom, hbox_rtop, hbox_rbottom = QHBoxLayout(), QHBoxLayout(), QHBoxLayout(), QHBoxLayout()
        hbox_middle.addLayout(hbox_mtop), hbox_middle.addLayout(hbox_mbottom)
        hbox_right.addLayout(hbox_rtop), hbox_right.addLayout(hbox_rbottom)

    # ----------------------------------------------------
        # draw additional "line" to separate parameters from plots and to separate navigation from rest
        vline = QFrame()
        vline.setFrameShape(QFrame.VLine | QFrame.Raised)
        vline.setLineWidth(2)
        hbox_m1.addWidget(vline)

        hline = QFrame()
        hline.setFrameShape(QFrame.HLine | QFrame.Raised)
        hline.setLineWidth(2)
        vbox_middle.addWidget(hline)

        # ------------------------------------------------------------------------------------
        # left side of main window (-> data treatment)
        vbox_top.addWidget(w)
        vbox_top.setContentsMargins(5, 10, 10, 5)

    # ---------------------------------------------------
        # PARAMETERS
        # general settings
        temperature_label, temperature_unit_label = QLabel(self), QLabel(self)
        temperature_label.setText('Temperature')
        temperature_unit_label.setText('degC')
        self.temperature_edit = QLineEdit(self)
        self.temperature_edit.setValidator(QDoubleValidator())
        self.temperature_edit.setAlignment(Qt.AlignRight)
        self.temperature_edit.setText('15.')

        pH_label = QLabel(self)
        pH_label.setText('pH range')
        self.phrange_edit = QLineEdit(self)
        self.phrange_edit.setValidator(QRegExpValidator())
        self.phrange_edit.setAlignment(Qt.AlignRight)
        self.phrange_edit.setText('0, 14, 0.01')

        smprate_label, smprate_unit_label = QLabel(self), QLabel(self)
        smprate_label.setText('sampling rate')
        smprate_unit_label.setText('s')
        self.smpgrate_edit = QLineEdit(self)
        self.smpgrate_edit.setValidator(QDoubleValidator())
        self.smpgrate_edit.setAlignment(Qt.AlignRight)
        self.smpgrate_edit.setText('1.')

        # pH sensor settings
        ph_t90_label, ph_t90_unit = QLabel(self), QLabel(self)
        ph_t90_label.setText('Response time')
        ph_t90_unit.setText('s')
        self.ph_t90_edit = QLineEdit(self)
        self.ph_t90_edit.setValidator(QDoubleValidator())
        self.ph_t90_edit.setAlignment(Qt.AlignRight)
        self.ph_t90_edit.setText('1e-2')

        ph_signal_label, ph_signal_unit = QLabel(self), QLabel(self)
        ph_signal_label.setText('Sensor potential')
        ph_signal_unit.setText('mV')
        self.ph_signal_edit = QLineEdit(self)
        self.ph_signal_edit.setValidator(QRegExpValidator())
        self.ph_signal_edit.setAlignment(Qt.AlignRight)
        self.ph_signal_edit.setText('5.0, 400.')

        ph_res_label, ph_res_unit = QLabel(self), QLabel(self)
        ph_res_label.setText('Sensor resolution')
        ph_res_unit.setText('mV')
        self.ph_res_edit = QLineEdit(self)
        self.ph_res_edit.setValidator(QDoubleValidator())
        self.ph_res_edit.setAlignment(Qt.AlignRight)
        self.ph_res_edit.setText('1e-5')

        ph_ref_label, ph_ref_unit = QLabel(self), QLabel(self)
        ph_ref_label.setText('Pot. reference electrode')
        ph_ref_unit.setText('mV')
        self.ph_ref_edit = QLineEdit(self)
        self.ph_ref_edit.setValidator(QDoubleValidator())
        self.ph_ref_edit.setAlignment(Qt.AlignRight)
        self.ph_ref_edit.setText('5e-3')

        # NH3 sensor settings
        nh3_t90_label, nh3_t90_unit = QLabel(self), QLabel(self)
        nh3_t90_label.setText('Response time')
        nh3_t90_unit.setText('s')
        self.nh3_t90_edit = QLineEdit(self)
        self.nh3_t90_edit.setValidator(QDoubleValidator())
        self.nh3_t90_edit.setAlignment(Qt.AlignRight)
        self.nh3_t90_edit.setText('0.1')

        nh3_signal_label, nh3_signal_unit = QLabel(self), QLabel(self)
        nh3_signal_label.setText('Sensor potential')
        nh3_signal_unit.setText('mV')
        self.nh3_signal_edit = QLineEdit(self)
        self.nh3_signal_edit.setValidator(QRegExpValidator())
        self.nh3_signal_edit.setAlignment(Qt.AlignRight)
        self.nh3_signal_edit.setText('0.02, 0.09')

        nh3_res_label, nh3_res_unit = QLabel(self), QLabel(self)
        nh3_res_label.setText('Sensor resolution')
        nh3_res_unit.setText('mV')
        self.nh3_res_edit = QLineEdit(self)
        self.nh3_res_edit.setValidator(QDoubleValidator())
        self.nh3_res_edit.setAlignment(Qt.AlignRight)
        self.nh3_res_edit.setText('1e-9')

        nh3_pka_label = QLabel(self)
        nh3_pka_label.setText('pKa')
        self.nh3_pka_edit = QLineEdit(self)
        self.nh3_pka_edit.setValidator(QDoubleValidator())
        self.nh3_pka_edit.setAlignment(Qt.AlignRight)
        self.nh3_pka_edit.setText('9.25')

        nh3_cGG_label, nh3_cGG_unit = QLabel(self), QLabel(self)
        nh3_cGG_label.setText('c(NH3) at pKa')
        nh3_cGG_unit.setText('ppm')
        self.nh3_cGG_edit = QLineEdit(self)
        self.nh3_cGG_edit.setValidator(QDoubleValidator())
        self.nh3_cGG_edit.setAlignment(Qt.AlignRight)
        self.nh3_cGG_edit.setText('230.')

        nh3_alpha_label, nh3_alpha_unit = QLabel(self), QLabel(self)
        nh3_alpha_label.setText('alpha(NH3)')
        nh3_alpha_unit.setText('%')
        self.nh3_alpha_edit = QLineEdit(self)
        self.nh3_alpha_edit.setValidator(QRegExpValidator())
        self.nh3_alpha_edit.setAlignment(Qt.AlignRight)
        self.nh3_alpha_edit.setText('0., 100., 0.01')

        # -----------------------
        # General navigation
        self.plot_button = QPushButton('Plot', self)
        self.plot_button.setFixedWidth(100)
        self.clear_button = QPushButton('Clear', self)
        self.clear_button.setFixedWidth(100)
        self.save_button = QPushButton('Save', self)
        self.save_button.setFixedWidth(100)

        # -------------------------------------------------------------------------------------------
        # GroupBoxes to structure the layout
        # General navigation
        navigation_group = QGroupBox("Navigation Tool Bar")
        grid_load = QGridLayout()
        navigation_group.setFixedHeight(70)
        vbox_bottom.addWidget(navigation_group)
        navigation_group.setLayout(grid_load)

        grid_load.addWidget(self.plot_button, 0, 0)
        grid_load.addWidget(self.clear_button, 0, 1)
        grid_load.addWidget(self.save_button, 0, 2)

        # ----------------------------------------------
        # create GroupBox to structure the layout
        general_group = QGroupBox("General Settings")
        general_group.setFixedWidth(250)
        grid_load = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_ltop.addWidget(general_group)
        general_group.setLayout(grid_load)
        grid_load.addWidget(temperature_label, 0, 0)
        grid_load.addWidget(self.temperature_edit, 0, 1)
        grid_load.addWidget(temperature_unit_label, 0, 2)
        grid_load.addWidget(pH_label, 1, 0)
        grid_load.addWidget(self.phrange_edit, 1, 1)
        grid_load.addWidget(smprate_label, 2, 0)
        grid_load.addWidget(self.smpgrate_edit, 2, 1)
        grid_load.addWidget(smprate_unit_label, 2, 2)

        general_group.setContentsMargins(1, 15, 15, 1)
        hbox_ltop.addSpacing(10)

    # -----------------------
        # pH Sensor Settings
        phsens_group = QGroupBox("pH Sensor Settings")
        phsens_group.setFixedWidth(300)
        grid_load = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_lmiddle.addWidget(phsens_group)
        phsens_group.setLayout(grid_load)

        grid_load.addWidget(ph_t90_label, 0, 0)
        grid_load.addWidget(self.ph_t90_edit, 0, 1)
        grid_load.addWidget(ph_t90_unit, 0, 2)
        grid_load.addWidget(ph_signal_label, 1, 0)
        grid_load.addWidget(self.ph_signal_edit, 1, 1)
        grid_load.addWidget(ph_signal_unit, 1, 2)
        grid_load.addWidget(ph_res_label, 2, 0)
        grid_load.addWidget(self.ph_res_edit, 2, 1)
        grid_load.addWidget(ph_res_unit, 2, 2)
        grid_load.addWidget(ph_ref_label, 3, 0)
        grid_load.addWidget(self.ph_ref_edit, 3, 1)
        grid_load.addWidget(ph_ref_unit, 3, 2)

        phsens_group.setContentsMargins(1, 15, 15, 1)
        hbox_lmiddle.addSpacing(10)

    # -----------------------
        # NH3 Sensor Settings
        nh3sens_group = QGroupBox("NH3 Sensor Settings")
        nh3sens_group.setFixedWidth(300)
        grid_load = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_lbottom.addWidget(nh3sens_group)
        nh3sens_group.setLayout(grid_load)

        grid_load.addWidget(nh3_t90_label, 0, 0)
        grid_load.addWidget(self.nh3_t90_edit, 0, 1)
        grid_load.addWidget(nh3_t90_unit, 0, 2)
        grid_load.addWidget(nh3_signal_label, 1, 0)
        grid_load.addWidget(self.nh3_signal_edit, 1, 1)
        grid_load.addWidget(nh3_signal_unit, 1, 2)
        grid_load.addWidget(nh3_res_label, 2, 0)
        grid_load.addWidget(self.nh3_res_edit, 2, 1)
        grid_load.addWidget(nh3_res_unit, 2, 2)
        grid_load.addWidget(nh3_pka_label, 3, 0)
        grid_load.addWidget(self.nh3_pka_edit, 3, 1)
        grid_load.addWidget(nh3_cGG_label, 4, 0)
        grid_load.addWidget(self.nh3_cGG_edit, 4, 1)
        grid_load.addWidget(nh3_cGG_unit, 4, 2)
        grid_load.addWidget(nh3_alpha_label, 5, 0)
        grid_load.addWidget(self.nh3_alpha_edit, 5, 1)
        grid_load.addWidget(nh3_alpha_unit, 5, 2)

        nh3sens_group.setContentsMargins(1, 15, 15, 1)
        hbox_lbottom.addSpacing(10)

        # for all parameters - connect LineEdit with function
        self.temperature_edit.editingFinished.connect(self.print_temperature)
        self.phrange_edit.editingFinished.connect(self.print_phrange)
        self.smpgrate_edit.editingFinished.connect(self.print_samplingrate)
        self.ph_t90_edit.editingFinished.connect(self.print_ph_t90)
        self.ph_signal_edit.editingFinished.connect(self.print_ph_signal)
        self.ph_res_edit.editingFinished.connect(self.print_ph_resolution)
        self.ph_ref_edit.editingFinished.connect(self.print_ph_reference_pot)
        self.nh3_t90_edit.editingFinished.connect(self.print_nh3_t90)
        self.nh3_signal_edit.editingFinished.connect(self.print_nh3_signal)
        self.nh3_res_edit.editingFinished.connect(self.print_nh3_resolution)
        self.nh3_pka_edit.editingFinished.connect(self.print_nh3_pka)
        self.nh3_cGG_edit.editingFinished.connect(self.print_nh3_concentration)
        self.nh3_alpha_edit.editingFinished.connect(self.print_nh3_alpha)

    # -------------------------------------------------------------------
        # connect button with function
        # !!!TODO: connect button with function - plot, clear, and save

    # --------------------------------------------------------------------------------------------
        # Calibration of pH
        self.fig_phcalib, self.ax_phcalib = plt.subplots()
        self.canvas_phcalib = FigureCanvasQTAgg(self.fig_phcalib)
        self.navi_phcalib = NavigationToolbar2QT(self.canvas_phcalib, w, coordinates=False)
        self.ax_phcalib.set_xlim(0, 15)
        self.ax_phcalib.set_xlabel('pH value', fontsize=fs)
        self.ax_phcalib.set_ylabel('Potential [mV]', fontsize=fs)
        self.fig_phcalib.subplots_adjust(left=0.15, right=0.95, bottom=0.29, top=0.85)
        sns.despine()

        # connect onclick event with function
    #     self.fig_phcalib.canvas.mpl_connect('button_press_event', self.onclick_timedrive)

        # create GroupBox to structure the layout
        phcalib_group = QGroupBox("pH Sensor Calibration")
        phcalib_group.setMinimumWidth(200)
        grid_phcalib = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_mtop.addWidget(phcalib_group)
        phcalib_group.setLayout(grid_phcalib)
        grid_phcalib.addWidget(self.canvas_phcalib)
        grid_phcalib.addWidget(self.navi_phcalib)

        # ---------------------------------------------------
        # Calibration of NH3
        self.fig_nh3calib, self.ax_nh3calib = plt.subplots()
        self.canvas_nh3calib = FigureCanvasQTAgg(self.fig_nh3calib)
        self.navi_nh3calib = NavigationToolbar2QT(self.canvas_nh3calib, w, coordinates=False)
        self.ax_nh3calib.set_xlim(0., 100.)
        self.ax_nh3calib.set_xlabel('alpha [%]', fontsize=fs)
        self.ax_nh3calib.set_ylabel('Potential [mV]', fontsize=fs)
        self.fig_nh3calib.subplots_adjust(left=0.15, right=0.95, bottom=0.29, top=0.85)
        sns.despine()

        # connect onclick event with function
        # self.fig_nh3calib.canvas.mpl_connect('button_press_event', self.onclick_timedrive)

        nh3calib_group = QGroupBox("NH3 Sensor Calibration")
        grid_nh3calib = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_mbottom.addWidget(nh3calib_group)
        nh3calib_group.setLayout(grid_nh3calib)
        grid_nh3calib.addWidget(self.canvas_nh3calib)
        grid_nh3calib.addWidget(self.navi_nh3calib)

        # ---------------------------------------------------
        # nh3+nh4 simulation
        self.fig_nh3sim, self.ax_nh3sim = plt.subplots()
        self.canvas_nh3sim = FigureCanvasQTAgg(self.fig_nh3sim)
        self.navi_nh3sim = NavigationToolbar2QT(self.canvas_nh3sim, w, coordinates=False)
        self.ax_nh3sim.set_xlabel('Time [s]', fontsize=fs)
        self.ax_nh3sim.set_ylabel('NH$_3$ / NH$_4^+$ [ppm]', fontsize=fs)
        self.fig_nh3sim.subplots_adjust(left=0.15, right=0.95, bottom=0.29, top=0.85)
        sns.despine()

        # connect onclick event with function
        # self.fig_nh3calib.canvas.mpl_connect('button_press_event', self.onclick_timedrive)

        nh3sim_group = QGroupBox("NH3 / NH4+ Simulation")
        nh3sim_group.setMinimumWidth(200)
        grid_nh3sim = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_rtop.addWidget(nh3sim_group)
        nh3sim_group.setLayout(grid_nh3sim)
        grid_nh3sim.addWidget(self.canvas_nh3sim)
        grid_nh3sim.addWidget(self.navi_nh3sim)

        # TAN simulation
        self.fig_tanim, self.ax_tansim = plt.subplots()
        self.canvas_tansim = FigureCanvasQTAgg(self.fig_tanim)
        self.navi_tansim = NavigationToolbar2QT(self.canvas_tansim, w, coordinates=False)
        self.ax_tansim.set_xlabel('Time [s]', fontsize=fs)
        self.ax_tansim.set_ylabel('TAN [ppm]', fontsize=fs)
        self.fig_tanim.subplots_adjust(left=0.15, right=0.95, bottom=0.29, top=0.85)
        sns.despine()

        # connect onclick event with function
        # self.fig_nh3calib.canvas.mpl_connect('button_press_event', self.onclick_timedrive)

        tansim_group = QGroupBox("Total Ammonia Simulation")
        tansim_group.setMinimumWidth(200)
        grid_tansim = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_rbottom.addWidget(tansim_group)
        tansim_group.setLayout(grid_tansim)
        grid_tansim.addWidget(self.canvas_tansim)
        grid_tansim.addWidget(self.navi_tansim)

    # ---------------------------------------------------------------
    #     # right side bottom
    #     vbox_right.addLayout(hbox_bottom)
    #
    #     vbox_bottom_left = QVBoxLayout()
    #     vbox_bottom_right = QVBoxLayout()
    #     hbox_bottom.addLayout(vbox_bottom_left)
    #     hbox_bottom.addLayout(vbox_bottom_right)
    #
    #     # Plot for histogram (left side in right half of main window)
    #     self.fig_histogram, self.ax_histogram = plt.subplots()
    #     self.canvas_histo = FigureCanvasQTAgg(self.fig_histogram)
    #     self.navi_histo = NavigationToolbar2QT(self.canvas_histo, w)
    #     self.ax_histogram = plt.gca()
    #     self.ax_histogram.set_xlim(0, 8)
    #     self.ax_histogram.set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])
    #     self.led_selection = ['380 nm', '403 nm', '438 nm', '453 nm', '472 nm', '526 nm', '593 nm', '640 nm']
    #     self.ax_histogram.set_xticklabels(self.led_selection, rotation=15)
    #     self.ax_histogram.set_xlabel('Wavelength [nm]')
    #     self.ax_histogram.set_ylabel('Relative signal intensity [rfu]')
    #     self.fig_histogram.subplots_adjust(left=0.1, right=0.9, bottom=0.18, top=0.85)
    #
    #     # Plot for score plot (right side in right half of main window)
    #     self.fig_scoreplot = plt.figure()
    #     if self.threedimensional_plot_checkbox.isChecked():
    #         # 3D Plot
    #         self.ax_scoreplot = self.fig_scoreplot.gca(projection='3d')
    #         self.ax_scoreplot.view_init(elev=16., azim=57)
    #         self.canvas_score = FigureCanvasQTAgg(self.fig_scoreplot)
    #         self.navi_score = NavigationToolbar2QT(self.canvas_score, w)
    #         self.ax_scoreplot = plt.gca()
    #         self.ax_scoreplot.set_xlabel('LDA 1', fontsize=13, labelpad=10)
    #         self.ax_scoreplot.set_ylabel('LDA 2', fontsize=13, labelpad=10)
    #         self.ax_scoreplot.set_zlabel('LDA 3', fontsize=13, labelpad=10)
    #         self.fig_scoreplot.subplots_adjust(left=0.1, right=0.9, bottom=0.18, top=0.85)
    #     else:
    #         # 2D Plot
    #         self.ax_scoreplot = self.fig_scoreplot.add_subplot(111)
    #         self.canvas_score = FigureCanvasQTAgg(self.fig_scoreplot)
    #         self.navi_score = NavigationToolbar2QT(self.canvas_score, w)
    #         self.ax_scoreplot = plt.gca()
    #         self.ax_scoreplot.set_xlabel('LDA 1', fontsize=13, labelpad=10)
    #         self.ax_scoreplot.set_ylabel('LDA 2', fontsize=13, labelpad=10)
    #         self.fig_scoreplot.subplots_adjust(left=0.10, right=0.9, bottom=0.18, top=0.85)
    #
    #     # create GroupBox to structure the layout
    #     result_histogram_group = QGroupBox("Histogram")
    #     grid_histogram = QGridLayout()
    #     result_scoreplot_group = QGroupBox("Score plot")
    #     grid_scoreplot = QGridLayout()
    #
    #     # add GroupBox to layout and load buttons in GroupBox
    #     vbox_bottom_left.addWidget(result_histogram_group)
    #     result_histogram_group.setLayout(grid_histogram)
    #     vbox_bottom_right.addWidget(result_scoreplot_group)
    #     result_scoreplot_group.setLayout(grid_scoreplot)
    #
    #     grid_histogram.addWidget(self.canvas_histo)
    #     grid_histogram.addWidget(self.navi_histo)
    #     grid_scoreplot.addWidget(self.canvas_score)
    #     grid_scoreplot.addWidget(self.navi_score)
    #
    #     #self.canvas_histo.setMinimumWidth(500)
    #     #self.canvas_score.setMinimumWidth(800)

        self.show()

###################################################################################################################
# Functions for analysis
###################################################################################################################
    # print parameter
    def print_temperature(self):
        print('Temperature: ', self.temperature_edit.text(), 'degC')

    def print_phrange(self):
        print('pH range: ', self.phrange_edit.text())

    def print_samplingrate(self):
        print('Sampling rate: ', self.smpgrate_edit.text(), 's')

    def print_ph_t90(self):
        print('pH sensor response: ', self.ph_t90_edit.text(), 's')

    def print_ph_signal(self):
        print('pH sensor signal (min, max): ', self.ph_signal_edit.text(), 'mV')

    def print_ph_resolution(self):
        print('pH sensor resolution: ', self.ph_res_edit.text(), 'mV')

    def print_ph_reference_pot(self):
        print('Potential reference electrode: ', self.ph_ref_edit.text(), 'mV')

    def print_nh3_t90(self):
        print('NH3 sensor response: ', self.nh3_t90_edit.text(), 's')

    def print_nh3_signal(self):
        print('nH3 sensor signal: ', self.nh3_signal_edit.text(), 'mV')

    def print_nh3_resolution(self):
        print('Sensor resolution: ', self.nh3_res_edit.text(), 'mV')

    def print_nh3_pka(self):
        print('pKa: ', self.nh3_pka_edit.text())

    def print_nh3_concentration(self):
        print('NH3 concentration at pKa: ', self.nh3_cGG_edit.text(), 'ppm')

    def print_nh3_alpha(self):
        print('NH3 proportion range: ', self.nh3_alpha_edit.text(), '%')



    # def open_sample(self):
    #     self.fname = QFileDialog.getOpenFileName(self, "Select a measurement file", 'measurement/')[0]
    #
    #     if not self.fname:
    #         return
    #     self.read_sample_name(self.fname)
    #
    # def clear_all(self):
    #     self.sample_edit.clear()
    #     self.blank_edit.clear()
    #     self.prescan_edit.clear()
    #     self.database_edit.clear()
    #     self.ex_correction_edit.clear()
    #     self.em_correction_edit.clear()
    #     if self.xcoords_prescan:
    #         self.xcoords_prescan = None
    #     if self.xcoords:
    #         self.xcoords.clear()
    #
    # def reportInput_blank_correction(self):
    #     print('blank correction: ', self.blank_cor_checkbox.isChecked())
    #     return self.blank_cor_checkbox.isChecked()
    #
    # def reportInput_external_blank(self):
    #     print('external blank: ', self.blank_externalfile_checkbox.isChecked())
    #     return self.blank_externalfile_checkbox.isChecked()
    #
    # def reportInput_correction(self):
    #     print('em-/ex-correction: ', self.correction_checkbox.isChecked())
    #
    # def reportInput_led_calibration_fit(self):
    #     print('led-calibration fit linear:', self.calibration_linearfit_checkbox.isChecked())
    #     return self.calibration_linearfit_checkbox.isChecked()
    #
    # def reportInput_single_peak(self):
    #     print('single peak evaluation: ', self.singlepeak_detect_checkbox.isChecked())
    #     return self.singlepeak_detect_checkbox.isChecked()
    #
    # def reportInput_priority(self):
    #     print('Priority for dinophyta: ', self.priority_checkbox.isChecked())
    #     return self.priority_checkbox.isChecked()
    #
    # def reportSampling_evaluation(self):
    #     print('Sample will be evaluated and plotted: ', self.sampling_checkbox.isChecked())
    #     return self.sampling_checkbox.isChecked()
    #
    # def reportInput_normalize(self):
    #     print('data normalization: ', self.normalize_checkbox.isChecked())
    #     return self.normalize_checkbox.isChecked()
    #
    # def reportInput_standardize(self):
    #     print('data standardization: ', self.standardize_checkbox.isChecked())
    #     return self.standardize_checkbox.isChecked()
    #
    # def reportInput_threedimensional(self):
    #     print('3D plot for score plot: ', self.threedimensional_plot_checkbox.isChecked())
    #     if self.threedimensional_plot_checkbox.isChecked() is True:
    #         self.twodimensional_plot_checkbox.setCheckState(False)
    #     return self.threedimensional_plot_checkbox.isChecked()
    #
    # def reportInput_twodimensional(self):
    #     print('2D plot for score plot: ', self.twodimensional_plot_checkbox.isChecked())
    #     if self.twodimensional_plot_checkbox.isChecked() is True:
    #         self.threedimensional_plot_checkbox.setCheckState(False)
    #     return self.twodimensional_plot_checkbox.isChecked()
    #
    # def reportInput_LoD(self):
    #     print('peak detection with LoD ', self.peak_detect_checkbox.isChecked())
    #     return self.peak_detect_checkbox.isChecked()
    #
    # def print_pumprate(self):
    #     print('Pumprate: ', self.pumprate_edit.text(), 'mL/min')
    #     return self.pumprate_edit.text()
    #
    # def print_separation(self):
    #     if self.separation_edit.text() == 'phylum' or self.separation_edit.text() == 'order':
    #         print('Separation level: ', self.separation_edit.text())
    #         return self.separation_edit.text()
    #     elif self.separation_edit.text() == 'family' or self.separation_edit.text() == 'order':
    #         print('Separation level: ', self.separation_edit.text())
    #         return self.separation_edit.text()
    #     else:
    #         separation_level_failed = QMessageBox()
    #         separation_level_failed.setIcon(QMessageBox.Information)
    #         separation_level_failed.setText("Invalid separation level!")
    #         separation_level_failed.setInformativeText("Choose either 'phylum', 'class', 'order' or 'family' ... ")
    #         separation_level_failed.setWindowTitle("Error!")
    #         separation_level_failed.exec_()
    #         return
    #
    # def print_limit(self):
    #     print('Limit for score plot: ', self.limit_edit.text())
    #     return self.limit_edit.text()
    #
    # def onclick_timedrive(self, event):
    #     modifiers = QApplication.keyboardModifiers()
    #     if modifiers != Qt.ControlModifier:  # change selected range
    #         return
    #     if len(self.xcoords) >= 4:
    #         return
    #     if event.xdata == None:
    #         if len(self.xcoords) < 2:
    #             event.xdata = self.loaded_data['l_corr'].index.max()
    #         else:
    #             event.xdata = 0
    #     self.xcoords.append(event.xdata)
    #     self.ax_timedrive.vlines(x=self.xcoords, ymin=self.loaded_data['l_corr'].min().min(),
    #                              ymax=self.loaded_data['l_corr'].max().max(), lw=0.5)
    #     if len(self.xcoords) == 2:
    #         self.ax_timedrive.axvspan(self.xcoords[0], self.xcoords[1], color='grey', alpha=0.3)
    #     elif len(self.xcoords) == 4:
    #         self.ax_timedrive.axvspan(self.xcoords[0], self.xcoords[1], color='grey', alpha=0.3)
    #         self.ax_timedrive.axvspan(self.xcoords[2], self.xcoords[3], color='grey', alpha=0.3)
    #     self.fig_timedrive.canvas.draw()
    #
    # def plot_timedrive(self, df, name, date, ax, f, unit):
    #     color_LED = []
    #     for i in df.columns:
    #         color_LED.append(led_color_dict[i])
    #
    #     # plotting spectra
    #     ylim_max = pd.DataFrame(np.zeros(shape=(len(df.columns), 1)), index=df.columns).T
    #     ylim_min = pd.DataFrame(np.zeros(shape=(len(df.columns), 1)), index=df.columns).T
    #
    #     for c, p in zip(color_LED, df):
    #         if df[p].dropna().empty is True:
    #             df[p][np.isnan(df[p])] = 0
    #             df[p].plot(ax=ax, color=c, linewidth=1.75)
    #         else:
    #             df[p].dropna().plot(ax=ax, color=c, label=p, linewidth=1.75)
    #             ylim_max[p] = df[p].dropna().max()
    #             ylim_min[p] = df[p].dropna().min()
    #
    #     # General layout-stuff
    #     ax.set_xlabel('Time [s]')
    #     ax.set_ylabel('Rel. Fluorescence intensity [{}]'.format(unit))
    #     ax.legend(loc=0, ncol=1, frameon=True, fancybox=True, framealpha=0.5,
    #               fontsize=9) # 'upper center', bbox_to_anchor=(1.08, 1.)
    #
    #     # Define plotting area. Default is 2 but if max value is higher it has to be rearranged
    #     x_max = df.index[-1] * 1.05
    #     x_min = np.abs(df.index[0]) - np.abs(df.index[-1])*0.05
    #     ax.set_xlim(x_min, x_max)
    #
    #     y_max = ylim_max.max(axis=1).values[0] * 1.05
    #     y_min = ylim_min.min(axis=1).values[0] * 1.05
    #     ax.set_ylim(y_min, y_max)
    #     ax.set_title("{} {}/{}/{} {}:{}h - \n Select time range for sample and "
    #                  "blank.".format(name, date[6:8], date[4:6], date[:4], date[8:10], date[10:12]),
    #                  fontsize=11)
    #
    #     f.tight_layout()
    #     f.canvas.draw()
    #
    # def plot_histogram(self, mean, ax, f):
    #     # plot histogram: Relative intensity @ different excitation LEDs.
    #     mean = mean.sort_index()
    #
    #     # prepare general information about the sample and the setup
    #     LED_color = []
    #     LED_wl = []
    #     for i in mean.index:
    #         if type(i) == str:
    #             if len(i.split(' ')) == 1:
    #                 # led as string without 'nm'
    #                 j = i + ' nm'
    #             else:
    #                 # led with 'nm'
    #                 j = i
    #         else:
    #             j = str(i) + ' nm'
    #         LED_wl.append(j)
    #         LED_color.append(led_color_dict[j])
    #
    #     mean.index = LED_wl
    #
    #     # normalize mean
    #     means = mean / mean.max()
    #     for i in means.index:
    #         if (means.ix[i, :].values[0]) < 0:
    #             means.ix[i, :].values[0] = 0
    #
    #     self.led_total = pd.DataFrame(self.led_selection)
    #     x = []
    #     for i in mean.index:
    #         x.append(self.led_total[self.led_total[0] == i].index[0])
    #
    #     for k, l in enumerate(mean.index):
    #         ax.bar(x[k]+0.5, means.ix[l, :], width=0.9, color=LED_color[k])
    #     ax = plt.gca()
    #     ax.set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])
    #     ax.set_xticklabels(mean.index)
    #     ax.xaxis.grid(False)
    #     plt.setp(ax.get_xticklabels(), rotation=45, ha='center')
    #
    #     ax.set_xlabel('Wavelength [nm]')
    #     ax.set_ylabel('Relative signal intensity [rfu]')
    #     f.tight_layout()
    #
    #     f.canvas.draw()
    #
    # def plot_distribution_2d(self, f, ax, d, df_score, color, alg_group, alg_phylum, phyl_group, separation):
    #     ax.cla()
    #     # preparation of figure plot
    #     if ax is None:
    #         f, ax = plt.subplots()
    #
    #     # plotting the centers of each algal class
    #     for el in d.index:
    #         ax.scatter(d.ix[el, 0], d.ix[el, 1], facecolor=color[el], edgecolor='k', s=60)
    #
    #     # calculate the standard deviation within one class to built up a solid (sphere with 3 different radii)
    #     # around the centroid using spherical coordinates
    #     for i in d.index:
    #         # replace sample name (phyl_group.index) with phylum name (in phyl_group['phylum_label'])
    #         phyl = alg_group[alg_group[separation] == i]['phylum'].values[0]
    #         if phyl in phyl_group['phylum_label'].values:
    #             pass
    #         else:
    #             phyl_group.ix[i, 'phylum_label'] = phyl
    #             phyl_group.ix[i, 'color'] = alg_phylum.ix[phyl, :].values
    #
    #         rx = np.sqrt(d['LDA1var'].ix[i])
    #         ry = np.sqrt(d['LDA2var'].ix[i])
    #         c_x = d.ix[i]['LDA1']
    #         c_y = d.ix[i]['LDA2']
    #         ells = Ellipse(xy=[c_x, c_y], width=rx, height=ry, angle=0, edgecolor=color[i], lw=1, facecolor=color[i],
    #                        alpha=0.6, label=i)
    #         ax.add_artist(ells)
    #
    #         ells2 = Ellipse(xy=[c_x, c_y], width=2*rx, height=2*ry, angle=0, edgecolor=color[i], lw=1,
    #                         facecolor=color[i], alpha=0.4)
    #         ax.add_artist(ells2)
    #
    #         ells3 = Ellipse(xy=[c_x, c_y], width=3*rx, height=3*ry, angle=0, edgecolor=color[i], lw=0.5,
    #                         facecolor=color[i], alpha=0.1)
    #         ax.add_artist(ells3)
    #
    #     # patch = pd.DataFrame(np.zeros(shape=(len(phyl_group), 2)), index=phyl_group.index)
    #     patch = []
    #     for i in phyl_group.index:
    #         patch.append(mpatches.Patch(color=phyl_group.ix[i, 'color'][0], label=phyl_group.ix[i, 'phylum_label']))
    #         ax.legend(handles=patch, loc="upper center", bbox_to_anchor=(1.2, 0.9), frameon=True, fontsize=11)
    #
    #     plt.setp(ax.get_xticklabels(), fontsize=13)
    #     plt.setp(ax.get_yticklabels(), fontsize=13)
    #     ax.set_xlabel('LDA1', fontsize=13, labelpad=5)
    #     ax.set_ylabel('LDA2', fontsize=13, labelpad=5)
    #     plt.title('')
    #
    #     # plotting the sample scores
    #     if self.sampling_checkbox.isChecked():
    #         for i in range(len(df_score.T)):
    #             ax.plot(df_score.ix[0, 0], df_score.ix[1, 0], marker='^', markersize=14, color='orangered', label='')
    #     f.subplots_adjust(left=0.1, right=0.75, bottom=0.18, top=0.85)
    #
    #     f.canvas.draw()
    #
    # def plot_distribution_3d(self, f, ax, d, df_score, color, alg_group, alg_phylum, phyl_group, separation):
    #     ax.cla()
    #     # preparation of figure plot
    #     if ax is None:
    #         f = plt.figure()
    #         ax = f.gca(projection='3d')
    #         ax.set_aspect('auto')
    #     # initial view of score plot to enhance the separation of algae and cyanos
    #     ax.view_init(elev=19., azim=-67)
    #
    #     # plotting the centers of each algal class
    #     for el in d.index:
    #         ax.scatter(d.ix[el, 0], d.ix[el, 1], d.ix[el, 2], marker='.', color='k', s=60)
    #
    #     # calculate the standard deviation within one class to built up a solid (sphere with 3 different radii)
    #     # around the centroid using spherical coordinates
    #     for i in d.index:
    #         # replace sample name (phyl_group.index) with phylum name (in phyl_group['phylum_label'])
    #         phyl = alg_group[alg_group[separation] == i]['phylum'].values[0]
    #         if phyl in phyl_group['phylum_label'].values:
    #             pass
    #         else:
    #             phyl_group.ix[i, 'phylum_label'] = phyl
    #
    #             phyl_group.ix[i, 'color'] = alg_phylum.ix[phyl, :].values
    #
    #         rx = np.sqrt(d['LDA1var'].ix[i])
    #         ry = np.sqrt(d['LDA2var'].ix[i])
    #         rz = np.sqrt(d['LDA3var'].ix[i])
    #         c_x = d.ix[i]['LDA1']
    #         c_y = d.ix[i]['LDA2']
    #         c_z = d.ix[i]['LDA3']
    #
    #         u, v = np.mgrid[0:2 * np.pi:10j, 0:np.pi:20j]
    #         x = rx * np.cos(u) * np.sin(v) + c_x
    #         y = ry * np.sin(u) * np.sin(v) + c_y
    #         z = rz * np.cos(v) + c_z
    #         ax.plot_wireframe(x, y, z, color=color[i], alpha=0.5, linewidth=1, label=phyl)
    #
    #         x1 = 2 * rx * np.cos(u) * np.sin(v) + c_x
    #         y1 = 2 * ry * np.sin(u) * np.sin(v) + c_y
    #         z1 = 2 * rz * np.cos(v) + c_z
    #         ax.plot_wireframe(x1, y1, z1, color=color[i], alpha=0.2, linewidth=1)
    #
    #         x2 = 3 * rx * np.cos(u) * np.sin(v) + c_x
    #         y2 = 3 * ry * np.sin(u) * np.sin(v) + c_y
    #         z2 = 3 * rz * np.cos(v) + c_z
    #         ax.plot_wireframe(x2, y2, z2, color=color[i], alpha=0.15, linewidth=0.5)
    #
    #     patch = []
    #     for i in phyl_group.index:
    #         patch.append(mpatches.Patch(color=phyl_group.ix[i, 'color'][0], label=phyl_group.ix[i, 'phylum_label']))
    #         ax.legend(handles=patch, loc="upper center", bbox_to_anchor=(0., 0.9), frameon=True, fancybox=True,
    #                   fontsize=10)
    #
    #     plt.setp(ax.get_xticklabels(), va='center', ha='left', fontsize=13)
    #     plt.setp(ax.get_yticklabels(), va='center', ha='left', fontsize=13)
    #     plt.setp(ax.get_zticklabels(), va='center', fontsize=13)
    #
    #     ax.set_xlabel('LDA1', fontsize=14, labelpad=16, rotation=-2)
    #     ax.set_ylabel('LDA2', fontsize=14, labelpad=14, rotation=18)
    #     ax.set_zlabel('LDA3', fontsize=14, labelpad=10, rotation=90)
    #     plt.title(' ')
    #
    #     # plotting the sample scores
    #     if self.sampling_checkbox.isChecked():
    #         for i in range(len(df_score.T)):
    #             ax.scatter(df_score.ix[0, 0], df_score.ix[1, 0], df_score.ix[2, 0], marker='^', s=300,
    #                        color='orangered', label='')
    #
    #     f.subplots_adjust(left=0.1, right=0.9, bottom=0.06, top=0.99)
    #
    #     f.canvas.draw()

#################################################################################################
#   Load data and select range for sample and blank
#################################################################################################
    def load_timedrive(self):
        self.message.clear()
        self.message.setText(' ')
        self.report.clear()
        self.report_sup.clear()
        self.ax_timedrive.cla()
        self.ax_timedrive.set_xlim(0, 10)
        self.ax_timedrive.set_xlabel('Time [s]')
        self.ax_timedrive.set_ylabel('Rel. intensity [pW]')
        self.fig_timedrive.canvas.draw()
        self.ax_histogram.cla()
        self.ax_histogram.set_xlim(0, 8)
        self.ax_histogram.set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])
        self.ax_histogram.set_xticklabels(self.led_selection, rotation=15)
        self.ax_histogram.set_xlabel('Wavelength [nm]')
        self.ax_histogram.set_ylabel('Relative signal intensity [rfu]')
        self.fig_histogram.canvas.draw()

        if self.threedimensional_plot_checkbox.isChecked():
            # 3D Plot redrawing
            self.ax_scoreplot.cla()
            self.fig_scoreplot.clear()
            self.ax_scoreplot = self.fig_scoreplot.gca(projection='3d')
            # self.ax_scoreplot.set_aspect('normal')
            self.ax_scoreplot.set_xlabel('LDA 1', fontsize=13, labelpad=10)
            self.ax_scoreplot.set_ylabel('LDA 2', fontsize=13, labelpad=10)
            self.ax_scoreplot.set_zlabel('LDA 3', fontsize=13, labelpad=10)
            self.fig_scoreplot.subplots_adjust(left=0.1, right=0.9, bottom=0.18, top=0.85)
            self.fig_scoreplot.canvas.draw()
        else:
            # 2D Plot activated
            self.ax_scoreplot.cla()
            self.fig_scoreplot.clear()
            self.ax_scoreplot = self.fig_scoreplot.add_subplot(111)
            self.ax_scoreplot.set_xlabel('LDA 1', fontsize=13, labelpad=10)
            self.ax_scoreplot.set_ylabel('LDA 2', fontsize=13, labelpad=10)
            self.fig_scoreplot.subplots_adjust(left=0.12, right=0.9, bottom=0.18, top=0.85)
            self.fig_scoreplot.canvas.draw()

        # select if pre-scanned data or raw data are analysed
        if self.sample_edit.toPlainText():
            # raw data analysis
            self.prescan_edit.clear()
            self.fname_prescan = None
            self.xcoords_prescan = None

        elif self.prescan_edit.toPlainText():
            # already pre-scanned data are analysed
            self.sample_edit.clear()
            self.fname_sample = None
            self.xcoords.clear()
            self.correction_checkbox.setCheckState(False)
            self.calibration_linearfit_checkbox.setCheckState(False)

        # find sample filename
        try:
            self.fname
        except:
            try:
                # if pre-scanned data are analysed
                self.fname_prescan
            except:
                run_sample_load_failed = QMessageBox()
                run_sample_load_failed.setIcon(QMessageBox.Information)
                run_sample_load_failed.setText("Sample file is missing!")
                run_sample_load_failed.setInformativeText("Choose a specific sample file for analysis or choose an "
                                                          "already analysed prescan file.")
                run_sample_load_failed.setWindowTitle("Error!")
                run_sample_load_failed.exec_()
                return

        # if data should be corrected
        if self.correction_checkbox.isChecked() is True:
            correction = True
            # if correction is chosen, you need the em- and ex- correction files!
            # device for emission correction
            try:
                self.number_em  # number 0, 1, 2, 3
            except:
                run_correction_em_failed = QMessageBox()
                run_correction_em_failed.setIcon(QMessageBox.Information)
                run_correction_em_failed.setText("File for emission correction is missing!")
                run_correction_em_failed.setInformativeText("Choose a correction file for analysis.")
                run_correction_em_failed.setWindowTitle("Error!")
                run_correction_em_failed.exec_()
                return

            # excitation correction
            # linear fit between 30-50mA or whole current between 9-97mA
            if self.calibration_linearfit_checkbox.isChecked() is True:
                full_calibration = False
            else:
                full_calibration = True

            if self.ex_correction_edit.toPlainText():
                # manual choice of correction file!
                # check if the correct file was chosen (spectrum; not reference)
                if self.fname_ex.split('/')[-1].split('_')[4] == 'reference':
                    self.message.append('Wrong file for LED correction was chosen! '
                                        'Using spectrum correction instead of reference correction ...')
                    fname_ex_corr = self.fname_ex.split('reference')[0] + 'spectrum' + \
                                    self.fname_ex.split('reference')[1]
                    kappa_spec = fname_ex_corr
                else:
                    kappa_spec = self.fname_ex

            else:
                # automatic choice of correction file
                self.message.append('No specific file for excitation correction is chosen. Choose automatically ...')
                kappa_spec = None
        else:
            # No excitation correction
            correction = False
            device = None
            kappa_spec = None
            self.number_em = 'device-1'
            full_calibration = False

        # pump rate either from file or from GUI
        if not self.pumprate_edit.text():
            self.message.append("No external pump rate defined. Take information from header!")
            pumprate = None
        else:
            pumprate = float(self.pumprate_edit.text())

        # no additional offset compensation when database is updated
        additional = False

        # single event evaluation
        # !! TODO: Include single peak evaluation!
        if self.singlepeak_detect_checkbox.isChecked() is True:
            single_events = True
        else:
            single_events = False

        # peak_detection means LoD is used
        if self.peak_detect_checkbox.isChecked() is True:
            peak_detection = True
        else:
            peak_detection = False

        # other parameter not used in the normal analysis
        ampli = None
        factor = 1

        # blank_corr of data
        if self.blank_cor_checkbox.isChecked() is True:
            blank_correction = True
        else:
            blank_correction = False

        # use external file instead of blank stored in the header
        if self.blank_externalfile_checkbox.isChecked() is True:
            try:
                self.fname_blank
            except:
                run_blank_external_failed = QMessageBox()
                run_blank_external_failed.setIcon(QMessageBox.Information)
                run_blank_external_failed.setText("File for external blank correction is missing!")
                run_blank_external_failed.setInformativeText("Choose an external blank file for analysis.")
                run_blank_external_failed.setWindowTitle("Error!")
                run_blank_external_failed.exec_()
                return
            # blank is not (ex- or em-)corrected so far
            [blank, header_bl, unit_bl] = alg.read_rawdata(filename=self.fname_blank, additional=additional, co=None,
                                                           factor=factor, blank_corr=blank_correction,
                                                           blank_mean_ex=None, blank_std_ex=None, plot_raw=False)
            # convert the blank in nW
            if unit_bl == 'nW':
                blank_mean_ex =blank.mean().tolist()
                blank_std_ex = blank.std().tolist()
            elif unit_bl == 'pW':
                blank_mean_ex =blank.mean().tolist() / 1000
                blank_std_ex = blank.std().tolist() / 1000
                unit_bl = 'nW'
            elif unit_bl == 'µW':
                blank_mean_ex =blank.mean().tolist() * 1000
                blank_std_ex = blank.std().tolist() * 1000
                unit_bl = 'nW'
            else:
                self.message.append('The blank is unusually high! Please check the blank ....')
                return
        else:
            blank_mean_ex = None
            blank_std_ex = None

        # selection of xcoords depending on sample type
        # -----------------------------------------------------------------------------------------------------------
        # RAW DATA ANALYSIS
        # -----------------------------------------------------------------------------------------------------------
        if self.sample_edit.toPlainText() and not self.prescan_edit.toPlainText():
            # raw data do not have a x-coord selection. the xcoord prescan should be empty
            if self.xcoords_prescan:
                self.xcoords_prescan = None
            if self.xcoords:
                self.xcoords.clear()

            print(self.name, self.number_em.split('-')[1], kappa_spec, pumprate, ampli)
            print(correction, full_calibration, blank_correction, factor)

            # Load data and correct them with prescan_load_file function.
            [l, l_corr, header, firstline, current, date, self.sample_name, blank_mean, blank_std, blank_corrected,
             rg9_sample, rg665_sample, volume, pumprate, unit, unit_corr, unit_bl,
             path] = algae_analysis.prescan_load_file(filename=self.fname, device=self.number_em.split('-')[1],
                                                      kappa_spec=kappa_spec, pumprate=pumprate, ampli=ampli,
                                                      correction=correction, full_calibration=full_calibration,
                                                      blank_corr=blank_correction, factor=factor,
                                                      blank_mean_ex=blank_mean_ex,
                                                      blank_std_ex=blank_std_ex, additional=additional)

            self.led_total = l.columns
            # Dataframe reduction for analysis according to selected LEDs
            self.led_used = algae_analysis.led_reduction(LED380_checkbox=self.LED380_checkbox.isChecked(),
                                                         LED403_checkbox=self.LED403_checkbox.isChecked(),
                                                         LED438_checkbox=self.LED438_checkbox.isChecked(),
                                                         LED453_checkbox=self.LED453_checkbox.isChecked(),
                                                         LED472_checkbox=self.LED472_checkbox.isChecked(),
                                                         LED526_checkbox=self.LED526_checkbox.isChecked(),
                                                         LED593_checkbox=self.LED593_checkbox.isChecked(),
                                                         LED640_checkbox=self.LED640_checkbox.isChecked())

            l_red = pd.DataFrame(np.zeros(shape=(len(l.index), 0)), index=l.index)
            l_corr_red = pd.DataFrame(np.zeros(shape=(len(l_corr.index), 0)), index=l_corr.index)
            for i in l.columns:
                if self.led_used.ix[0, i] == True:
                    l_red.ix[:, i] = l.ix[:, i]
            for i in l_corr.columns:
                if self.led_used.ix[0, i] == True:
                    l_corr_red.ix[:, i] = l_corr.ix[:, i]
            # sample datas are sorted
            l_corr_red = l_corr_red.sort_index(axis=1)

            # Store parameter for further evaluation
            self.loaded_data = {'l': l_red, 'l_corr': l_corr_red, 'header': header, 'firstline': firstline,
                                'current': current, 'date': date, 'name': self.sample_name, 'blank_mean': blank_mean,
                                'blank_std': blank_std, 'rg9_sample': rg9_sample, 'rg665_sample': rg665_sample,
                                'volume': volume, 'path': path, 'full_calibration': full_calibration,
                                'kappa_spec': kappa_spec, 'correction': correction,
                                'device': self.number_em.split('-')[1], 'blank_corr': blank_correction,
                                'peak_detection': peak_detection, 'additional': additional, 'unit': unit_corr,
                                'pumprate': pumprate, 'unit_blank': unit_bl}

            # Plotting corrected data to select the time-range for sample and blank. Time-range stored in xcoords.
            self.plot_timedrive(df=self.loaded_data['l_corr'], name=self.loaded_data['name'],
                                date=self.loaded_data['date'],
                                f=self.fig_timedrive, ax=self.ax_timedrive, unit=unit)

        # -----------------------------------------------------------------------------------------------------------
        # PRESCANNED DATA ANALYSIS
        # -----------------------------------------------------------------------------------------------------------
        elif self.prescan_edit.toPlainText() and not self.sample_edit.toPlainText():
            # Use histogram
            self.message.append('Already pre-scanned data chosen. Apply LDA...')

            [led_mean, self.sample_name, self.cur, self.ampli, volume, LoD, counted_cells, unit,
             date] = algae_analysis.processed_data_load_file(self.fname_prescan)

            # Dataframe reduction for analysis according to selected LEDs
            self.led_used = algae_analysis.led_reduction(LED380_checkbox=self.LED380_checkbox.isChecked(),
                                                         LED403_checkbox=self.LED403_checkbox.isChecked(),
                                                         LED438_checkbox=self.LED438_checkbox.isChecked(),
                                                         LED453_checkbox=self.LED453_checkbox.isChecked(),
                                                         LED472_checkbox=self.LED472_checkbox.isChecked(),
                                                         LED526_checkbox=self.LED526_checkbox.isChecked(),
                                                         LED593_checkbox=self.LED593_checkbox.isChecked(),
                                                         LED640_checkbox=self.LED640_checkbox.isChecked())

            # sample reduction
            self.led_mean_red = pd.DataFrame(np.zeros(shape=(0, len(led_mean.columns))),
                                             columns=led_mean.columns)
            for i in self.led_used.columns:
                if self.led_used.ix[0, i] == True:
                    self.led_mean_red.ix[i, self.led_mean_red.columns[0]] =\
                        led_mean.ix[i[:3], led_mean.columns[0]]
                else:
                    pass

            sample_name = self.fname_prescan.split('.')[0].split('/')[-1]
            path = self.fname_prescan.split('.')[0][:-len(sample_name)]
            self.plot_histogram(mean=self.led_mean_red, f=self.fig_histogram, ax=self.ax_histogram)

            self.loaded_data = {'full_calibration': full_calibration, 'kappa_spec': kappa_spec, 'volume': volume,
                                'correction': correction, 'device': self.number_em.split('-')[1], 'date': date,
                                'blank_corr': blank_correction, 'peak_detection': peak_detection, 'LoD': LoD,
                                'unit': unit, 'additional': additional, 'pumprate': pumprate, 'path': path,
                                'counted_cells': counted_cells}

        elif self.prescan_edit.toPlainText() and self.sample_edit.toPlainText():
            run_both_sample_failed = QMessageBox()
            run_both_sample_failed.setIcon(QMessageBox.Information)
            run_both_sample_failed.setText("Select just one sample file!")
            run_both_sample_failed.setInformativeText("Either raw data with time drive or already processed sample "
                                                      "file")
            run_both_sample_failed.setWindowTitle("Too many input!")
            run_both_sample_failed.exec_()
        else:
            run_missing_sample_failed = QMessageBox()
            run_missing_sample_failed.setIcon(QMessageBox.Information)
            run_missing_sample_failed.setText("Choose sample file!")
            run_missing_sample_failed.setInformativeText("Either raw data or already processed sample file")
            run_missing_sample_failed.setWindowTitle("Missing input!")
            run_missing_sample_failed.exec_()
            return

        self.run_button.setEnabled(True)


#################################################################################################
#   Run analysis
#################################################################################################
    def run_analysis(self):
        self.message.clear()
        self.report.clear()
        self.report_sup.clear()
        if self.xcoords:
            # xcoords from time drive selected -> currently no histogram
            self.ax_histogram.clear()
            self.ax_histogram.set_xlim(0, 8)
            self.ax_histogram.set_xticks([0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5])
            self.ax_histogram.set_xticklabels(self.led_selection, rotation=15)
            self.ax_histogram.set_xlabel('Wavelength [nm]')
            self.ax_histogram.set_ylabel('Relative signal intensity [rfu]')
        if self.threedimensional_plot_checkbox.isChecked():
            # 3D Plot redrawing
            self.ax_scoreplot.cla()
            self.fig_scoreplot.clear()
            self.ax_scoreplot = self.fig_scoreplot.gca(projection='3d')
            self.ax_scoreplot.set_xlabel('LDA 1', fontsize=13, labelpad=10)
            self.ax_scoreplot.set_ylabel('LDA 2', fontsize=13, labelpad=10)
            self.ax_scoreplot.set_zlabel('LDA 3', fontsize=13, labelpad=10)
            self.fig_scoreplot.subplots_adjust(left=0.1, right=0.9, bottom=0.18, top=0.85)
            self.fig_scoreplot.canvas.draw()
        else:
            # 2D Plot activated
            self.ax_scoreplot.cla()
            self.fig_scoreplot.clear()
            self.ax_scoreplot = self.fig_scoreplot.add_subplot(111)
            self.ax_scoreplot.set_xlabel('LDA 1', fontsize=13, labelpad=10)
            self.ax_scoreplot.set_ylabel('LDA 2', fontsize=13, labelpad=10)
            self.fig_scoreplot.subplots_adjust(left=0.12, right=0.9, bottom=0.18, top=0.85)
            self.fig_scoreplot.canvas.draw()

        if not self.xcoords and not self.xcoords_prescan:
            self.message.append('No time range selected - use whole timeline as sample')
            self.xcoords = [self.loaded_data['l_corr'].index[0], self.loaded_data['l_corr'].index[-1], 0, 0]
        if len(self.xcoords) == 2:
            self.message.append('Only the sample range was selected - use blank which is stored in header')
            self.xcoords = [self.xcoords[0], self.xcoords[1], 0, 0]

        if self.priority_checkbox.isChecked() is True:
            priority = True
        else:
            priority = False

        if self.prescan_edit.toPlainText() and not self.sample_edit.toPlainText():
            # already evaluated data
            self.mean_corr = self.led_mean_red
        elif not self.prescan_edit.toPlainText() and self.sample_edit.toPlainText():
            # raw data evaluation
            self.mean_corr = self.loaded_data['l_corr'].mean()

        try:
            self.fname_database
        except:
            run_missing_database = QMessageBox()
            run_missing_database.setIcon(QMessageBox.Information)
            run_missing_database.setText("Database is missing!")
            run_missing_database.setInformativeText("Select a folder for database!")
            run_missing_database.setWindowTitle("Missing input!")
            run_missing_database.exec_()
            return

        trainingsdata = self.fname_database
        device_training = 1

        # Load trainings data, which is already corrected
        self.message.append('Start loading reference database... ')

        # reduce trainings database to selected LEDs
        [self.training_corr_red_sort, training,
         self.training_red] = algae_analysis.training_database(trainings_path=trainingsdata, led_used=self.led_used)

        # message that trainings matrix is loaded
        self.message.append('Reference data base is loaded: check ✓')

        # Load additional information, e.g. color classes and genus names
        if self.threedimensional_plot_checkbox.isChecked() is True:
            self.score_type = 3
        elif self.twodimensional_plot_checkbox.isChecked() is True:
            self.score_type = 2
        else:
            self.message.append('Error! Choose a distinct plot dimension')
            return


# =====================================================================================================================
if __name__ == '__main__':
    app = QApplication(sys.argv)
    view = MainWindow()
    view.setGeometry(50, 70, 1300, 750)
    sys.exit(app.exec_())
