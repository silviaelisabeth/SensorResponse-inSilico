__author__ = 'szieger'
__project__ = 'in silico study for sensor response'

import basics_sensorSignal as bs
import sys
from PyQt5.QtGui import *
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QMainWindow, QPushButton, QAction, qApp,
                             QGridLayout, QLabel, QLineEdit, QGroupBox, QFileDialog, QFrame, QMessageBox)
from PyQt5.QtGui import QIcon, QDoubleValidator
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT, FigureCanvasQTAgg
import numpy as np
import seaborn as sns
import os

# .....................................................................................................................
# global parameter
dcolor = dict({'pH': '#1CC49A', 'sig pH': '#4B5258', 'NH3': '#196E94', 'sig NH3': '#314945', 'NH4': '#DCA744',
               'sig NH4': '#89621A', 'TAN': '#A86349'})
ls = dict({'target': '-.', 'simulation': ':'})
save_type = ['png', 'svg']
sns.set_context('paper', font_scale=1.)


# .....................................................................................................................
class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.initUI()
        self.dic_sens_record = None
        self.setWindowIcon(QIcon('icon.png'))

    def initUI(self):
        # creating main window (GUI)
        w = QWidget()
        self.setCentralWidget(w)
        self.setWindowTitle('Sensor response')

        # ---------------------------------------------------------------------------------------
        # Menu bar - Load data, Save report, Save all, Exit
        exitAction = QAction('&Exit', self)
        exitAction.setShortcut('Ctrl+Q')
        exitAction.setStatusTip('Exit application')
        exitAction.triggered.connect(qApp.quit)
        self.statusBar()

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(exitAction)

        # ---------------------------------------------------------------------------------------
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

        # ----------------------------------------------------
        # left side of main window (-> data treatment)
        vbox_top.addWidget(w)
        vbox_top.setContentsMargins(5, 10, 10, 5)

        # ---------------------------------------------------------------------------------------------------------
        # PARAMETERS
        # general settings
        temperature_label, temperature_unit_label = QLabel(self), QLabel(self)
        temperature_label.setText('Temperature')
        temperature_unit_label.setText('degC')
        self.temperature_edit = QLineEdit(self)
        self.temperature_edit.setValidator(QDoubleValidator())
        self.temperature_edit.setAlignment(Qt.AlignRight)
        self.temperature_edit.setText('25.')

        pH_label = QLabel(self)
        pH_label.setText('pH range')
        self.phrange_edit = QLineEdit(self)
        self.phrange_edit.setValidator(QRegExpValidator())
        self.phrange_edit.setAlignment(Qt.AlignRight)
        self.phrange_edit.setText('0, 14')

        tsteady_label, tsteady_unit = QLabel(self), QLabel(self)
        tsteady_label.setText('Plateau time'), tsteady_unit.setText('s')
        self.tsteady_edit = QLineEdit(self)
        self.tsteady_edit.setValidator(QDoubleValidator())
        self.tsteady_edit.setAlignment(Qt.AlignRight)
        self.tsteady_edit.setText('10.')

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
        self.ph_t90_edit.setText('0.1')

        ph_signal_label, ph_signal_unit = QLabel(self), QLabel(self)
        ph_signal_label.setText('Background signal')
        ph_signal_unit.setText('mV')
        self.ph_signal_edit = QLineEdit(self)
        self.ph_signal_edit.setValidator(QRegExpValidator())
        self.ph_signal_edit.setAlignment(Qt.AlignRight)
        self.ph_signal_edit.setText('5.0')

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
        self.nh3_t90_edit.setText('0.5')

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
        self.nh3_cGG_edit.setText('230.0')

        nh3_alpha_label, nh3_alpha_unit = QLabel(self), QLabel(self)
        nh3_alpha_label.setText('alpha(NH3)')
        nh3_alpha_unit.setText('%')
        self.nh3_alpha_edit = QLineEdit(self)
        self.nh3_alpha_edit.setValidator(QRegExpValidator())
        self.nh3_alpha_edit.setAlignment(Qt.AlignRight)
        self.nh3_alpha_edit.setText('0., 100., 0.01')

        # -----------------------
        # General navigation
        self.load_button = QPushButton('Load', self)
        self.load_button.setFixedWidth(100)
        self.inputFileLineEdit = QLineEdit(self)
        self.inputFileLineEdit.setValidator(QDoubleValidator())
        self.inputFileLineEdit.setMaximumWidth(300)
        self.inputFileLineEdit.setAlignment(Qt.AlignRight)
        self.plot_button = QPushButton('Plot', self)
        self.plot_button.setFixedWidth(100)
        self.clearP_button = QPushButton('Clear parameter', self)
        self.clearP_button.setFixedWidth(150)
        self.clearF_button = QPushButton('Clear plots', self)
        self.clearF_button.setFixedWidth(100)
        self.save_button = QPushButton('Save all', self)
        self.save_button.setFixedWidth(100)
        self.saveR_button = QPushButton('Save report', self)
        self.saveR_button.setFixedWidth(100)

        # draw additional "line" to separate load from plot and plot from save
        vline1 = QFrame()
        vline1.setFrameShape(QFrame.VLine | QFrame.Raised)
        vline1.setLineWidth(2)
        vline2 = QFrame()
        vline2.setFrameShape(QFrame.VLine | QFrame.Raised)
        vline2.setLineWidth(2)

        # -------------------------------------------------------------------------------------------
        # GroupBoxes to structure the layout
        # General navigation
        navigation_group = QGroupBox("Navigation Tool Bar")
        grid_load = QGridLayout()
        navigation_group.setFixedHeight(70)
        vbox_bottom.addWidget(navigation_group)
        navigation_group.setLayout(grid_load)

        grid_load.addWidget(self.load_button, 0, 0)
        grid_load.addWidget(self.inputFileLineEdit, 0, 1)
        grid_load.addWidget(vline1, 0, 2)
        grid_load.addWidget(self.plot_button, 0, 3)
        grid_load.addWidget(self.clearP_button, 0, 4)
        grid_load.addWidget(self.clearF_button, 0, 5)
        grid_load.addWidget(vline2, 0, 6)
        grid_load.addWidget(self.save_button, 0, 7)
        grid_load.addWidget(self.saveR_button, 0, 8)

        # ----------------------------------------------
        # create GroupBox to structure the layout
        general_group = QGroupBox("General Settings")
        general_group.setFixedWidth(250)
        grid_load = QGridLayout()
        grid_load.setSpacing(5), grid_load.setVerticalSpacing(1)

        # add GroupBox to layout and load buttons in GroupBox
        hbox_ltop.addWidget(general_group)
        general_group.setLayout(grid_load)
        grid_load.addWidget(temperature_label, 0, 0)
        grid_load.addWidget(self.temperature_edit, 0, 1)
        grid_load.addWidget(temperature_unit_label, 0, 2)
        grid_load.addWidget(pH_label, 1, 0)
        grid_load.addWidget(self.phrange_edit, 1, 1)
        grid_load.addWidget(tsteady_label, 2, 0)
        grid_load.addWidget(self.tsteady_edit, 2, 1)
        grid_load.addWidget(tsteady_unit, 2, 2)
        grid_load.addWidget(smprate_label, 3, 0)
        grid_load.addWidget(self.smpgrate_edit, 3, 1)
        grid_load.addWidget(smprate_unit_label, 3, 2)

        general_group.setContentsMargins(1, 15, 15, 1)
        hbox_ltop.addSpacing(10)

        # -----------------------
        # pH Sensor Settings
        phsens_group = QGroupBox("pH Sensor Settings")
        phsens_group.setFixedWidth(300)
        grid_load = QGridLayout()
        grid_load.setSpacing(5), grid_load.setVerticalSpacing(1)

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
        grid_load.setSpacing(5), grid_load.setVerticalSpacing(1)

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
        self.tsteady_edit.editingFinished.connect(self.print_tsteady)
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

        # ----------------------------------------------------------------------------------------------------------------
        # connect buttons in navigation manager with functions
        self.load_button.clicked.connect(self.load_data)
        self.plot_button.clicked.connect(self.run_simulation)
        self.clearP_button.clicked.connect(self.clear_parameters)
        self.clearF_button.clicked.connect(self.clear_phclaib)
        self.clearF_button.clicked.connect(self.clear_nh3calib)
        self.clearF_button.clicked.connect(self.clear_nh3timedrive)
        self.clearF_button.clicked.connect(self.clear_tantimdrive)
        self.save_button.clicked.connect(self.save)
        self.saveR_button.clicked.connect(self.save_report)

        # ----------------------------------------------------------------------------------------------------------------
        # Calibration of pH
        self.fig_phcalib, self.ax_phcalib = plt.subplots()
        self.canvas_phcalib = FigureCanvasQTAgg(self.fig_phcalib)
        self.navi_phcalib = NavigationToolbar2QT(self.canvas_phcalib, w, coordinates=False)
        self.ax_phcalib.set_xlim(-0, 15)
        self.ax_phcalib.set_xlabel('pH value')
        self.ax_phcalib.set_ylabel('Potential [mV]')
        self.fig_phcalib.subplots_adjust(left=0.15, right=0.95, bottom=0.2, top=0.9)
        sns.despine()

        # create GroupBox to structure the layout
        phcalib_group = QGroupBox("pH Sensor Calibration")
        phcalib_group.setMinimumWidth(220), phcalib_group.setMinimumHeight(320)
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
        self.ax_nh3calib.set_xlim(-5, 105.)
        self.ax_nh3calib.set_xlabel('alpha(NH$_3$) [%]')#, fontsize=fs)
        self.ax_nh3calib.set_ylabel('Potential [mV]')#, fontsize=fs)
        self.fig_nh3calib.subplots_adjust(left=0.15, right=0.95, bottom=0.2, top=0.9)
        sns.despine()

        nh3calib_group = QGroupBox("NH3 Sensor Calibration")
        nh3calib_group.setMinimumWidth(220), nh3calib_group.setMinimumHeight(320)
        grid_nh3calib = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_mbottom.addWidget(nh3calib_group)
        nh3calib_group.setLayout(grid_nh3calib)
        grid_nh3calib.addWidget(self.canvas_nh3calib)
        grid_nh3calib.addWidget(self.navi_nh3calib)

        # ---------------------------------------------------
        # nh3+nh4 simulation
        self.fig_nh3sim, self.ax_nh3sim = plt.subplots()
        self.ax_nh3sim_ph = self.ax_nh3sim.twinx()

        self.canvas_nh3sim = FigureCanvasQTAgg(self.fig_nh3sim)
        self.navi_nh3sim = NavigationToolbar2QT(self.canvas_nh3sim, w, coordinates=False)
        self.ax_nh3sim.set_xlabel('Time [s]')#, fontsize=fs)
        self.ax_nh3sim.set_ylabel('NH$_3$ / NH$_4^+$ [ppm]')#, fontsize=fs)

        self.fig_nh3sim.subplots_adjust(left=0.15, right=0.95, bottom=0.2, top=0.9)
        sns.despine()

        nh3sim_group = QGroupBox("NH3 / NH4+ Simulation")
        nh3sim_group.setMinimumWidth(220), nh3sim_group.setMinimumHeight(320)
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
        self.ax_tansim_ph = self.ax_tansim.twinx()
        self.ax_tansim.set_xlabel('Time [s]')#, fontsize=fs)
        self.ax_tansim.set_ylabel('TAN [ppm]')#, fontsize=fs)
        self.ax_tansim_ph.set_ylabel('pH', color='gray')#, fontsize=fs)
        self.fig_tanim.subplots_adjust(left=0.15, right=0.85, bottom=0.2, top=0.9)
        self.ax_tansim.spines['top'].set_visible(False), self.ax_tansim_ph.spines['top'].set_visible(False)

        tansim_group = QGroupBox("Total Ammonia Simulation")
        tansim_group.setMinimumWidth(220), tansim_group.setMinimumHeight(320)
        grid_tansim = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_rbottom.addWidget(tansim_group)
        tansim_group.setLayout(grid_tansim)
        grid_tansim.addWidget(self.canvas_tansim)
        grid_tansim.addWidget(self.navi_tansim)

        # -------------------------------------------------------------------------------------------------------------
        self.show()

    # -----------------------------------------------------------------------------------------------------------------
    # Functions for analysis
    # ------------------------------------------------------
    # print parameter
    def print_temperature(self):
        print('Temperature: ', self.temperature_edit.text(), 'degC')

    def print_phrange(self):
        print('pH range: ', self.phrange_edit.text())

    def print_tsteady(self):
        print('plateau time for step function: ', self.tsteady_edit.text())

    def print_samplingrate(self):
        print('Sampling rate: ', self.smpgrate_edit.text(), 's')

    def print_ph_t90(self):
        print('pH sensor response: ', self.ph_t90_edit.text(), 's')

    def print_ph_signal(self):
        print('pH sensor signal (min): ', self.ph_signal_edit.text(), 'mV')

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

    # ---------------------------------------------------
    def clear_parameters(self):
        # re-write default parameters
        self.temperature_edit.setText('25.')
        self.phrange_edit.setText('0, 14')
        self.tsteady_edit.setText('10.')
        self.smpgrate_edit.setText('1.')
        self.ph_t90_edit.setText('0.1')
        self.ph_signal_edit.setText('5.0')
        self.ph_res_edit.setText('1e-5')
        self.ph_ref_edit.setText('5e-3')
        self.nh3_t90_edit.setText('0.5')
        self.nh3_signal_edit.setText('0.02, 0.09')
        self.nh3_res_edit.setText('1e-9')
        self.nh3_pka_edit.setText('9.25')
        self.nh3_cGG_edit.setText('230.0')
        self.nh3_alpha_edit.setText('0., 100., 0.01')

    def clear_phclaib(self):
        self.ax_phcalib.cla()
        self.ax_phcalib.set_xlim(-0, 15)
        self.ax_phcalib.set_xlabel('pH value')#, fontsize=fs)
        self.ax_phcalib.set_ylabel('Potential [mV]')#, fontsize=fs)
        self.fig_phcalib.subplots_adjust(left=0.15, right=0.95, bottom=0.2, top=0.9)
        sns.despine()
        self.fig_phcalib.canvas.draw()

    def clear_nh3calib(self):
        self.ax_nh3calib.cla()
        self.ax_nh3calib.set_xlim(-5, 105.)
        self.ax_nh3calib.set_xlabel('alpha(NH$_3$) [%]')#, fontsize=fs)
        self.ax_nh3calib.set_ylabel('Potential [mV]')#, fontsize=fs)
        self.fig_nh3calib.subplots_adjust(left=0.15, right=0.95, bottom=0.2, top=0.9)
        sns.despine()
        self.fig_nh3calib.canvas.draw()

    def clear_nh3timedrive(self):
        self.ax_nh3sim.cla()
        self.ax_nh3sim.set_xlabel('Time [s]')#, fontsize=fs)
        self.ax_nh3sim.set_ylabel('NH$_3$ / NH$_4^+$ [ppm]')#, fontsize=fs)
        self.fig_nh3sim.subplots_adjust(left=0.15, right=0.95, bottom=0.2, top=0.9)
        sns.despine()
        self.fig_nh3sim.canvas.draw()

    def clear_tantimdrive(self):
        self.ax_tansim.cla()
        self.ax_tansim_ph.cla()
        self.ax_tansim.set_xlabel('Time [s]')#, fontsize=fs)
        self.ax_tansim.set_ylabel('TAN [ppm]')#, fontsize=fs)
        self.ax_tansim_ph.set_ylabel('pH', color='gray')#, fontsize=fs)
        self.fig_tanim.subplots_adjust(left=0.15, right=0.85, bottom=0.2, top=0.9)
        self.ax_tansim.spines['top'].set_visible(False), self.ax_tansim_ph.spines['top'].set_visible(False)
        self.ax_tansim_ph.spines['right'].set_visible(True)

        self.fig_tanim.canvas.draw()

    # ---------------------------------------------------
    def load_data(self):
        # opens a dialog window in the current path
        fname, filter = QFileDialog.getOpenFileName(self, "Select specific txt file for temperature compensation",
                                                    "", "Text files (*.txt *.csv *xls)")
        if fname:
            self.inputFileLineEdit.setText(fname)
            print('now do something with this file')
            df_general, df_ph, df_nh3 = bs.load_data(fname)

            print(df_general)
            # set parameter to run the simulation
            self.temperature_edit.setText(df_general.loc['Temperature', 'values'])
            s = df_nh3.loc['pH range', 'values'][1:-1]
            self.phrange_edit.setText(s)
            self.tsteady_edit.setText(df_general.loc['Plateau time', 'values'])
            self.smpgrate_edit.setText(df_general.loc['sampling rate'].values[0])

            self.ph_t90_edit.setText(df_ph.loc['t90', 'values'])
            self.ph_signal_edit.setText(df_ph.loc['background signal', 'values'])
            self.ph_res_edit.setText(df_ph.loc['resolution', 'values'])
            self.ph_ref_edit.setText(df_ph.loc['E0', 'values'])

            self.nh3_t90_edit.setText(df_nh3.loc['response time', 'values'])
            s = df_nh3.loc['signal min', 'values'] + ', ' + df_nh3.loc['signal max', 'values']
            self.nh3_signal_edit.setText(s)
            self.nh3_res_edit.setText(df_nh3.loc['resolution', 'values'])
            self.nh3_pka_edit.setText(df_nh3.loc['pKa', 'values'])
            self.nh3_cGG_edit.setText(df_general.loc['GGW concentration', 'values'])
            self.nh3_alpha_edit.setText(df_nh3.loc['nh3 range', 'values'][1:-1])

    def save(self):
        # opens window in current path. User input to define file name
        fname_save = QFileDialog.getSaveFileName(self, 'Save File')[0]

        if fname_save:
            if '.txt' in fname_save or '.csv' in fname_save:
                pass
            else:
                fname_save = fname_save + '.txt'

            # check whether we have data to save
            try:
                if self.dic_sens_record:
                    # save output now
                    output = bs.save_report(para_meas=self.para_meas, sensor_ph=self.sensor_ph, dtarget=self.dic_target,
                                            sensor_nh3=self.sensor_nh3, dsens_record=self.dic_sens_record)
                    output.to_csv(fname_save, sep='\t', header=None)

                    # save figures in separate folder
                    for f in self.dic_figures.keys():
                        figure_name = fname_save + '_Graph-' + '-'.join(f.split(' ')) + '.'
                        for t in save_type:
                            self.dic_figures[f].savefig(figure_name + t, dpi=300)
                else:
                    msgBox = QMessageBox()
                    msgBox.setIcon(QMessageBox.Information)
                    msgBox.setText("Simulate before saving")
                    msgBox.setWindowTitle("Warning")
                    msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

                    returnValue = msgBox.exec()
                    if returnValue == QMessageBox.Ok:
                        pass
            except NameError:
                msgBox = QMessageBox()
                msgBox.setIcon(QMessageBox.Information)
                msgBox.setText("Simulate before saving")
                msgBox.setWindowTitle("Warning")
                msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

                returnValue = msgBox.exec()
                if returnValue == QMessageBox.Ok:
                    pass

    def save_report(self):
        # opens window in current path. User input to define file name
        fname_save = QFileDialog.getSaveFileName(self, 'Save File')[0]

        if fname_save:
            if '.txt' in fname_save or '.csv' in fname_save:
                pass
            else:
                fname_save = fname_save + '.txt'

            # check whether we have data to save
            try:
                if self.dic_sens_record:
                    # save output now
                    output = bs.save_report(para_meas=self.para_meas, sensor_ph=self.sensor_ph, dtarget=self.dic_target,
                                            sensor_nh3=self.sensor_nh3, dsens_record=self.dic_sens_record)
                    output.to_csv(fname_save, sep='\t', header=None)

                else:
                    msgBox = QMessageBox()
                    msgBox.setIcon(QMessageBox.Information)
                    msgBox.setText("Simulate before saving")
                    msgBox.setWindowTitle("Warning")
                    msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

                    returnValue = msgBox.exec()
                    if returnValue == QMessageBox.Ok:
                        pass
            except NameError:
                msgBox = QMessageBox()
                msgBox.setIcon(QMessageBox.Information)
                msgBox.setText("Simulate before saving")
                msgBox.setWindowTitle("Warning")
                msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

                returnValue = msgBox.exec()
                if returnValue == QMessageBox.Ok:
                    pass

    # ---------------------------------------------------
    def run_simulation(self):
        self.sensor_ph = dict({'E0': float(self.ph_ref_edit.text()), 't90': float(self.ph_t90_edit.text()),
                               'resolution': float(self.ph_res_edit.text()), 'sensitivity': 2,
                               'time steps': float(self.smpgrate_edit.text()) / 1000,
                               'background signal': float(self.ph_signal_edit.text())})

        pH_range = [float(i) for i in self.phrange_edit.text().split(',')]
        self.sensor_nh3 = dict({'pH range': pH_range, 'sensitivity': 2,  'pKa': float(self.nh3_pka_edit.text()),
                                'response time': float(self.nh3_t90_edit.text()),
                                'resolution': float(self.nh3_res_edit.text()),
                                'time steps': float(self.smpgrate_edit.text()) / 1000,
                                'nh3 range': [float(i) for i in self.nh3_alpha_edit.text().split(',')],
                                'signal min': float(self.nh3_signal_edit.text().split(',')[0]),
                                'signal max': float(self.nh3_signal_edit.text().split(',')[1])})

        self.para_meas = dict({'Temperature': float(self.temperature_edit.text()),
                               'Plateau time': float(self.tsteady_edit.text()),
                               'pH steps': list(np.arange(pH_range[0], pH_range[1] + 1)),
                               'sampling rate': float(self.smpgrate_edit.text()),
                               'GGW concentration': float(self.nh3_cGG_edit.text())})

        # ------------------------------------------------------------------------------
        # TAN system - pH curve of NH3 and NH4 in percentage
        df_alpha = bs._tan_simulation(c_nh4=float(self.nh3_cGG_edit.text()), phmin=self.sensor_nh3['pH range'][0],
                                      phmax=self.sensor_nh3['pH range'][1], step_ph=1.,
                                      ph_deci=self.sensor_ph['sensitivity'], pKa=self.sensor_nh3['pKa'])

        # ------------------------------------------------------------------------------
        # pH sensor modelation
        dfpH_target, dfSig_calib, dfph_re = bs.pH_sensor(sensor_ph=self.sensor_ph, para_meas=self.para_meas,
                                                         df_alpha=df_alpha)

        # ------------------------------------------------------------------------------
        # NH3 / NH4+ sensor modulation
        [df_tan_target, dfconc_target, para_nh3,
         df_record, df_tan] = bs.NH3_sensor(df_alpha=df_alpha, dfph_re=dfph_re, para_meas=self.para_meas,
                                            dfpH_target=dfpH_target, sensor_nh3=self.sensor_nh3)

        # --------------------------------------------------------------------------------------------------------------
        # plotting part
        # single sensor calibration
        fig1 = plot_calibration_ph(dfSig_calib=dfSig_calib, fig=self.fig_phcalib, ax=self.ax_phcalib)

        conc_nh3 = np.arange(self.sensor_nh3['nh3 range'][0], self.sensor_nh3['nh3 range'][1],
                             step=self.sensor_nh3['nh3 range'][2])
        fig2 = plot_calibration_nh3(conc_nh3=conc_nh3, para_nh3=para_nh3, fig=self.fig_nh3calib, ax=self.ax_nh3calib)

        # final model
        fig3, fig4 = plot_tanModel(dfconc_target=dfconc_target, df_tan_target=df_tan_target, df_record=df_record,
                                   df_tan=df_tan, dfph_re=dfph_re, phmax=pH_range[1], ph_target=dfpH_target,
                                   fig1=self.fig_nh3sim, ax1=self.ax_nh3sim, fig2=self.fig_tanim, ax2=self.ax_tansim,
                                   ax22=self.ax_tansim_ph, ax12=self.ax_nh3sim_ph)

        # ------------------------------------------------------------------------------------------------------------------
        # collect for result output (save data)
        self.dic_target = dict({'TAN': df_tan_target, 'NH3 simulation': df_alpha, 'target conc nh3': dfconc_target})
        self.dic_sens_calib = dict({'pH': dfSig_calib, 'NH3': para_nh3})
        self.dic_sens_record = dict({'tan': df_tan, 'NH3': df_record, 'pH': dfph_re})
        self.dic_figures = dict({'pH calib': fig1, 'NH3 calib': fig2, 'model': fig3, 'TAN': fig4})


# .....................................................................................................................
def plot_calibration_ph(dfSig_calib, fig=None, ax=None):
    ax.cla()
    # preparation of figure plot
    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect('auto')
    ax.set_xlabel('pH value'), ax.set_ylabel('Potential [mV]')
    ax.plot([int(i) for i in dfSig_calib.index], dfSig_calib, lw=1., color='k')
    ax.set_xlim(0, 15)

    sns.despine(), fig.subplots_adjust(left=0.18, right=0.95, bottom=0.2, top=0.9)
    fig.canvas.draw()
    return fig


def plot_calibration_nh3(conc_nh3, para_nh3, fig=None, ax=None):
    ax.cla()
    # preparation of figure plot
    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect('auto')
    ax.set_xlabel('alpha(NH$_3$) [%]'), ax.set_ylabel('Potential [mV]')

    ax.plot(conc_nh3, para_nh3[0] * conc_nh3 + para_nh3[1], lw=1., color='k')

    sns.despine(), fig.subplots_adjust(left=0.18, right=0.95, bottom=0.2, top=0.9)
    fig.canvas.draw()
    return fig


def plot_tanModel(dfconc_target, df_tan_target, df_record, df_tan, dfph_re, phmax, ph_target, fig1=None, ax1=None,
                  ax12=None, fig2=None, ax2=None, ax22=None):
    ax1.cla(), ax2.cla(), ax22.cla()
    # preparation of figure plot
    if ax1 is None:
        fig1, ax1 = plt.subplots()
        ax1.set_aspect('auto')
        ax12 = ax1.twinx()
    if ax2 is None:
        fig2, ax2 = plt.subplots()
        ax2.set_aspect('auto')
        ax22 = ax2.twinx()

    ax1.spines['top'].set_visible(False), ax1.spines['right'].set_visible(False)
    ax22.spines['right'].set_visible(True), ax12.spines['right'].set_visible(True)

    ax1.set_xlabel('Time [s]'), ax2.set_xlabel('Time [s]')
    ax1.set_ylabel('NH$_3$ / NH$_4^+$ [ppm]'), ax12.set_ylabel('pH', color='gray')
    ax22.set_ylabel('pH', color='gray'), ax2.set_ylabel('TAN [ppm]', color=dcolor['TAN'])

    # top plot
    ax1.plot(df_record['nh3 / ppm'], lw=1., color=dcolor['NH3'], label='NH$_3$')
    ax1.plot(df_record['nh4 / ppm'], lw=1., color=dcolor['NH4'], label='NH$_4^+$')
    ax1.legend(frameon=True, fancybox=True, loc=0)
    ax1.plot(dfconc_target['nh3 / ppm'], lw=1., ls=ls['target'], color='k')
    ax1.plot(dfconc_target['nh4 / ppm'], lw=1., ls=ls['target'], color='gray')
    ax1.set_ylim(-10, dfconc_target['nh4 / ppm'].max()*1.1)
    ax12.plot(ph_target['pH'], lw=1., ls='-.', color='gray')

    # bottom plot
    ax2.plot(df_tan_target, lw=1., ls=ls['target'], color='k')
    ax2.plot(df_tan, lw=1., color=dcolor['TAN'])
    ax22.plot(ph_target['pH'], lw=1., ls='-.', color='gray')
    ax22.plot(dfph_re, lw=1., ls=':', color='k')
    ax2.set_ylim(df_tan['TAN'].loc[5:].min()*0.95, df_tan['TAN'].loc[5:].max()*1.05), ax22.set_ylim(-0.5, phmax * 1.05)

    fig1.subplots_adjust(left=0.15, right=0.95, bottom=0.2, top=0.9)
    fig2.subplots_adjust(left=0.15, right=0.85, bottom=0.2, top=0.9)
    fig1.canvas.draw(), fig2.canvas.draw()
    return fig1, fig2


# .....................................................................................................................
if __name__ == '__main__':
    app = QApplication(sys.argv)
    # where to display the GUI (which monitor in case there are several)
    view = MainWindow()
    # set size of the monitor (frame)
    view.setGeometry(50, 70, 1300, 750)
    sys.exit(app.exec_())
