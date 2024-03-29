__author__ = 'szieger'
__project__ = 'in silico study for sensor response'

import basics_sensorSignal_v3 as bs
import sys
import re
from PyQt5.QtGui import *
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QMainWindow, QPushButton, QAction, qApp,
                             QGridLayout, QLabel, QLineEdit, QGroupBox, QFileDialog, QFrame, QMessageBox, QCheckBox,
                             QDialog, QTableWidget, QTableWidgetItem)
from PyQt5.QtGui import QIcon, QDoubleValidator
from PyQt5.QtCore import Qt
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT, FigureCanvasQTAgg
import numpy as np
from scipy import integrate
import seaborn as sns
import pandas as pd
import os

# .....................................................................................................................
# global parameter
dcolor = dict({'pH': '#1CC49A', 'sig pH': '#4B5258', 'NH3': '#196E94', 'sig NH3': '#314945', 'NH4': '#DCA744',
               'sig NH4': '#89621A', 'TAN': '#A86349'})
ls = dict({'target': '-.', 'simulation': ':'})
save_type = ['png', 'svg']
sns.set_context('paper', font_scale=1.)

# fixed parameter depending on literature research
step_ph = 0.01
ph_deci = 2                # decimals for sensor sensitivity
ph_res = 1e-5              # resolution of the pH sensor

sig_max = 400              # maximal signal response at maximal NH4+ concentration in mV
sig_bgd = .5               # background signal / offset in mV at 0M NH4+
E0 = 0.43                  # zero potential of the reference electrode
tsteps = 1e-3              # time steps for pH and NH3 sensor (theory)

# electrochemical NH3 sensor
tan_calib = 1                # concentration point for TAN at which the potential was measured
sigNH3_max = 0.09            # maximal signal at maximal NH3 concentration in mV (pH=1)
sigNH3_bgd = 0.02            # background signal / offset in mV at 0M NH3
nh3_res = 1e-9               # resolution of the NH3 sensor
sbgd_nhx = 0.03              # background signal

# !!! TODO: clean up all out commented parts
# !!! TODO: finalize software to run it on all systems
# !!! TODO: double-check calculation in mg/L! does it make sense?
# !!! TODO: update save -> dataframe and include integral


# .....................................................................................................................
class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.initUI()
        self.df_res, self.dres = None, dict()
        self.setWindowIcon(QIcon('icon.png'))

        # -----------------------------------------------------------
        # for all parameters - connect LineEdit with function
        self.ph_t90_edit.returnPressed.connect(self.print_ph_t90)
        self.nh3_t90_edit.returnPressed.connect(self.print_nh3_t90)
        self.NH3_cbox.stateChanged.connect(self.NH3clickBox)
        self.NH4_cbox.stateChanged.connect(self.NH4clickBox)

        # connect buttons in navigation manager with functions
        self.load_button.clicked.connect(self.load_data)
        self.plot_button.clicked.connect(self.check_parameter)
        self.int_button.clicked.connect(self.calc_integral)
        self._integral_counter = 0
        self.clearP_button.clicked.connect(self.clear_parameters)
        self.clearF_button.clicked.connect(self.clear_phsim)
        self.clearF_button.clicked.connect(self.clear_nh3timedrive)
        self.clearF_button.clicked.connect(self.clear_tantimdrive)
        self.save_button.clicked.connect(self.save)
        self.saveR_button.clicked.connect(self.save_report)

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

        # ---------------------------------------------------------------------------------------
        # (Invisible) structure of main window (=grid)
        mlayout = QVBoxLayout(w)

        # 1st layer: box with vertical alignment to split in top and bottom
        # bottom: Navigation, middle: line, top: everything
        vbox_top, vbox_middle, vbox_bottom = QHBoxLayout(), QHBoxLayout(), QHBoxLayout()
        mlayout.addLayout(vbox_top), mlayout.addLayout(vbox_middle), mlayout.addLayout(vbox_bottom)

        # 2nd layer: box with vertical alignment to split in left and right (left: parameters, m1: line, right: plots)
        hbox_left, hbox_m1, hbox_right = QVBoxLayout(), QVBoxLayout(), QVBoxLayout()
        vbox_top.addLayout(hbox_left), vbox_top.addLayout(hbox_m1), vbox_top.addLayout(hbox_right)

        # 3rd layer: left box with parameter settings (top, middle, bottom)
        hbox_ltop, hbox_lmiddle, hbox_lbottom = QHBoxLayout(), QHBoxLayout(), QHBoxLayout()
        hbox_ltop.setContentsMargins(5, 10, 10, 5)
        hbox_left.addLayout(hbox_ltop), hbox_left.addLayout(hbox_lmiddle), hbox_left.addLayout(hbox_lbottom)

        # 4th layer: right box split horizontally (top: individual sensors bottom: TAN)
        hbox_rtop, hbox_rbottom = QHBoxLayout(), QHBoxLayout()
        hbox_right.addLayout(hbox_rtop), hbox_right.addLayout(hbox_rbottom)

        # 5th layer: right top box split vertically (left: pH right: NH3/NH4+)
        hbox_tright, hbox_tleft = QVBoxLayout(), QVBoxLayout()
        hbox_rtop.addLayout(hbox_tright), hbox_rtop.addLayout(hbox_tleft)

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
        tsteady_label, tsteady_unit = QLabel(self), QLabel(self)
        tsteady_label.setText('Plateau time'), tsteady_unit.setText('s')
        self.tsteady_edit = QLineEdit(self)
        self.tsteady_edit.setValidator(QDoubleValidator())
        self.tsteady_edit.setAlignment(Qt.AlignRight)
        self.tsteady_edit.setText('75.')

        # pH sensor settings
        pH_label = QLabel(self)
        pH_label.setText('pH concentration(s)')
        self.ph_edit = QLineEdit(self)
        self.ph_edit.setValidator(QRegExpValidator())
        self.ph_edit.setAlignment(Qt.AlignRight)
        self.ph_edit.setText('8.4, 10.9')

        ph_t90_label, ph_t90_unit = QLabel(self), QLabel(self)
        ph_t90_label.setText('Response time')
        ph_t90_unit.setText('s')
        self.ph_t90_edit = QLineEdit(self)
        self.ph_t90_edit.setValidator(QDoubleValidator())
        self.ph_t90_edit.setAlignment(Qt.AlignRight)
        self.ph_t90_edit.setText('30.')

        self.ph_drift_box, ph_drift_unit = QCheckBox('Sensor drift', self), QLabel(self)
        ph_drift_unit.setText('mV/s')
        self.ph_drift_edit = QLineEdit(self)
        self.ph_drift_edit.setValidator(QDoubleValidator())
        self.ph_drift_edit.setAlignment(Qt.AlignRight)
        self.ph_drift_edit.setText('0.01')

        # ---------------------
        # NH3 sensor settings
        nhx_label, nhx_unit = QLabel(self), QLabel(self)
        nhx_label.setText('c(TAN)'), nhx_unit.setText('mg/L')
        self.TAN_edit = QLineEdit(self)
        self.TAN_edit.setFixedWidth(100)
        self.TAN_edit.setValidator(QRegExpValidator())
        self.TAN_edit.setAlignment(Qt.AlignRight)
        self.TAN_edit.setText('100.0')

        # self.analyte = None
        self.NH3_cbox, self.NH4_cbox = QCheckBox('NH3', self), QCheckBox('NH4+', self)
        self.NH4_cbox.setChecked(True)

        nh3_t90_label, nh3_t90_unit = QLabel(self), QLabel(self)
        nh3_t90_label.setText('Response time')
        nh3_t90_unit.setText('s')
        self.nh3_t90_edit = QLineEdit(self)
        self.nh3_t90_edit.setFixedWidth(100)
        self.nh3_t90_edit.setValidator(QDoubleValidator())
        self.nh3_t90_edit.setAlignment(Qt.AlignRight)
        self.nh3_t90_edit.setText('60.')
        nh3_pka_label = QLabel(self)
        nh3_pka_label.setText('pKa')
        self.nh3_pka_edit = QLineEdit(self)
        self.nh3_pka_edit.setFixedWidth(100)
        self.nh3_pka_edit.setValidator(QDoubleValidator())
        self.nh3_pka_edit.setAlignment(Qt.AlignRight)
        self.nh3_pka_edit.setText('9.25')

        self.nh3_drift_box, nh3_drift_unit = QCheckBox('Sensor drift', self), QLabel(self)
        nh3_drift_unit.setText('mV/s')
        self.nh3_drift_edit = QLineEdit(self)
        self.nh3_drift_edit.setFixedWidth(100)
        self.nh3_drift_edit.setValidator(QDoubleValidator())
        self.nh3_drift_edit.setAlignment(Qt.AlignRight)
        self.nh3_drift_edit.setText('-0.1')

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
        self.int_button = QPushButton('Integral', self)
        self.int_button.setFixedWidth(100)
        self.int_button.setEnabled(False)
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
        vline3 = QFrame()
        vline3.setFrameShape(QFrame.VLine | QFrame.Raised)
        vline3.setLineWidth(2)

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
        grid_load.addWidget(self.int_button, 0, 4)
        grid_load.addWidget(vline2, 0, 5)
        grid_load.addWidget(self.clearP_button, 0, 6)
        grid_load.addWidget(self.clearF_button, 0, 7)
        grid_load.addWidget(vline3, 0, 8)
        grid_load.addWidget(self.save_button, 0, 9)
        grid_load.addWidget(self.saveR_button, 0, 10)

        # ----------------------------------------------
        # create GroupBox to structure the layout
        general_group = QGroupBox("General Settings")
        general_group.setFixedWidth(300)
        general_group.setFixedHeight(100)
        grid_load = QGridLayout()
        grid_load.setSpacing(5), grid_load.setVerticalSpacing(2)

        # add GroupBox to layout and load buttons in GroupBox
        hbox_ltop.addWidget(general_group)
        general_group.setLayout(grid_load)
        grid_load.addWidget(tsteady_label, 1, 0)
        grid_load.addWidget(self.tsteady_edit, 1, 1)
        grid_load.addWidget(tsteady_unit, 1, 2)

        general_group.setContentsMargins(1, 15, 15, 1)
        hbox_ltop.addSpacing(10)

        # -----------------------
        # pH Sensor Settings
        phsens_group = QGroupBox("pH Sensor Settings")
        phsens_group.setFixedWidth(300)
        phsens_group.setFixedHeight(150)
        grid_load = QGridLayout()
        grid_load.setSpacing(5), grid_load.setVerticalSpacing(1)

        # add GroupBox to layout and load buttons in GroupBox
        hbox_lmiddle.addWidget(phsens_group)
        phsens_group.setLayout(grid_load)

        grid_load.addWidget(pH_label, 0, 0)
        grid_load.addWidget(self.ph_edit, 0, 1)
        grid_load.addWidget(ph_t90_label, 1, 0)
        grid_load.addWidget(self.ph_t90_edit, 1, 1)
        grid_load.addWidget(ph_t90_unit, 1, 2)
        grid_load.addWidget(self.ph_drift_box, 2, 0)
        grid_load.addWidget(self.ph_drift_edit, 2, 1)
        grid_load.addWidget(ph_drift_unit, 2, 2)

        phsens_group.setContentsMargins(1, 15, 15, 1)
        hbox_lmiddle.addSpacing(10)

        # -----------------------
        # NH3 Sensor Settings
        nh3sens_group = QGroupBox("NH3/NH4+ Sensor Settings")
        nh3sens_group.setFixedWidth(300)
        nh3sens_group.setFixedHeight(200)
        grid_load = QGridLayout()
        grid_load.setSpacing(5), grid_load.setVerticalSpacing(1)

        # add GroupBox to layout and load buttons in GroupBox
        hbox_lbottom.addWidget(nh3sens_group)
        nh3sens_group.setLayout(grid_load)

        grid_load.addWidget(nhx_label, 0, 0)
        grid_load.addWidget(self.TAN_edit, 0, 1)
        grid_load.addWidget(nhx_unit, 0, 2)
        grid_load.addWidget(self.NH3_cbox, 1, 1)
        grid_load.addWidget(self.NH4_cbox, 1, 2)
        grid_load.addWidget(nh3_t90_label, 2, 0)
        grid_load.addWidget(self.nh3_t90_edit, 2, 1)
        grid_load.addWidget(nh3_t90_unit, 2, 2)
        grid_load.addWidget(nh3_pka_label, 3, 0)
        grid_load.addWidget(self.nh3_pka_edit, 3, 1)
        grid_load.addWidget(self.nh3_drift_box, 4, 0)
        grid_load.addWidget(self.nh3_drift_edit, 4, 1)
        grid_load.addWidget(nh3_drift_unit, 4, 2)
        
        nh3sens_group.setContentsMargins(1, 15, 15, 1)
        hbox_lbottom.addSpacing(10)

        # ----------------------------------------------------------------------------------------------------------------
        # pH Simulation
        self.fig_phsim, self.ax_phsim = plt.subplots()
        self.canvas_phsim = FigureCanvasQTAgg(self.fig_phsim)
        self.ax_phsim.set_xlabel('Time / s')
        self.ax_phsim.set_ylabel('pH value')
        self.fig_phsim.tight_layout(pad=1, rect=(0.05, 0.1, 1, 0.98))
        sns.despine()

        # create GroupBox to structure the layout
        phsim_group = QGroupBox("pH Sensor Simulation")
        grid_phsim = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_tright.addWidget(phsim_group)
        phsim_group.setLayout(grid_phsim)
        grid_phsim.addWidget(self.canvas_phsim)

        # ---------------------------------------------------
        # NH3+NH4 simulation
        self.fig_nh3sim, self.ax_nh3sim = plt.subplots()
        self.ax1_nh3sim = self.ax_nh3sim.twinx()
        self.canvas_nh3sim = FigureCanvasQTAgg(self.fig_nh3sim)
        # self.navi_nh3sim = NavigationToolbar2QT(self.canvas_nh3sim, w, coordinates=False)
        self.ax_nh3sim.set_xlabel('Time / s')
        self.ax_nh3sim.set_ylabel('NH$_4^+$ / mg/L', color=dcolor['NH4'])
        self.ax1_nh3sim.set_ylabel('NH$_3$ / mg/L', color=dcolor['NH3'])
        self.ax_nh3sim.spines['top'].set_visible(False), self.ax1_nh3sim.spines['top'].set_visible(False)
        self.fig_nh3sim.tight_layout(pad=1, rect=(0.05, 0.1, 0.95, 0.98))

        nh3sim_group = QGroupBox("NH3 / NH4+ Simulation")
        grid_nh3sim = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_tleft.addWidget(nh3sim_group)
        nh3sim_group.setLayout(grid_nh3sim)
        grid_nh3sim.addWidget(self.canvas_nh3sim)
        # grid_nh3sim.addWidget(self.navi_nh3sim)

        # TAN simulation
        self.fig_tansim, self.ax_tansim = plt.subplots()
        self.canvas_tansim = FigureCanvasQTAgg(self.fig_tansim)
        # self.navi_tansim = NavigationToolbar2QT(self.canvas_tansim, w, coordinates=False)
        self.ax_tansim.set_xlabel('Time / s')
        self.ax_tansim.set_ylabel('TAN / mg/L')
        self.fig_tansim.tight_layout(pad=.75, rect=(0, 0.1, 1, .95))
        sns.despine()

        tansim_group = QGroupBox("Total Ammonia Simulation")
        grid_tansim = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_rbottom.addWidget(tansim_group)
        tansim_group.setLayout(grid_tansim)
        grid_tansim.addWidget(self.canvas_tansim)
        # grid_tansim.addWidget(self.navi_tansim)

        # -------------------------------------------------------------------------------------------------------------
        self.show()

    # ---------------------------------------------------
    def print_ph_t90(self):
        print('pH sensor response: ', self.ph_t90_edit.text(), 's')

        # in case set response time equals 0, return warning
        if float(self.ph_t90_edit.text()) == 0:
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Information)
            msgBox.setText("Please enter a valid sensor response time for the pH sensor, notably a number > 0")
            msgBox.setFont(QFont('Helvetica Neue', 11))
            msgBox.setWindowTitle("Warning")
            msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

            returnValue = msgBox.exec()
            if returnValue == QMessageBox.Ok:
                pass

    def print_nh3_t90(self):
        print('NH3 sensor response: ', self.nh3_t90_edit.text(), 's')
        # in case set response time equals 0, return warning
        if float(self.nh3_t90_edit.text()) == 0:
            msgBox = QMessageBox()
            msgBox.setIcon(QMessageBox.Information)
            msgBox.setText("Please enter a valid sensor response time for the TAN sensor, notably a number > 0")
            msgBox.setFont(QFont('Helvetica Neue', 11))
            msgBox.setWindowTitle("Warning")
            msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

            returnValue = msgBox.exec()
            if returnValue == QMessageBox.Ok:
                pass

    # ---------------------------------------------------
    def clear_parameters(self):
        # re-write default parameters
        self.tsteady_edit.setText('75.')
        self.ph_edit.setText('8.4, 10.9')
        self.ph_t90_edit.setText('30.')
        self.TAN_edit.setText('100.0')
        self.NH3_cbox.setChecked(False)
        self.NH4_cbox.setChecked(True)
        self.nh3_t90_edit.setText('60.')
        self.nh3_pka_edit.setText('9.25')
        self.ph_drift_edit.setText('0.01')
        self.nh3_drift_edit.setText('-0.1')
        self.inputFileLineEdit.setText('')

    def clear_phsim(self):
        self.ax_phsim.cla()
        self.ax_phsim.set_xlabel('Time / s')
        self.ax_phsim.set_ylabel('pH value')
        self.fig_phsim.tight_layout(pad=0.6, rect=(0., 0.015, 1, 0.98))
        sns.despine()
        self.fig_phsim.canvas.draw()

    def clear_nh3timedrive(self):
        self.ax_nh3sim.cla()
        self.ax1_nh3sim.cla()
        self.ax_nh3sim.set_xlabel('Time / s')
        self.ax_nh3sim.set_ylabel('NH$_4^+$ / mg/L', color=dcolor['NH4'])
        self.ax1_nh3sim.set_ylabel('NH$_3$ / mg/L', color=dcolor['NH3'])
        self.fig_nh3sim.tight_layout(pad=0., rect=(0.015, 0.05, 0.98, 0.97))
        self.ax_nh3sim.spines['top'].set_visible(False), self.ax1_nh3sim.spines['top'].set_visible(False)
        self.fig_nh3sim.canvas.draw()

    def clear_tantimdrive(self):
        self.ax_tansim.cla()
        self.ax_tansim.set_xlabel('Time / s')
        self.ax_tansim.set_ylabel('TAN / mg/L')
        self.fig_tansim.tight_layout(pad=.75, rect=(0.01, 0.01, 1, .95))
        sns.despine()

        self.fig_tansim.canvas.draw()

    def NH3clickBox(self, state):
        if state == Qt.Checked:
            self.NH4_cbox.setCheckState(False)

    def NH4clickBox(self, state):
        if state == Qt.Checked:
            self.NH3_cbox.setCheckState(False)

    def align_concentrations(self, ls_ph, ls_cTAN_ppm):
        if len(ls_cTAN_ppm) == len(ls_ph):
            pass
        elif len(ls_cTAN_ppm) == 1 and len(ls_ph) > 1:
            ls_cTAN_ppm = ls_cTAN_ppm * len(ls_ph)
        elif len(ls_cTAN_ppm) > 1 and len(ls_ph) == 1:
            ls_ph = ls_ph * len(ls_cTAN_ppm)
        else:
            if max(len(ls_ph), len(ls_cTAN_ppm)) % min(len(ls_ph), len(ls_cTAN_ppm)) == 0:
                print('double the length of each other')
                if len(ls_ph) < len(ls_cTAN_ppm):
                    ls_ph = ls_ph * int(len(ls_cTAN_ppm) / 2)
                else:
                    ls_cTAN_ppm = ls_cTAN_ppm * int(len(ls_ph) / 2)
            else:
                if len(ls_ph) < len(ls_cTAN_ppm):
                    ls_ph.append(ls_ph[-1])
                else:
                    ls_cTAN_ppm.append(ls_cTAN_ppm[-1])
        return ls_ph, ls_cTAN_ppm

    def parameter_prep(self, t_plateau, ls_ph, ls_tan, pKa):
        # target fluctuation of analytes
        target_ph = bs.target_fluctuation(ls_conc=ls_ph, tstart=0, tstop=t_plateau * 2, nP=1, analyte='pH')
        target_tan = bs.target_fluctuation(ls_conc=ls_tan, tstart=0, tstop=t_plateau * 2, nP=1, analyte='TAN')

        # target concentration NH3/NH4 based on target pH
        nh3 = pd.DataFrame(bs.henderson_TAN4NH3(pH=target_ph['signal pH'].to_numpy(),
                                                c_tan=target_tan['signal TAN'].to_numpy(), pKa=pKa),
                           index=target_ph.index, columns=['signal NH3'])
        nh4 = pd.DataFrame(bs.henderson_TAN4NH4(pH=target_ph['signal pH'].to_numpy(),
                                                c_tan=target_tan['signal TAN'].to_numpy(), pKa=pKa),
                           index=target_ph.index, columns=['signal NH4'])

        # calculate target concentration of individual NHx parameters
        df_res = pd.concat([target_ph, target_tan, nh3, nh4], axis=1)

        return df_res

    # ---------------------------------------------------
    def load_data(self):
        # opens a dialog window in the current path
        fname, filter = QFileDialog.getOpenFileName(self, "Select specific txt file for temperature compensation",
                                                    "", "Text files (*.txt *.csv *xls)")

        if fname:
            self.inputFileLineEdit.setText(fname)
            df_general, df_ph, df_nh3 = bs.load_data(fname)

            # set parameter to run the simulation
            # self.temperature_edit.setText(df_general.loc['temperature', 'values'])
            self.tsteady_edit.setText(df_general.loc['plateau time', 'values'])

            # pH sensor
            self.ph_t90_edit.setText(df_ph.loc['t90', 'values'])
            self.ph_drift_edit.setText(df_ph.loc['drift', 'values'])
            ph_target = df_ph.loc['pH target', 'values'][1:-1]
            self.ph_edit.setText(ph_target)

            # NHx sensor
            self.nh3_t90_edit.setText(df_nh3.loc['t90', 'values'])
            self.nh3_pka_edit.setText(df_nh3.loc['pKa', 'values'])
            self.nh3_drift_edit.setText(df_nh3.loc['drift', 'values'])

            self.TAN_edit.setText(df_nh3.loc['TAN target'].values[0][0][1:-1])
            if df_nh3.loc['analyte', 'values'] == 'NH3':
                self.NTAN_editH3_cbox.setCheckState(True)
                self.NH4_cbox.setCheckState(False)
            else:
                self.NH4_cbox.setCheckState(True)
                self.NH3_cbox.setCheckState(False)

            # select status of drift
            if float(self.nh3_drift_edit.text()) != 0:
                self.nh3_drift_box.setChecked(True)
            if float(self.ph_drift_edit.text()) != 0:
                self.ph_drift_box.setChecked(True)

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
                if self.df_res is None:
                    msgBox = QMessageBox()
                    msgBox.setIcon(QMessageBox.Information)
                    msgBox.setText("Simulate before saving. Simulation might have been unsuccessful.")
                    msgBox.setWindowTitle("Warning")
                    msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

                    returnValue = msgBox.exec()
                    if returnValue == QMessageBox.Ok:
                        pass
                else:
                    # save output now
                    output = bs.save_report(para_meas=self.para_meas, sensor_ph=self.sensor_ph, df_res=self.df_res,
                                            sensor_nh3=self.sensor_nh3, dres=self.dres)
                    with open(fname_save, 'w') as f:
                        output[0].to_csv(f)
                    with open(fname_save, 'a') as f:
                        output[1].to_csv(f, header=False)

                    # save figures in separate folder
                    for f in self.dic_figures.keys():
                        figure_name = fname_save.split('.')[0] + '_Graph-' + '-'.join(f.split(' ')) + '.'
                        for t in save_type:
                            self.dic_figures[f].savefig(figure_name + t, dpi=300)

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
                if self.df_res is None:
                    msgBox = QMessageBox()
                    msgBox.setIcon(QMessageBox.Information)
                    msgBox.setText("Simulate before saving")
                    msgBox.setWindowTitle("Warning")
                    msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

                    returnValue = msgBox.exec()
                    if returnValue == QMessageBox.Ok:
                        pass
                else:
                    # convert to DataFrame
                    int_out = pd.DataFrame.from_dict(self.dres)
                    int_out.index = ['Integration range (s)', 'target pH', 'target TAN', 'integral target TAN',
                                     'integral observed TAN', 'error']
                    # save output
                    output = bs.save_report(para_meas=self.para_meas, sensor_ph=self.sensor_ph, df_res=self.df_res,
                                            sensor_nh3=self.sensor_nh3, dres=self.dres)
                    # output.to_csv(fname_save, sep='\t', header=None)
                    with open(fname_save, 'w') as f:
                        output[0].to_csv(f)
                    with open(fname_save, 'a') as f:
                        output[1].to_csv(f, header=False)

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
    def _linEdit2list(self, line):
        ls = list()
        if '.' in line and ',' in line:
            ls = [float(i) for i in line.split(',')]
        elif '.' in line:
            # double-check only one '.' in str
            if len(re.findall(r"\.", line)) == 1:
                ls.append(float(line))
            else:
                if len(line.split(' ')) != 1:
                    ls = [float(i) for i in line.split(' ')]
                else:
                    raise ValueError('Typo error. Please check your pH entry: ' + line)
        elif ',' in line:
            ls = [float(i) for i in line.split(',')]
        else:
            ls.append(float(line))
        return ls

    def check_parameter(self):
        # identify NHx analyte
        if self.NH3_cbox.isChecked():
            analyte, additional = 'NH3', 'NH4'
        elif self.NH4_cbox.isChecked():
            analyte, additional = 'NH4', 'NH3'

        # get TAN concentration(s) and pH value(s) and align list of fluctuation points
        ls_cTAN_ppm = self._linEdit2list(line=self.TAN_edit.text())
        ls_ph = self._linEdit2list(line=self.ph_edit.text())
        ls_ph, ls_cTAN_ppm = self.align_concentrations(ls_ph=ls_ph, ls_cTAN_ppm=ls_cTAN_ppm)

        # determine all required target concentrations
        df_target = self.parameter_prep(ls_ph=ls_ph, ls_tan=ls_cTAN_ppm, pKa=float(self.nh3_pka_edit.text()),
                                        t_plateau=float(self.tsteady_edit.text()))

        # check-box status
        if self.ph_drift_box.isChecked():
            drift1 = float(self.ph_drift_edit.text())
        else:
            drift1 = 0
        if self.nh3_drift_box.isChecked():
            drift2 = float(self.nh3_drift_edit.text())
        else:
            drift2 = 0

        # collect all relevant parameter
        self.sensor_ph = dict({'E0': E0, 't90': float(self.ph_t90_edit.text()), 'resolution': ph_res, 'drift': drift1,
                               'time steps': tsteps, 'start pH': ls_ph[0], 'sensitivity': ph_deci,
                               'pH target': df_target['signal pH'], 'set values': ls_ph})
        self.sensor_nh3 = dict({'sensitivity': ph_deci, 'pKa': float(self.nh3_pka_edit.text()), 'time steps': tsteps,
                                't90': float(self.nh3_t90_edit.text()), 'drift': drift2, 'resolution': nh3_res,
                                'NH3 target': df_target['signal NH3'], 'NH4 target': df_target['signal NH4'],
                                'TAN target': df_target['signal TAN'], 'calib TAN': tan_calib, 'analyte': analyte,
                                'additional': additional, 'signal min': sigNH3_bgd, 'signal max': sigNH3_max,
                                'background signal': sbgd_nhx, 'set values': ls_cTAN_ppm})
        self.para_meas = dict({#'temperature': float(self.temperature_edit.text()),
                               'plateau time': float(self.tsteady_edit.text())})

        # -----------------------------------------------------------------------------------------------------------
        # check pH sensor
        try:
            if float(self.sensor_ph['t90']) == 0:
                msgBox = QMessageBox()
                msgBox.setIcon(QMessageBox.Critical)
                msgBox.setText("Please enter a valid sensor response time for the pH sensor, in particular a number "
                               "> 0")
                msgBox.setFont(QFont('Helvetica Neue', 11))
                msgBox.setWindowTitle("Error")
                msgBox.exec_()
            else:
                try:
                    # check NH3/NH4 sensor
                    if float(self.sensor_nh3['t90']) == 0:
                        msgBox = QMessageBox()
                        msgBox.setIcon(QMessageBox.Critical)
                        msgBox.setText("Please enter a valid sensor response time for the TAN sensor, in particular a "
                                       "number > 0")
                        msgBox.setFont(QFont('Helvetica Neue', 11))
                        msgBox.setWindowTitle("Error")
                        msgBox.exec_()
                    else:
                        # in case parameter check was successful, start run_simulation
                        self.run_simulation()
                except:
                    pass
        except:
            pass

    def run_simulation(self):
        # clear figures and xcorrds
        self.clear_phsim(), self.clear_nh3timedrive(), self.clear_tantimdrive()
        self.ls_xcoords = []

        # target concentrations
        df_res = pd.concat([self.sensor_ph['pH target'], self.sensor_nh3['NH3 target'], self.sensor_nh3['NH4 target'],
                            self.sensor_nh3['TAN target']], axis=1)
        df_res.columns = ['target pH', 'target_mg/L NH3', 'target_mg/L NH4', 'target_mg/L TAN']

        # individual sensor - reduce by individual plotting (different function)
        [df_res, cplateaupH, cplateauTAN,
         para_nhx] = bs._alignSensorSettings(sensor_ph=self.sensor_ph, sensor_nh3=self.sensor_nh3, df_target=df_res,
                                             para_meas=self.para_meas)
        xnew = [round(u, 2) for u in df_res.index]
        df_res.index = xnew
        df_res = df_res.dropna()

        # pH sensor response
        [df_pHrec, df_pHdrift, df_pHcalc] = bs.pH_sensor(cplateau=cplateaupH, sensor_ph=self.sensor_ph,
                                                         para_meas=self.para_meas)
        xnew = [round(i, 2) for i in df_pHdrift.index]
        df_pHdrift.index, df_pHcalc.index = xnew, xnew

        # include sensor response and drift to result dataframe
        df_res = pd.concat([df_res, df_pHdrift.loc[df_res.index], df_pHcalc.loc[df_res.index]],
                           axis=1).sort_index(axis=1)

        # NHx sensor - calculate individual sensor response as well as TAN
        df_res = bs.NHxSensorcalc(df_res=df_res, analyte=self.sensor_nh3['analyte'], cplateauTAN=cplateauTAN,
                                  df_pHcalc=df_pHcalc, sensor_nh3=self.sensor_nh3, para_meas=self.para_meas,
                                  sensor_ph=self.sensor_ph)
        self.df_res = df_res.sort_index(axis=1)

        # --------------------------------------------------------------------------------------------------------------
        # plotting part | individual sensor - target vs record
        fig_pH = plot_phsensor(df_res=self.df_res, fig=self.fig_phsim, ax=self.ax_phsim)
        fig_NHx = plot_NHxsensor(df_res=self.df_res, fig=self.fig_nh3sim, ax=self.ax_nh3sim, ax1=self.ax1_nh3sim)

        # final TAN model
        fig_tan = plot_tanModel(df_res=self.df_res, fig1=self.fig_tansim, ax1=self.ax_tansim)
        # allow click events for TAN for integration (simpson for discrete measurement data)
        self.fig_tansim.canvas.mpl_connect('button_press_event', self.onclick_integral)

        # --------------------------------------------------------------------------------------------------------------
        # collect for result output (save data)
        self.dic_figures = dict({'pH': fig_pH, 'NH3': fig_NHx, 'TAN': fig_tan})

    def onclick_integral(self, event):
        modifiers = QApplication.keyboardModifiers()
        if modifiers != Qt.ControlModifier:  # change selected range
            return
        if len(self.ls_xcoords) >= 2:
            self.ls_xcoords = list()

        if event.xdata is None:
            if len(self.ls_xcoords) < 2:
                event.xdata = self.df_res.index.max()
            else:
                event.xdata = 0

        self.ls_xcoords.append(event.xdata)
        self.ax_tansim.vlines(x=self.ls_xcoords, ymin=self.df_res['TAN calc'].min(), ymax=self.df_res['TAN calc'].max(),
                              color='k', lw=0.15)
        if len(self.ls_xcoords) == 2:
            self.ax_tansim.axvspan(self.ls_xcoords[0], self.ls_xcoords[1], color='grey', alpha=0.3)

        # update figure plot
        self.fig_tansim.canvas.draw()

        # allow integral calculation as soon as xcoords have been collected
        self.int_button.setEnabled(True)

    def calc_integral(self):
        self._integral_counter += 1
        try:
            self.ls_xcoords
        except:
            self.ls_xcoords = list()

        # calculate integral for calculated/target TAN using simpson's rule for discrete data
        # subtract calculated TAN by target TAN
        if len(self.ls_xcoords) == 2:
            self.df4int = self.df_res[['TAN calc', 'target_mg/L TAN']].loc[min(self.ls_xcoords):max(self.ls_xcoords)]
            dfint_calc = integrate.simpson(self.df4int['TAN calc'].to_numpy(), x=self.df4int.index, dx=0.05, axis=- 1)
            dfint_target = integrate.simpson(self.df4int['target_mg/L TAN'].to_numpy(), x=self.df4int.index, dx=0.05,
                                             axis=- 1)

            # Integration range (s), target pH, target TAN, integral target TAN, integral TAN, error
            result = list([(round(min(self.ls_xcoords),2), round(max(self.ls_xcoords)),2),
                           self.df_res['target pH'].mean(), self.df_res['target_mg/L TAN'].mean(), dfint_calc,
                           dfint_target, dfint_calc - dfint_target])

            # add current results to dictionary (where all selected peak information are stored)
            self.dres[self._integral_counter] = result

            # open a pop up window with options to select what shall be saved
            global wInt
            wInt = IntegralWindow(self.dres, self.para_meas, self.sensor_ph, self.sensor_nh3, self._integral_counter)
            if wInt.isVisible() is False:
                wInt.show()


class IntegralWindow(QDialog):
    def __init__(self, res_int, para_meas, sensor_ph, sensor_nh3, count):
        super().__init__()
        self.result, self.count = res_int, count
        self.para_meas, self.sensor_ph, self.sensor_nh3 = para_meas, sensor_ph, sensor_nh3
        self.initUI()

        # when checkbox selected, save information in registered field
        self.reset_button.clicked.connect(self.reset)
        self.save_button.clicked.connect(self.save_integral)
        self.close_button.clicked.connect(self.close_window)

        # execute function as soon as the window pups up or as soon as integral_button is clicked
        self.res2table()

    def initUI(self):
        self.setWindowTitle("calculation error via integration")
        self.setGeometry(150, 180, 700, 300) # x, y, width, height

        # close window button
        self.close_button = QPushButton('OK', self)
        self.close_button.setFixedWidth(100), self.close_button.setFont(QFont('Helvetica Neue', 11))
        self.reset_button = QPushButton('Reset', self)
        self.reset_button.setFixedWidth(100), self.reset_button.setFont(QFont('Helvetica Neue', 11))
        self.save_button = QPushButton('Save', self)
        self.save_button.setFixedWidth(100), self.save_button.setFont(QFont('Helvetica Neue', 11))

        # create table to store data
        self.tab_report = QTableWidget(self)
        self.tab_report.setColumnCount(6), self.tab_report.setRowCount(1)
        self.tab_report.setHorizontalHeaderLabels(['Integration range (s)', 'target pH', 'target TAN',
                                                   'integral target TAN', 'integral observed TAN', 'error'])
        self.tab_report.resizeColumnsToContents()
        self.tab_report.resizeRowsToContents()

        # creating window layout
        mlayout2 = QVBoxLayout()
        vbox2_top, vbox2_middle, vbox2_bottom = QHBoxLayout(), QHBoxLayout(), QHBoxLayout()
        mlayout2.addLayout(vbox2_top), mlayout2.addLayout(vbox2_middle), mlayout2.addLayout(vbox2_bottom)

        # panel for integration results
        table_grp = QGroupBox("Results")
        grid_res = QGridLayout()
        table_grp.setFont(QFont('Helvetica Neue', 10))
        vbox2_top.addWidget(table_grp)
        table_grp.setLayout(grid_res)

        grid_res.addWidget(self.tab_report)
        self.tab_report.adjustSize()

        # vertical line
        hline = QFrame()
        hline.setFrameShape(QFrame.HLine | QFrame.Raised)
        hline.setLineWidth(2)
        vbox2_middle.addWidget(hline)

        # navigation panel
        navi_grp = QGroupBox("Navigation")
        grid_navi = QGridLayout()
        navi_grp.setFont(QFont('Helvetica Neue', 10))
        vbox2_bottom.addWidget(navi_grp)
        navi_grp.setLayout(grid_navi)

        # include widgets in the layout
        grid_navi.addWidget(self.close_button, 0, 0)
        grid_navi.addWidget(self.reset_button, 0, 1)
        grid_navi.addWidget(self.save_button, 0, 2)

        # add everything to the window layout
        self.setLayout(mlayout2)

    def res2table(self):
        # get the current row by counting the clicks on the integral button
        x0 = self.count-1

        # check whether row is empty (first column)
        item = self.tab_report.item(x0, 0)
        x = x0 if not item or not item.text() else x0 + 1

        # add the number of rows according to keys in dictionary
        self.tab_report.setRowCount(len(self.result.keys()))

        # go through the dictionary and fill up the table (again)
        for c in self.result.keys():
            # columns: Integration range (s), target pH, target TAN, integral target TAN, integral TAN, error
            for en in enumerate(self.result[c]):
                if en[0] == 0:
                    l = str(round(en[1][0], 2)) + ' - ' + str(round(en[1][1], 2))
                else:
                    l = str(round(en[1], 2))

                item = QTableWidgetItem(l)
                item.setTextAlignment(Qt.AlignRight)

                # item structure: row, table, content
                self.tab_report.setItem(c-1, en[0], item)
        self.tab_report.resizeColumnsToContents(), self.tab_report.resizeRowsToContents()

    def reset(self):
        self.tab_report.clear()
        self.tab_report.setHorizontalHeaderLabels(['Integration range (s)', 'target pH', 'target TAN',
                                                   'integral target TAN', 'integral observed TAN', 'error'])
        self.tab_report.setRowCount(0)

    def save_integral(self):
        # opens window in current path. User input to define file name
        fname_save = QFileDialog.getSaveFileName(self, 'Save File')[0]

        if fname_save:
            if '.txt' in fname_save or '.csv' in fname_save:
                fname_save = fname_save.split('.')[0] + '_integral.txt'
            else:
                fname_save = fname_save + '_integral.txt'

            # check whether we have data to save
            try:
                if self.result is None:
                    msgBox = QMessageBox()
                    msgBox.setIcon(QMessageBox.Information)
                    msgBox.setText("Simulate before saving. Simulation might have been unsuccessful.")
                    msgBox.setWindowTitle("Warning")
                    msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

                    returnValue = msgBox.exec()
                    if returnValue == QMessageBox.Ok:
                        pass
                else:
                    # save output now
                    output = bs.save_integral(para_meas=self.para_meas, sensor_ph=self.sensor_ph, dres=self.result,
                                              sensor_nh3=self.sensor_nh3)
                    with open(fname_save, 'w') as f:
                        output.to_csv(f)

            except NameError:
                msgBox = QMessageBox()
                msgBox.setIcon(QMessageBox.Information)
                msgBox.setText("Simulate before saving")
                msgBox.setWindowTitle("Warning")
                msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

                returnValue = msgBox.exec()
                if returnValue == QMessageBox.Ok:
                    pass

    def close_window(self):
        self.hide()


# .....................................................................................................................
def plot_phsensor(df_res, fig=None, ax=None):
    # preparation of figure plot
    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect('auto')
    else:
        ax.cla()
    ax.set_xlabel('Time / s'), ax.set_ylabel('pH value')

    ax.plot(df_res['target pH'].dropna(), ls='-.', lw=1., color='gray')
    ax.plot(df_res['pH calc'].dropna(), color=dcolor['pH'])

    ax.set_xlim(-0.5, df_res['pH calc'].dropna().index[-1]*1.05)
    sns.despine(), fig.tight_layout(pad=0.6, rect=(0., 0.015, 1, 0.98))
    fig.canvas.draw()
    return fig


def plot_NHxsensor(df_res, fig=None, ax=None, ax1=None):
    # preparation of figure plot
    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect('auto')
        ax1 = ax.twinx()
    else:
        ax.cla(), ax1.cla()
    ax.set_xlabel('Time / s')
    ax.set_ylabel('NH$_4^+$ / mg/L', color=dcolor['NH4'])
    ax1.set_ylabel('NH$_3$ / mg/L', color=dcolor['NH3'])

    ax1.plot(df_res['target_mg/L NH3'].dropna(), lw=1., ls=':', color='gray', label='NH$_3$ target')
    ax1.plot(df_res['NH3 calc'].dropna(), lw=1., color=dcolor['NH3'], label='NH$_3$')

    ax.plot(df_res['target_mg/L NH4'].dropna(), lw=1., ls='--', color='grey', label='NH$_4^+$ target')
    ax.plot(df_res['NH4 calc'].dropna(), lw=1., color=dcolor['NH4'], label='NH$_4^+$')

    ax.set_xlim(-0.5, df_res['NH4 calc'].dropna().index[-1] * 1.05)
    sns.despine(), fig.tight_layout(pad=0., rect=(0.015, 0.05, 0.98, 0.97))
    fig.canvas.draw()
    return fig


def plot_tanModel(df_res, fig1=None, ax1=None):
    ax1.cla()
    # preparation of figure plot
    if ax1 is None:
        fig1, ax1 = plt.subplots()
        ax1.set_aspect('auto')

    sns.despine()
    ax1.set_xlabel('Time / s'), ax1.set_ylabel('TAN / mg/L')

    ax1.plot(df_res['target_mg/L TAN'].dropna(), lw=1., ls=ls['target'], color='k', label='TAN target')
    ax1.plot(df_res['TAN calc'].dropna(), lw=1., color=dcolor['TAN'], label='TAN')

    ax1.set_xlim(-0.5, df_res['TAN calc'].dropna().index[-1] * 1.05)

    fig1.tight_layout(pad=.75, rect=(0.01, 0.005, 1, .95))
    fig1.canvas.draw()
    return fig1


# .....................................................................................................................
if __name__ == '__main__':
    app = QApplication(sys.argv)
    path = os.path.join('/Users/au652733/Python/Project_Fabi/', 'icon.png')
    app.setWindowIcon(QIcon(path))

    # where to display the GUI (which monitor in case there are several)
    view = MainWindow()

    # screen Size adjustment
    screen = app.primaryScreen()
    rect = screen.availableGeometry()
    view.setMaximumHeight(int(rect.height() * 0.9))
    view.setGeometry(50, 10, 1300, 750)

    sys.exit(app.exec_())
