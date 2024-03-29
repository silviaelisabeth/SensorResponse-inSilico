__author__ = 'szieger'
__project__ = 'in silico study for sensor response'

import basics_sensorSignal_v4 as bs
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns
import re
from PyQt5 import QtGui, QtCore, QtWidgets
from PyQt5.QtGui import *
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QMessageBox, QGridLayout,
                             QLabel, QLineEdit, QGroupBox, QFileDialog, QFrame, QDialog, QTableWidget, QTableWidgetItem,
                             QWizard, QWizardPage)
from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas,
                                                NavigationToolbar2QT as NavigationToolbar)
from scipy import integrate
import pandas as pd
import os

# .....................................................................................................................
# global parameter
dcolor = dict({'pH': '#1CC49A', 'sig pH': '#4B5258', 'para2': '#196E94', 'sig para2': '#314945', 'para3': '#DCA744',
               'sig para3': '#89621A', 'total para': '#A86349', 'background': '#383e42', 'font': 'white'})
ls = dict({'target': '-.', 'simulation': ':'})
save_type = ['png', 'svg']
fs_font = 10
sns.set_context('paper', font_scale=1.)

# fixed parameter depending on literature research
step_ph = 0.01
ph_deci = 2                # decimals for sensor sensitivity
ph_res = 1e-5              # resolution of the pH sensor

sig_max = 400              # maximal signal response at maximal acid concentration in mV
sig_bgd = .5               # background signal / offset in mV at 0M acid
E0 = 0.43                  # zero potential of the reference electrode
tsteps = 1e-3              # time steps for pH and base sensor (theory)

# electrochemical sensor · 2
sum_calib = 1              # concentration point for sum parameter at which the potential was measured
sigBase_max = 0.09         # maximal signal at maximal analyte concentration in mV (pH=1)
sigBase_bgd = 0.02         # background signal / offset in mV at 0M base
base_res = 1e-9            # resolution of the base sensor
sbgd_nhx = 0.03            # background signal

wizard_page_index = {"IntroPage": 0, "SimPage": 1}

# !!! TODO: finalize software to run it on all systems

# .....................................................................................................................
class MagicWizard(QWizard):
    def __init__(self):
        super(MagicWizard, self).__init__()
        self.introPage = IntroPage()
        self.setPage(wizard_page_index["IntroPage"], self.introPage)
        self.simPage = SimPage()
        self.setPage(wizard_page_index["SimPage"], self.simPage)

        # set start page
        self.setStartId(wizard_page_index["IntroPage"])

        # GUI layout
        self.setWindowTitle("In silico simulation of pH-dependent systems")
        # define Wizard style and certain options
        self.setWizardStyle(QWizard.ClassicStyle)
        self.setOptions(QtWidgets.QWizard.NoCancelButtonOnLastPage | QtWidgets.QWizard.IgnoreSubTitles)
        self.setStyleSheet("color: white; background-color: #383e42")

        # add a background image
        path = os.path.join('/Users/au652733/Python/Project_Fabi/', 'icon.png')
        pixmap = QtGui.QPixmap(path)
        pixmap = pixmap.scaled(200, 200, QtCore.Qt.KeepAspectRatio, QtCore.Qt.SmoothTransformation)
        self.setPixmap(QWizard.BackgroundPixmap, pixmap)


class IntroPage(QWizardPage):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setTitle("IntroPage | define parameter for initial plot")

        # create layout
        self.initUI()

        # connect checkbox and load file button with a function
        self.para1_conc_edit.editingFinished.connect(self.checkRange)
        self.para2_name_edit.editingFinished.connect(self.checkRange)
        self.load_button.clicked.connect(self.load_data)
        self.save_button.clicked.connect(self.save_path)

        # when all conditions are met, enable NEXT button:
        self.fname, self.fsave = QLineEdit(), QLineEdit()
        self.registerField("Data", self.fname)
        self.registerField("Storage path", self.fsave)
        self.registerField("parameter2*", self.para2_name_edit)
        self.registerField("conc. parameter1*", self.para1_conc_edit)
        self.registerField("conc. sum parameter*", self.paraSum_conc_edit)
        self.registerField("t90 parameter1*", self.para1_respT_edit)
        self.registerField("t90 parameter2*", self.para2_respT_edit)
        self.registerField("pKa parameter2*", self.para2_pKa_edit)
        self.registerField("plateau time", self.plateau_time_edit)

    def initUI(self):
        # path for measurement file (csv)
        self.load_button = QPushButton('Load meas. file', self)
        self.load_button.setFixedWidth(150), self.load_button.setFont(QFont('Helvetica Neue', fs_font))
        self.inputFileLineEdit = QLineEdit(self)
        self.inputFileLineEdit.setValidator(QtGui.QDoubleValidator())
        self.inputFileLineEdit.setFixedWidth(375), self.inputFileLineEdit.setAlignment(Qt.AlignRight)
        self.inputFileLineEdit.setFont(QFont('Helvetica Neue', int(0.9 * fs_font)))

        # directory to store files
        self.save_button = QPushButton('Storage path', self)
        self.save_button.setFixedWidth(150), self.save_button.setFont(QFont('Helvetica Neue', fs_font))
        self.inputSaveLineEdit = QLineEdit(self)
        self.inputSaveLineEdit.setValidator(QtGui.QDoubleValidator())
        self.inputSaveLineEdit.setFixedWidth(375), self.inputSaveLineEdit.setAlignment(Qt.AlignRight)
        self.inputSaveLineEdit.setFont(QFont('Helvetica Neue', int(0.9 * fs_font)))

        # create label for parameters | name, concentration, response time, and pKa (if applicable)
        # while the first parameter will always be pH, the second varies
        para1_lbl, para1_name_lbl, para2_lbl, para2_name_lbl = QLabel(self), QLabel(self), QLabel(self), QLabel(self)
        para1_lbl.setText('Sensor · 1'), para1_name_lbl.setText('Name')
        para2_lbl.setText('Sensor · 2'), para2_name_lbl.setText('Name')
        para1_name, self.para2_name_edit = QLabel(self), QLineEdit(self)
        para1_name.setText('pH'), para1_name.setAlignment(Qt.AlignRight)
        self.para2_name_edit.setAlignment(Qt.AlignRight)

        para1_conc_lbl, paraSum_conc_lbl = QLabel(self), QLabel(self)
        para1_conc_unit, paraSum_conc_unit = QLabel(self), QLabel(self)
        para1_conc_lbl.setText('pH value(s)'), paraSum_conc_lbl.setText('Concentration(s) sum parameter')
        para1_conc_unit.setText(''), paraSum_conc_unit.setText('mg/L')
        self.para1_conc_edit, self.paraSum_conc_edit = QLineEdit(self), QLineEdit(self)
        self.para1_conc_edit.setValidator(QRegExpValidator()), self.para1_conc_edit.setAlignment(Qt.AlignRight)
        self.paraSum_conc_edit.setValidator(QRegExpValidator()), self.paraSum_conc_edit.setAlignment(Qt.AlignRight)

        para1_respT_lbl, para2_respT_lbl = QLabel(self), QLabel(self)
        para1_respT_lbl.setText('Response time'), para2_respT_lbl.setText('Response time')
        para1_respT_unit, para2_respT_unit = QLabel(self), QLabel(self)
        para1_respT_unit.setText('s'), para2_respT_unit.setText('s')
        self.para1_respT_edit, self.para2_respT_edit = QLineEdit(self), QLineEdit(self)
        self.para1_respT_edit.setValidator(QRegExpValidator()), self.para1_respT_edit.setAlignment(Qt.AlignRight)
        self.para2_respT_edit.setValidator(QRegExpValidator()), self.para2_respT_edit.setAlignment(Qt.AlignRight)

        para2_pKa_lbl, self.para2_pKa_edit = QLabel(self), QLineEdit(self)
        para2_pKa_lbl.setText('pKa'), self.para2_pKa_edit.setValidator(QRegExpValidator())
        self.para2_pKa_edit.setAlignment(Qt.AlignRight)

        # create plateau time
        plateau_time_lbl, plateau_time_unit = QLabel(self), QLabel(self)
        plateau_time_lbl.setText('Plateau time'), plateau_time_unit.setText('s')
        self.plateau_time_edit = QLineEdit(self)
        self.plateau_time_edit.setValidator(QDoubleValidator()), self.plateau_time_edit.setAlignment(Qt.AlignRight)
        self.plateau_time_edit.setText('75.')

        divider, divider1 = QLabel(self), QLabel(self)
        divider.setText('   '), divider1.setText('   ')

        # -----------------------------------------------------------------------
        # creating main window (GUI)
        w = QWidget()
        # create layout grid
        mlayout = QVBoxLayout(w)
        vbox_top, vbox_middle, vbox_bottom = QHBoxLayout(), QHBoxLayout(), QHBoxLayout()
        mlayout.addLayout(vbox_top), mlayout.addLayout(vbox_middle), mlayout.addLayout(vbox_bottom)

        meas_settings = QGroupBox("Parameter Settings for Simulation")
        meas_settings.setStyleSheet("QGroupBox { font-weight: bold; } ")
        grid_load = QGridLayout()
        meas_settings.setFixedHeight(250), meas_settings.setFixedWidth(1000)
        meas_settings.setFont(QFont('Helvetica Neue', 12))
        vbox_top.addWidget(meas_settings)
        meas_settings.setLayout(grid_load)

        # include widgets in the layout
        # information for global system settings
        grid_load.addWidget(plateau_time_lbl, 0, 0)
        grid_load.addWidget(self.plateau_time_edit, 0, 2)
        grid_load.addWidget(plateau_time_unit, 0, 3)
        grid_load.addWidget(paraSum_conc_lbl, 1, 0)
        grid_load.addWidget(self.paraSum_conc_edit, 1, 2)
        grid_load.addWidget(paraSum_conc_unit, 1, 3)
        grid_load.addWidget(divider, 2, 2)

        # information for parameter 1
        grid_load.addWidget(para1_lbl, 3, 0)
        grid_load.addWidget(para1_name_lbl, 3, 1)
        grid_load.addWidget(para1_name, 3, 2)
        grid_load.addWidget(para1_conc_lbl, 4, 1)
        grid_load.addWidget(self.para1_conc_edit, 4, 2)
        grid_load.addWidget(para1_conc_unit, 4, 3)
        grid_load.addWidget(para1_respT_lbl, 5, 1)
        grid_load.addWidget(self.para1_respT_edit, 5, 2)
        grid_load.addWidget(para1_respT_unit, 5, 3)
        grid_load.addWidget(divider1, 6, 2)

        # information for parameter 2
        grid_load.addWidget(para2_lbl, 7, 0)
        grid_load.addWidget(para2_name_lbl, 7, 1)
        grid_load.addWidget(self.para2_name_edit, 7, 2)
        grid_load.addWidget(para2_respT_lbl, 8, 1)
        grid_load.addWidget(self.para2_respT_edit, 8, 2)
        grid_load.addWidget(para2_respT_unit, 8, 3)
        grid_load.addWidget(para2_pKa_lbl, 9, 1)
        grid_load.addWidget(self.para2_pKa_edit, 9, 2)

        divide_group = QGroupBox()
        divide_group.setStyleSheet("QGroupBox { border: 0px solid #383e42;}")
        grid_divide = QGridLayout()
        divide_group.setFixedHeight(100)
        vbox_middle.addWidget(divide_group)
        divide_group.setLayout(grid_divide)

        meas_file = QGroupBox("Define Directories")
        meas_file.setStyleSheet("QGroupBox { font-weight: bold; } ")
        grid_file = QGridLayout()
        meas_file.setFixedHeight(100), meas_file.setFixedWidth(1000)
        meas_file.setFont(QFont('Helvetica Neue', 12))
        vbox_bottom.addWidget(meas_file)
        meas_file.setLayout(grid_file)

        # include widgets in the layout
        grid_file.addWidget(self.load_button, 0, 0)
        grid_file.addWidget(self.inputFileLineEdit, 0, 1)
        grid_file.addWidget(self.save_button, 1, 0)
        grid_file.addWidget(self.inputSaveLineEdit, 1, 1)

        self.setLayout(mlayout)
        # margin order left, top, right, bottom
        self.layout().setContentsMargins(50, 50, 170, 170)
        self.layout().addStretch()

    def checkRange(self):
        if self.para1_conc_edit.text() and self.para2_name_edit.text():
            # check whether we work with H2S sensor / TDS project
            condA = ('H2S' in self.para2_name_edit.text() or 'HS-' in self.para2_name_edit.text())
            if condA or 'HS' in self.para2_name_edit.text():
                for ph in [float(i.strip()) for i in self.para1_conc_edit.text().split(',')]:
                    if ph > 9:
                        text = "Given pH value out of simulation range.  The H2S sensor simulation is only valid up " \
                               "to pH 9."
                        msgBox = error_message(text=text, type='Error')
                        returnValue = msgBox.exec()
                        if returnValue == QMessageBox.Ok:
                            pass

            elif 'TAN' in self.para2_name_edit.text() or 'TDS' in self.para2_name_edit.text():
                msgBox = error_message(text="Please enter the name of the sensor you want to simulate,  not a total "
                                            "parameter.", type='Error')
                returnValue = msgBox.exec()
                if returnValue == QMessageBox.Ok:
                    pass

    def load_data(self):
        # opens a dialog window in the current path
        fname, filter = QFileDialog.getOpenFileName(self, "Select specific txt file for temperature compensation",
                                                    "", "Text files (*.txt *.csv *xls)")

        if fname:
            self.inputFileLineEdit.setText(fname)
            df_general, df_ph, para2, df_para, df_total = bs.load_data(fname)

            # set parameter to run the simulation
            self.plateau_time_edit.setText(df_general.loc['plateau time (s)'].to_numpy()[0])

            # pH sensor
            self.para1_respT_edit.setText(df_ph.loc['t90 (s)', df_ph.columns[0]])
            ph_target = df_ph.loc['target values', df_ph.columns[0]]
            self.para1_conc_edit.setText(ph_target)

            # info for sensor 2 and total parameter concentration
            self.para2_name_edit.setText(para2)
            self.para2_respT_edit.setText(df_para.loc['t90 (s)', df_para.columns[0]])
            self.para2_pKa_edit.setText(df_para.loc['pKa', df_para.columns[0]])
            self.paraSum_conc_edit.setText(df_total.loc['target values'].to_numpy()[0])

    def save_path(self):
        fsave = QtWidgets.QFileDialog.getExistingDirectory(self, 'Select Folder')
        self.fsave.setText(fsave)
        if fsave:
            self.inputSaveLineEdit.setText(fsave)
            self.fsave.setText(fsave)

    def nextId(self) -> int:
        return wizard_page_index["SimPage"]


# .....................................................................................................................
class SimPage(QWizardPage):
    def __init__(self, parent=None):
        super(SimPage, self).__init__(parent)
        self.setTitle("in silico Simulation")

        # define certain parameter and gather from intro page
        self.df_res, self.dres = None, dict()

        # general layout
        self.initUI()

        # -----------------------------------------------------------
        # for all parameters - connect LineEdit with function
        self.initializePage()
        self.ph_t90_edit.returnPressed.connect(self.print_ph_t90)
        self.para2_t90_edit.returnPressed.connect(self.print_para2_t90)

        # connect buttons in navigation manager with functions
        self.plot_button.clicked.connect(self.check_parameter)
        self.int_button.clicked.connect(self.calc_integral)
        self._integral_counter = 0
        self.clearP_button.clicked.connect(self.clear_parameters)
        self.clearF_button.clicked.connect(self.clear_phsim)
        self.clearF_button.clicked.connect(self.clear_para2timedrive)
        self.clearF_button.clicked.connect(self.clear_sumtimdrive)
        self.save_button.clicked.connect(self.save)
        self.saveR_button.clicked.connect(self.save_report)

    def initUI(self):
        # parameters
        tsteady_lbl, tsteady_unit = QLabel(self), QLabel(self)
        tsteady_lbl.setText('Plateau time'), tsteady_unit.setText('s')
        self.tsteady_edit = QLineEdit(self)
        self.tsteady_edit.setValidator(QDoubleValidator()), self.tsteady_edit.setAlignment(Qt.AlignRight)

        # pH sensor settings
        pH_lbl = QLabel(self)
        pH_lbl.setText('pH value(s)')
        self.ph_edit = QLineEdit(self)
        self.ph_edit.setValidator(QRegExpValidator()), self.ph_edit.setAlignment(Qt.AlignRight)

        ph_t90_label, ph_t90_unit = QLabel(self), QLabel(self)
        ph_t90_label.setText('Response time'), ph_t90_unit.setText('s')
        self.ph_t90_edit = QLineEdit(self)
        self.ph_t90_edit.setValidator(QDoubleValidator()), self.ph_t90_edit.setAlignment(Qt.AlignRight)

        # second sensor settings
        self.para2_lbl, para2_unit = QLabel(self), QLabel(self)
        para2_unit.setText('mg/L')
        self.paraSum_conc_edit = QLineEdit(self)
        self.paraSum_conc_edit.setFixedWidth(100), self.paraSum_conc_edit.setValidator(QRegExpValidator())
        self.paraSum_conc_edit.setAlignment(Qt.AlignRight)

        para2_t90_lbl, para2_t90_unit = QLabel(self), QLabel(self)
        para2_t90_lbl.setText('Response time'), para2_t90_unit.setText('s')
        self.para2_t90_edit = QLineEdit(self)
        self.para2_t90_edit.setFixedWidth(100), self.para2_t90_edit.setValidator(QDoubleValidator())
        self.para2_t90_edit.setAlignment(Qt.AlignRight)

        para2_pka_lbl = QLabel(self)
        para2_pka_lbl.setText('pKa')
        self.para2_pka_edit = QLineEdit(self)
        self.para2_pka_edit.setFixedWidth(100), self.para2_pka_edit.setValidator(QDoubleValidator())
        self.para2_pka_edit.setAlignment(Qt.AlignRight)

        # general navigation
        self.plot_button = QPushButton('Plot', self)
        self.plot_button.setFixedWidth(100)
        self.int_button = QPushButton('Integral', self)
        self.int_button.setFixedWidth(100), self.int_button.setEnabled(False)
        self.clearP_button = QPushButton('Clear parameter', self)
        self.clearP_button.setFixedWidth(150)
        self.clearF_button = QPushButton('Clear plots', self)
        self.clearF_button.setFixedWidth(100)
        self.save_button = QPushButton('Save all', self)
        self.save_button.setFixedWidth(100)
        self.saveR_button = QPushButton('Save report', self)
        self.saveR_button.setFixedWidth(100)

        # ---------------------------------------------------------------------------------------
        w = QWidget()
        mlayout = QVBoxLayout(w)

        # 1st layer: box with vertical alignment to split in top and bottom
        vbox_top, vbox_middle, vbox_bottom = QHBoxLayout(), QHBoxLayout(), QHBoxLayout()
        mlayout.addLayout(vbox_top), mlayout.addLayout(vbox_middle), mlayout.addLayout(vbox_bottom)

        # 2nd layer: box with vertical alignment to split in left and right (left: parameters, m1: line, right: plots)
        hbox_left, hbox_m1, hbox_right = QVBoxLayout(), QVBoxLayout(), QVBoxLayout()
        vbox_top.addLayout(hbox_left), vbox_top.addLayout(hbox_m1), vbox_top.addLayout(hbox_right)

        # 3rd layer: left box with parameter settings (top, middle, bottom)
        hbox_ltop, hbox_lmiddle, hbox_lbottom = QHBoxLayout(), QHBoxLayout(), QHBoxLayout()
        hbox_ltop.setContentsMargins(5, 10, 10, 5)
        hbox_left.addLayout(hbox_ltop), hbox_left.addLayout(hbox_lmiddle), hbox_left.addLayout(hbox_lbottom)

        # 4th layer: right box split horizontally (top: individual sensors bottom: Sum)
        hbox_rtop, hbox_rbottom = QHBoxLayout(), QHBoxLayout()
        hbox_right.addLayout(hbox_rtop), hbox_right.addLayout(hbox_rbottom)

        # 5th layer: right top box split vertically (left: pH right: Base/Acid)
        hbox_tright, hbox_tleft = QVBoxLayout(), QVBoxLayout()
        hbox_rtop.addLayout(hbox_tright), hbox_rtop.addLayout(hbox_tleft)

        # draw additional "line" to separate parameters from plots and to separate navigation from rest
        vline = QFrame()
        vline.setFrameShape(QFrame.VLine | QFrame.Raised), vline.setLineWidth(2)
        hbox_m1.addWidget(vline)

        hline = QFrame()
        hline.setFrameShape(QFrame.HLine | QFrame.Raised), hline.setLineWidth(2)
        vbox_middle.addWidget(hline)

        # left side of main window (-> data treatment)
        vbox_top.addWidget(w)
        vbox_top.setContentsMargins(5, 10, 10, 5)

        # ---------------------------------------------------------------------------------------------------------
        # draw additional "line" to separate load from plot and plot from save
        vline1 = QFrame()
        vline1.setFrameShape(QFrame.VLine | QFrame.Raised), vline1.setLineWidth(2)
        vline2 = QFrame()
        vline2.setFrameShape(QFrame.VLine | QFrame.Raised), vline2.setLineWidth(2)
        vline3 = QFrame()
        vline3.setFrameShape(QFrame.VLine | QFrame.Raised), vline3.setLineWidth(2)

        # -------------------------------------------------------------------------------------------
        # GroupBoxes to structure the layout
        navigation_group = QGroupBox("Navigation Tool Bar")
        navigation_group.setStyleSheet("QGroupBox { font - weight: bold; border: 0px solid gray;}")
        grid_load = QGridLayout()
        navigation_group.setFixedHeight(70)
        vbox_bottom.addWidget(navigation_group)
        navigation_group.setLayout(grid_load)

        grid_load.addWidget(self.plot_button, 0, 1)
        grid_load.addWidget(self.int_button, 0, 2)
        grid_load.addWidget(vline2, 0, 3)
        grid_load.addWidget(self.clearP_button, 0, 4)
        grid_load.addWidget(self.clearF_button, 0, 5)
        grid_load.addWidget(vline3, 0, 6)
        grid_load.addWidget(self.save_button, 0, 7)
        grid_load.addWidget(self.saveR_button, 0, 8)

        # create GroupBox to structure the layout
        general_group = QGroupBox("General Settings")
        general_group.setStyleSheet("QGroupBox { font - weight: bold; border: 0px solid gray;}")
        general_group.setFixedHeight(75)
        grid_load = QGridLayout()
        grid_load.setSpacing(5), grid_load.setVerticalSpacing(2)

        # add GroupBox to layout and load buttons in GroupBox
        hbox_ltop.addWidget(general_group)
        general_group.setLayout(grid_load)
        grid_load.addWidget(tsteady_lbl, 1, 0)
        grid_load.addWidget(self.tsteady_edit, 1, 1)
        grid_load.addWidget(tsteady_unit, 1, 2)

        general_group.setContentsMargins(1, 15, 15, 1)
        hbox_ltop.addSpacing(10)

        # -----------------------
        # pH Sensor Settings
        phsens_group = QGroupBox("pH Sensor Settings")
        phsens_group.setStyleSheet("QGroupBox { font - weight: bold; border: 0px solid gray;}")
        phsens_group.setFixedHeight(110)
        grid_load = QGridLayout()
        grid_load.setSpacing(5), grid_load.setVerticalSpacing(1)

        # add GroupBox to layout and load buttons in GroupBox
        hbox_lmiddle.addWidget(phsens_group)
        phsens_group.setLayout(grid_load)

        grid_load.addWidget(pH_lbl, 0, 0)
        grid_load.addWidget(self.ph_edit, 0, 1)
        grid_load.addWidget(ph_t90_label, 1, 0)
        grid_load.addWidget(self.ph_t90_edit, 1, 1)
        grid_load.addWidget(ph_t90_unit, 1, 2)

        phsens_group.setContentsMargins(1, 15, 15, 1)
        hbox_lmiddle.addSpacing(10)

        # -----------------------
        # Settings of the other sensor
        self.para23sens_group = QGroupBox()
        self.para23sens_group.setStyleSheet("QGroupBox { font - weight: bold; border: 0px solid gray;}")
        self.para23sens_group.setFixedHeight(150)
        grid_load = QGridLayout()
        grid_load.setSpacing(5), grid_load.setVerticalSpacing(1)

        # add GroupBox to layout and load buttons in GroupBox
        hbox_lbottom.addWidget(self.para23sens_group)
        self.para23sens_group.setLayout(grid_load)

        grid_load.addWidget(self.para2_lbl, 0, 0)
        grid_load.addWidget(self.paraSum_conc_edit, 0, 1)
        grid_load.addWidget(para2_unit, 0, 2)
        grid_load.addWidget(para2_t90_lbl, 1, 0)
        grid_load.addWidget(self.para2_t90_edit, 1, 1)
        grid_load.addWidget(para2_t90_unit, 1, 2)
        grid_load.addWidget(para2_pka_lbl, 2, 0)
        grid_load.addWidget(self.para2_pka_edit, 2, 1)

        self.para23sens_group.setContentsMargins(1, 15, 15, 1)
        hbox_lbottom.addSpacing(10)

        # --------------------------------------------------------------------------------------------------------------
        # pH Simulation
        self.fig_phsim = plt.figure(facecolor=dcolor['background'])
        self.ax_phsim = self.fig_phsim.add_subplot()
        self.ax_phsim.set_facecolor(dcolor['background'])
        self.canvas_phsim = FigureCanvas(self.fig_phsim)
        self.ax_phsim.set_xlabel('Time / s', color=dcolor['font'])
        self.ax_phsim.set_ylabel('pH value', color=dcolor['font'])
        self.ax_phsim.tick_params(colors=dcolor['font'])
        [spine.set_edgecolor(dcolor['font']) for spine in self.ax_phsim.spines.values()]
        self.fig_phsim.tight_layout(pad=1, rect=(0.05, 0.1, 1, 0.98))
        sns.despine()

        # create GroupBox to structure the layout
        phsim_group = QGroupBox("pH Sensor Simulation")
        phsim_group.setMinimumWidth(450), phsim_group.setMinimumHeight(250)
        grid_phsim = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_tright.addWidget(phsim_group)
        phsim_group.setLayout(grid_phsim)
        grid_phsim.addWidget(self.canvas_phsim)

        # parameter 2-3 simulation
        self.fig_para2sim = plt.figure(facecolor=dcolor['background'])
        self.ax_para2sim = self.fig_para2sim.add_subplot()
        self.ax_para3sim = self.ax_para2sim.twinx()
        self.ax_para2sim.set_facecolor(dcolor['background'])
        self.ax_para2sim.set_xlabel('Time / s', color=dcolor['font'])
        self.canvas_para2sim = FigureCanvas(self.fig_para2sim)

        self.ax_para2sim.spines['top'].set_visible(False), self.ax_para3sim.spines['top'].set_visible(False)
        self.ax_para2sim.tick_params(colors=dcolor['font']), self.ax_para3sim.tick_params(colors=dcolor['font'])
        [spine.set_edgecolor(dcolor['font']) for spine in self.ax_para2sim.spines.values()]
        [spine.set_edgecolor(dcolor['font']) for spine in self.ax_para3sim.spines.values()]
        self.fig_para2sim.tight_layout(pad=1, rect=(0.05, 0.1, 0.95, 0.98))

        self.para2sim_group = QGroupBox("Single Parameter Simulation")
        self.para2sim_group.setMinimumWidth(450), self.para2sim_group.setMinimumHeight(250)
        grid_para2sim = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_tleft.addWidget(self.para2sim_group)
        self.para2sim_group.setLayout(grid_para2sim)
        grid_para2sim.addWidget(self.canvas_para2sim)

        # total parameter simulation
        self.fig_totalsim = plt.figure(facecolor=dcolor['background'])
        self.ax_totalsim = self.fig_totalsim.add_subplot()
        self.ax_totalsim.set_facecolor(dcolor['background'])
        self.canvas_totalsim = FigureCanvas(self.fig_totalsim)
        self.ntb = NavigationToolbar(self.canvas_totalsim, self, coordinates=False)
        self.ntb.setMinimumWidth(300)
        self.ntb.setStyleSheet("QToolBar { border: 0px, QToolButton{ background-color: white}}")

        self.ax_totalsim.set_xlabel('Time / s', color=dcolor['font'])
        self.ax_totalsim.tick_params(colors=dcolor['font'])
        [spine.set_edgecolor(dcolor['font']) for spine in self.ax_totalsim.spines.values()]

        self.fig_totalsim.tight_layout(pad=.75, rect=(0, 0.1, 1, .95))
        sns.despine()

        self.totalsim_group = QGroupBox()
        self.totalsim_group.setMinimumWidth(900), self.totalsim_group.setMinimumHeight(300)
        grid_totalsim = QGridLayout()

        # add GroupBox to layout and load buttons in GroupBox
        hbox_rbottom.addWidget(self.totalsim_group)
        self.totalsim_group.setLayout(grid_totalsim)
        grid_totalsim.addWidget(self.canvas_totalsim)
        grid_totalsim.addWidget(self.ntb)

        # -------------------------------------------------------------------------------------------------------------
        self.setLayout(mlayout)

    # ---------------------------------------------------
    def print_ph_t90(self):
        print('pH sensor response: ', self.para1_t90_edit.text(), 's')

        # in case set response time equals 0, return warning
        if float(self.para1_t90_edit.text()) == 0:
            msgBox = error_message(text="Please enter a valid sensor response time for the pH sensor, notably a "
                                        "number > 0")
            returnValue = msgBox.exec()
            if returnValue == QMessageBox.Ok:
                pass

    def print_para2_t90(self):
        print('{} sensor response: '.format(self.para2_name_edit.text()), self.para2_t90_edit.text(), 's')
        # in case set response time equals 0, return warning
        if float(self.para2_t90_edit.text()) == 0:
            msgBox = error_message(text="Please enter a valid sensor response time for the total sensor, notably a "
                                        "number > 0")
            returnValue = msgBox.exec()
            if returnValue == QMessageBox.Ok:
                pass

    def initializePage(self):
        self.save_path = self.field("Storage path")
        self.para1_conc, self.para2_conc = self.field("conc. parameter1"), self.field("conc. sum parameter")
        if isinstance(self.field("t90 parameter1"), type(None)):
            self.para1_respT = self.field("t90 parameter1")
        else:
            self.para1_respT = str(self.field("t90 parameter1").replace(',', '.'))
        if isinstance(self.field("t90 parameter2"), type(None)):
            self.para2_respT = self.field("t90 parameter2")
        else:
            self.para2_respT = str(self.field("t90 parameter2").replace(',', '.'))
        if isinstance(self.field("pKa parameter2"), type(None)):
            self.para2_pka = self.field("pKa parameter2")
        else:
            self.para2_pka = str(self.field("pKa parameter2").replace(',', '.'))
        if isinstance(self.field("plateau time"), type(None)):
            self.plateau_time = self.field("plateau time")
        else:
            self.plateau_time = float(self.field("plateau time").replace(',', '.'))

        self.para2_name = self.field("parameter2")
        # identify the system according to given parameter
        if self.para2_name:
            if 'NH' in self.para2_name:
                self.total_para, self.total_para_grp = '$TAN$', 'TAN'
                self.para2, self.para2_grp, self.para3, self.para3_grp = '$NH_3$', 'NH3', '$NH_4^+$', 'NH4+'
                # define base as the analyte if it is not stated differently
                if '3' in self.para2_name or 'TAN' in self.para2_name:
                    self.analyte = 'NH3'
                else:
                    self.analyte = 'NH4+'
            elif 'S' in self.para2_name:
                self.total_para, self.total_para_grp = '$TDS$', 'TDS'
                self.para2, self.para2_grp, self.para3, self.para3_grp = '$H_2S$', 'H2S', '$HS^-$', 'HS-'
                # define base as the analyte if it is not stated differently
                if '2' in self.para2_name:
                    self.analyte = 'H2S'
                else:
                    self.analyte = 'HS-'
        else:
            self.para2_grp, self.para3_grp, self.total_para_grp = 'Base', 'Acid', 'Sum Parmeter'
            self.para2, self.para3, self.total_para = 'base', 'acid', 'sum parmeter'
            self.analyte = 'sum parmeter'

        # fill in parameters
        if self.plateau_time:
            self.tsteady_edit.setText(str(self.plateau_time))
        self.ph_edit.setText(self.para1_conc)
        self.ph_t90_edit.setText(self.para1_respT)

        self.paraSum_conc_edit.setText(self.para2_conc)
        self.para2_lbl.setText('c({})'.format(self.total_para_grp))
        self.para2_t90_edit.setText(self.para2_respT)
        self.para2_pka_edit.setText(self.para2_pka)

        self.ax_para2sim.set_ylabel('{} / mg/L'.format(self.para2), color=dcolor['para2'])
        self.ax_para3sim.set_ylabel('{} / mg/L'.format(self.para3), color=dcolor['para3'])
        self.ax_totalsim.set_ylabel('{} / mg/L'.format(self.total_para), color=dcolor['font'])

        self.totalsim_group.setTitle("{} Simulation".format(self.total_para_grp))
        self.para2sim_group.setTitle("{}/{} Simulation".format(self.para2_grp, self.para3_grp))
        self.para23sens_group.setTitle("{}/{} Sensor Settings".format(self.para2_grp, self.para3_grp))

        self.fig_para2sim.canvas.draw(), self.fig_totalsim.canvas.draw()

        return

    # ---------------------------------------------------
    def clear_parameters(self):
        # re-write default parameters
        self.tsteady_edit.setText('75.')
        self.ph_edit.setText('8.4, 10.9')
        self.ph_t90_edit.setText('30.')
        self.paraSum_conc_edit.setText('100.0')
        self.para2_t90_edit.setText('60.')
        self.para2_pka_edit.setText('9.25')

    def clear_phsim(self):
        self.ax_phsim.cla()
        self.ax_phsim.set_xlabel('Time / s', color=dcolor['font'])
        self.ax_phsim.set_ylabel('pH value', color=dcolor['font'])
        self.ax_phsim.tick_params(colors=dcolor['font'])
        self.fig_phsim.tight_layout(pad=0, rect=(0.015, 0.05, 0.98, 0.97))
        sns.despine()
        self.fig_phsim.canvas.draw()

    def clear_para2timedrive(self):
        self.ax_para2sim.cla(), self.ax_para3sim.cla()
        self.ax_para2sim.set_xlabel('Time / s', color=dcolor['font'])
        self.ax_para2sim.set_ylabel('{} / mg/L'.format(self.para2), color=dcolor['para2'])
        self.ax_para3sim.set_ylabel('{} / mg/L'.format(self.para3), color=dcolor['para3'])

        self.fig_para2sim.tight_layout(pad=0., rect=(0.015, 0.05, 0.98, 0.97))
        self.ax_para2sim.spines['top'].set_visible(False), self.ax_para3sim.spines['top'].set_visible(False)
        self.ax_para2sim.tick_params(colors=dcolor['font']), self.ax_para3sim.tick_params(colors=dcolor['font'])
        self.fig_para2sim.canvas.draw()

    def clear_sumtimdrive(self):
        self.ax_totalsim.cla()
        self.ax_totalsim.set_xlabel('Time / s', color=dcolor['font'])
        self.ax_totalsim.set_ylabel('{} / mg/L'.format(self.total_para), color=dcolor['font'])

        self.fig_totalsim.tight_layout(pad=.75, rect=(0.01, 0.01, 1, .95))
        self.ax_totalsim.tick_params(colors=dcolor['font'])
        sns.despine()
        self.fig_totalsim.canvas.draw()

    def align_concentrations(self, ls_ph, ls_cSum_gmL):
        if len(ls_cSum_gmL) == len(ls_ph):
            pass
        elif len(ls_cSum_gmL) == 1 and len(ls_ph) > 1:
            ls_cSum_gmL = ls_cSum_gmL * len(ls_ph)
        elif len(ls_cSum_gmL) > 1 and len(ls_ph) == 1:
            ls_ph = ls_ph * len(ls_cSum_gmL)
        else:
            if max(len(ls_ph), len(ls_cSum_gmL)) % min(len(ls_ph), len(ls_cSum_gmL)) == 0:
                print('double the length of each other')
                if len(ls_ph) < len(ls_cSum_gmL):
                    ls_ph = ls_ph * int(len(ls_cSum_gmL) / 2)
                else:
                    ls_cSum_gmL = ls_cSum_gmL * int(len(ls_ph) / 2)
            else:
                if len(ls_ph) < len(ls_cSum_gmL):
                    ls_ph.append(ls_ph[-1])
                else:
                    ls_cSum_gmL.append(ls_cSum_gmL[-1])
        return ls_ph, ls_cSum_gmL

    def parameter_prep(self, t_plateau, ls_ph, ls_total, pKa):
        # target fluctuation of analytes
        target_ph = bs.target_fluctuation(ls_conc=ls_ph, tstart=0, tstop=t_plateau * 2, nP=1, analyte='pH')
        target_sum = bs.target_fluctuation(ls_conc=ls_total, tstart=0, tstop=t_plateau * 2, nP=1, analyte='Sum')

        # target concentration para2/para3 based on target pH (para2... base, para3... acid)
        para2 = pd.DataFrame(bs.henderson_Sum4Base(pKa=pKa, pH=target_ph['signal pH'].to_numpy(),
                                                   c_sum=target_sum['signal Sum'].to_numpy()), index=target_ph.index,
                             columns=['signal Base'])
        para3 = pd.DataFrame(bs.henderson_Sum4Acid(pKa=pKa, pH=target_ph['signal pH'].to_numpy(),
                                                   c_sum=target_sum['signal Sum'].to_numpy()), index=target_ph.index,
                             columns=['signal Acid'])

        # calculate target concentration of individual NHx parameters
        df_res = pd.concat([target_ph, target_sum, para2, para3], axis=1)
        return df_res

    # ---------------------------------------------------
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
                    msgBox = error_message(text="Simulate before saving. Simulation might have been unsuccessful.")
                    returnValue = msgBox.exec()
                    if returnValue == QMessageBox.Ok:
                        pass
                else:
                    # save all reports
                    self.save2txt(fname_save)

                    # save figures in separate folder
                    for f in self.dic_figures.keys():
                        figure_name = fname_save.split('.')[0] + '_Graph-' + '-'.join(f.split(' ')) + '.'
                        for t in save_type:
                            self.dic_figures[f].savefig(figure_name + t, dpi=300, transparent=True)

            except NameError:
                msgBox = error_message(text="Simulate before saving", type='Error')
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
                    msgBox = error_message(text="Simulate before saving", type='Error')
                    returnValue = msgBox.exec()
                    if returnValue == QMessageBox.Ok:
                        pass
                else:
                    # save the reports
                    self.save2txt(fname_save)

            except NameError:
                msgBox = error_message(text="Simulate before saving", type='Error')
                returnValue = msgBox.exec()
                if returnValue == QMessageBox.Ok:
                    pass

    def save2txt(self, fname_save):
        # if integral available, convert to DataFrame
        if len(self.dres.keys()) != 0:
            int_out = pd.DataFrame.from_dict(self.dres)
            int_out.index = ['integration range (s)', 'target pH', 'target Sum', 'integral target Sum',
                             'integral observed Sum', 'error']
        else:
            int_out = None

        # save output
        output = bs.save_report(para_meas=self.para_meas, sensor_ph=self.sensor_ph, df_res=self.df_res, dres=self.dres,
                                sensor_para2=self.sensor_para2, total_lbl=self.total_para_grp)

        with open(fname_save, 'w') as f:
            output[0].to_csv(f, sep='\t', decimal='.')
        with open(fname_save, 'a') as f:
            output[1].to_csv(f, sep='\t', decimal='.', header=False)
        if isinstance(int_out, type(None)):
            pass
        else:
            with open(fname_save.split('.')[0] + '_integral.txt', 'w') as fint:
                int_out.to_csv(fint, sep='\t', decimal='.', header=False)

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
        # get sum concentration(s) and pH value(s) and align list of fluctuation points
        ls_cSum_gmL = self._linEdit2list(line=self.paraSum_conc_edit.text())
        ls_ph = self._linEdit2list(line=self.ph_edit.text())
        ls_ph, ls_cSum_gmL = self.align_concentrations(ls_ph=ls_ph, ls_cSum_gmL=ls_cSum_gmL)

        # determine all required target concentrations
        df_target = self.parameter_prep(ls_ph=ls_ph, ls_total=ls_cSum_gmL, pKa=float(self.para2_pka_edit.text()),
                                        t_plateau=float(self.tsteady_edit.text()))

        # calculate target concentrations of sensor 2
        if self.analyte in ['HS-', 'hs-', 'HS', 'NH3', 'nh3']:
            df_target_sens2 = pd.DataFrame(df_target['signal Base'])
        elif self.analyte in ['H2S', 'h2s', 'NH4+', 'NH4', 'nh4', 'nh4+']:
            df_target_sens2 = pd.DataFrame(df_target['signal Acid'])
        else:
            df_target_sens2 = pd.DataFrame(df_target['signal Base'])
        ls_target_sens2 = list(df_target_sens2[df_target_sens2.columns[0]].drop_duplicates().values)

        # collect all relevant parameter
        self.sensor_ph = dict({'E0': E0, 't90': float(self.ph_t90_edit.text()), 'resolution': ph_res,
                               'start pH': ls_ph[0], 'time steps': tsteps, 'sensitivity': ph_deci,
                               'pH target': df_target['signal pH'], 'set values': ls_ph})

        self.sensor_para2 = dict({'sensitivity': ph_deci, 'analyte': self.analyte, 'background signal': sbgd_nhx,
                                  'pKa': float(self.para2_pka_edit.text()),  'signal min': sigBase_bgd,
                                  'signal max': sigBase_max, 'resolution': base_res, 'time steps': tsteps,
                                  't90': float(self.para2_t90_edit.text()), 'Base target': df_target['signal Base'],
                                  'Acid target': df_target['signal Acid'], 'Sum target': df_target['signal Sum'],
                                  'calib Total': sum_calib, 'set values': ls_cSum_gmL, 'set sensor2': ls_target_sens2})
        self.para_meas = dict({'plateau time': float(self.tsteady_edit.text())})

        # -----------------------------------------------------------------------------------------------------------
        # check pH sensor
        try:
            if float(self.sensor_ph['t90']) == 0:
                msgBox = error_message(text="Please enter a valid sensor response time for the pH sensor, in "
                                            "particular a number > 0", type='Error')
                msgBox.exec_()
            else:
                try:
                    # check Base/Acid sensor
                    if float(self.sensor_para2['t90']) == 0:
                        msgBox = error_message(text="Please enter a valid sensor response time for the pH sensor, in "
                                                    "particular a number > 0", type='Error')
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
        self.clear_phsim(), self.clear_para2timedrive(), self.clear_sumtimdrive()
        self.ls_xcoords = []

        # target concentrations
        df_res = pd.concat([self.sensor_ph['pH target'], self.sensor_para2['Base target'],
                            self.sensor_para2['Acid target'], self.sensor_para2['Sum target']], axis=1)
        df_res.columns = ['target pH', 'target_mg/L Base', 'target_mg/L Acid', 'target_mg/L Sum']

        # individual sensor - reduce by individual plotting (different function)
        [df_res, cplateaupH, cplateauTotal,
         para_para] = bs._alignSensorSettings(sensor_ph=self.sensor_ph, sensor_base=self.sensor_para2, df_target=df_res,
                                              para_meas=self.para_meas)
        df_res.index = [round(u, 2) for u in df_res.index]
        df_res = df_res.dropna()

        # pH sensor response
        [df_pHrec, df_pHcalc] = bs.pH_sensor(cplateau=cplateaupH, sensor_ph=self.sensor_ph, para_meas=self.para_meas)
        df_pHrec.index = [round(i, 2) for i in df_pHrec.index]

        # include sensor response and drift to result dataframe
        df_res = pd.concat([df_res, df_pHrec.loc[df_res.index], df_pHcalc.loc[df_res.index]], axis=1).sort_index(axis=1)

        # para2 sensor - calculate individual sensor response as well as total parameter
        df_res = bs.para2Sensorcalc(df_res=df_res, analyte=self.sensor_para2['analyte'], cplateauSum=cplateauTotal,
                                    df_pHcalc=df_pHcalc, sensor_Base=self.sensor_para2, para_meas=self.para_meas,
                                    sensor_ph=self.sensor_ph)
        self.df_res = df_res.sort_index(axis=1)

        # --------------------------------------------------------------------------------------------------------------
        # plotting part | individual sensor - target vs record
        _ = plot_phsensor(df_res=self.df_res, color=dcolor['font'], fig=self.fig_phsim, ax=self.ax_phsim)
        _ = plot_parasensor(df_res=self.df_res, para_base=self.para2, para_acid=self.para3, fig=self.fig_para2sim,
                            ax=self.ax_para2sim, ax1=self.ax_para3sim, color=dcolor['font'])

        # final total parameter model
        _ = plot_totalModel(df_res=self.df_res, para_sum=self.total_para, color=dcolor['font'], fig1=self.fig_totalsim,
                            ax1=self.ax_totalsim)

        # re-draw figures for saving with transparent/white background
        fig_total = plot_totalModel(df_res=self.df_res, para_sum=self.total_para, color='k')
        fig_pH = plot_phsensor(df_res=self.df_res, color='k')
        fig_base = plot_parasensor(df_res=self.df_res, para_base=self.para2, para_acid=self.para3, color='k')

        # allow click events for total parameter for integration (simpson for discrete measurement data)
        self.fig_totalsim.canvas.mpl_connect('button_press_event', self.onclick_integral)

        # --------------------------------------------------------------------------------------------------------------
        # collect for result output (save data)
        self.dic_figures = dict({'pH': fig_pH, 'Base': fig_base, 'Total': fig_total})

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
        self.ax_totalsim.vlines(x=self.ls_xcoords, ymin=self.df_res['Sum calc'].min(),
                                ymax=self.df_res['Sum calc'].max(), color='k', lw=0.15)
        if len(self.ls_xcoords) == 2:
            self.ax_totalsim.axvspan(self.ls_xcoords[0], self.ls_xcoords[1], color='grey', alpha=0.3)

        # update figure plot
        self.fig_totalsim.canvas.draw()

        # allow integral calculation as soon as xcoords have been collected
        self.int_button.setEnabled(True)

    def calc_integral(self):
        self._integral_counter += 1
        try:
            self.ls_xcoords
        except:
            self.ls_xcoords = list()

        # calculate integral for calculated/target total parameter using simpson's rule for discrete data
        # subtract calculated Sum by target Sum
        if len(self.ls_xcoords) == 2:
            self.df4int = self.df_res[['Sum calc', 'target_mg/L Sum']].loc[min(self.ls_xcoords):max(self.ls_xcoords)]
            dfint_calc = integrate.simps(self.df4int['Sum calc'].to_numpy(), x=self.df4int.index)
            dfint_target = integrate.simps(self.df4int['target_mg/L Sum'].to_numpy(), x=self.df4int.index)

            # Integration range (s), target pH, target Sum, integral target Sum, integral Sum, error
            result = list([(round(min(self.ls_xcoords), 2), round(max(self.ls_xcoords)), 2),
                           self.df_res['target pH'].loc[min(self.ls_xcoords):max(self.ls_xcoords)].mean(),
                           self.df_res['target_mg/L Sum'].mean(), dfint_calc, dfint_target, dfint_calc - dfint_target])

            # add current results to dictionary (where all selected peak information are stored)
            self.dres[self._integral_counter] = result

            # open a pop up window with options to select what shall be saved
            global wInt
            wInt = IntegralWindow(self.dres, self.para_meas, self.sensor_ph, self.sensor_para2, self.total_para_grp,
                                  self._integral_counter)
            if wInt.isVisible() is False:
                wInt.show()


class IntegralWindow(QDialog):
    def __init__(self, res_int, para_meas, sensor_ph, sensor_para2, sum_lbl, count):
        super().__init__()
        self.result, self.count, self.sum_lbl = res_int, count, sum_lbl
        self.para_meas, self.sensor_ph, self.sensor_para2 = para_meas, sensor_ph, sensor_para2
        self.initUI()

        # when checkbox selected, save information in registered field
        self.reset_button.clicked.connect(self.reset)
        self.save_button.clicked.connect(self.save_integral)
        self.close_button.clicked.connect(self.close_window)

        # execute function as soon as the window pups up or as soon as integral_button is clicked
        self.res2table()

    def initUI(self):
        self.setWindowTitle("error calculation via Simpson's integration")
        # define geometry of integration window x, y, width, height
        self.setGeometry(50, 100, 700, 300)

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
        self.tab_report.setHorizontalHeaderLabels(['Integration range (s)', 'Target pH', 'Target ' + self.sum_lbl,
                                                   'Integral target ' + self.sum_lbl, 'Integral observed ' + self.sum_lbl,
                                                   'Error'])
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
            # columns: Integration range (s), target pH, target Sum, integral target Sum, integral Sum, error
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
        self.tab_report.setHorizontalHeaderLabels(['integration range (s)', 'target pH', 'target Sum',
                                                   'integral target Sum', 'integral observed Sum', 'error'])
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
                    msgBox = error_message(text="Simulate before saving. Simulation might have been unsuccessful.")
                    returnValue = msgBox.exec()
                    if returnValue == QMessageBox.Ok:
                        pass
                else:
                    # save output now
                    output = bs.save_integral(para_meas=self.para_meas, sensor_ph=self.sensor_ph, dres=self.result,
                                              sensor_para2=self.sensor_para2, total_lbl=self.sum_lbl)

                    with open(fname_save, 'w') as f:
                        output.to_csv(f, sep='\t', decimal='.', header=False)

            except NameError:
                msgBox = error_message(text="Simulate before saving.")
                returnValue = msgBox.exec()
                if returnValue == QMessageBox.Ok:
                    pass

    def close_window(self):
        self.hide()


# .....................................................................................................................
def plot_phsensor(df_res, color, fig=None, ax=None):
    # preparation of figure plot
    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect('auto')
    else:
        ax.cla()
    ax.set_xlabel('Time / s'), ax.set_ylabel('pH value')

    ax.plot(df_res['target pH'].dropna(), ls='-.', lw=1., color=color)
    ax.plot(df_res['pH calc'].dropna(), color=dcolor['pH'])

    ax.set_xlim(-0.5, df_res['pH calc'].dropna().index[-1]*1.05)
    sns.despine(), fig.tight_layout(pad=0.6, rect=(0., 0.015, 1, 0.98))
    fig.canvas.draw()
    return fig


def plot_parasensor(df_res, para_base, para_acid, color, fig=None, ax=None, ax1=None):
    # preparation of figure plot
    if ax is None:
        fig, ax = plt.subplots()
        ax.set_aspect('auto')
        ax1 = ax.twinx()
        ax.spines['top'].set_visible(False), ax.spines['right'].set_visible(False), ax.tick_params(colors=color)
        ax1.spines['top'].set_visible(False), ax1.spines['right'].set_visible(True), ax1.tick_params(colors=color)
    else:
        ax.cla(), ax1.cla()
    ax.set_xlabel('Time / s', color=color)
    ax.set_ylabel('{} / mg/L'.format(para_acid), color=dcolor['para3'])
    ax1.set_ylabel('{} / mg/L'.format(para_base), color=dcolor['para2'])

    ax1.plot(df_res['target_mg/L Base'].dropna(), lw=1., ls=':', color=color, label=para_base + ' target')
    ax1.plot(df_res['Base calc'].dropna(), lw=1., color=dcolor['para2'], label=para_base)

    ax.plot(df_res['target_mg/L Acid'].dropna(), lw=1., ls='--', color=color, label=para_acid + ' target')
    ax.plot(df_res['Acid calc'].dropna(), lw=1., color=dcolor['para3'], label=para_acid)

    ax.set_xlim(-0.5, df_res['Acid calc'].dropna().index[-1] * 1.05)

    fig.tight_layout(pad=0., rect=(0.015, 0.05, 0.98, 0.97))
    fig.canvas.draw()
    return fig


def plot_totalModel(df_res, para_sum, color, fig1=None, ax1=None):
    # preparation of figure plot
    if ax1 is None:
        fig1, ax1 = plt.subplots()
        ax1.set_aspect('auto')
    else:
        ax1.cla()

    ax1.set_xlabel('Time / s', color=color), ax1.set_ylabel('{} / mg/L'.format(para_sum))
    ax1.tick_params(colors=color)

    ax1.plot(df_res['target_mg/L Sum'].dropna(), lw=1, ls=ls['target'], color=color, label=para_sum + ' target')
    ax1.plot(df_res['Sum calc'].dropna(), lw=1., color=dcolor['total para'], label=para_sum)

    ax1.set_xlim(-0.5, df_res['Sum calc'].dropna().index[-1] * 1.05)
    sns.despine(), fig1.tight_layout(pad=.75, rect=(0.01, 0.005, 1, .95))
    fig1.canvas.draw()
    return fig1


# .....................................................................................................................
def error_message(text, type='Warning'):
    msgBox = QMessageBox()
    msgBox.setIcon(QMessageBox.Information)
    msgBox.setText(text)
    msgBox.setFont(QFont('Helvetica Neue', 11))
    msgBox.setWindowTitle(type)
    msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)

    return msgBox


# .....................................................................................................................
if __name__ == '__main__':
    import sys

    app = QtWidgets.QApplication(sys.argv)
    path = os.path.join('/Users/au652733/Python/Project_Fabi/', 'icon.png')
    app.setWindowIcon(QIcon(path))

    # options available: 'Breeze', 'Oxygen', 'QtCurve', 'Windows', 'Fusion'
    app.setStyle('Fusion')

    # screen Size adjustment
    Wizard = MagicWizard()
    screen = app.primaryScreen()
    rect = screen.availableGeometry()
    Wizard.setMaximumHeight(int(rect.height() * 0.9))
    # position (screen width, screen height), width, height
    Wizard.setGeometry(int(rect.width() * 0.075), int(rect.height() * 0.1), 710, 550)

    # show wizard
    Wizard.show()
    sys.exit(app.exec_())
