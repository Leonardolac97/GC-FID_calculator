# ██████████████████████████████████████████████████████████████████████
# █                         GC-FID calculator                          █
# █                                                                    █
# █  Debug AI-assisted: OpenAI ChatGPT (GPT-5.3)                       █
# █  04.2026                                                           █
# ██████████████████████████████████████████████████████████████████████


# -----------------------------
# Libraries
# -----------------------------

import sys
import os
import re
import pandas as pd
import fitz
import numpy as np
import periodictable as pt
import seaborn as sns
import json
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from matplotlib.colors import (LinearSegmentedColormap, TwoSlopeNorm)


from PyQt6.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout,
    QListWidget, QPushButton, QLineEdit, QFileDialog,
    QLabel, QSpinBox, QDoubleSpinBox, QTabWidget,
    QComboBox, QTableWidget, QHeaderView, QTableWidgetItem, QMessageBox
)
from PyQt6.QtCore import Qt, QT_VERSION_STR

from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar
)
from matplotlib.figure import Figure

# -----------------------------------------------------------------------------------------------------------------
# Variables preamble: these variables can be edited depending on your reaction:
# - self.retention_times_df (possible to upload new file)
# - self.key_products (Hardcoded - to be changed directly in the script)
# - self.groups are the only (Hardcoded - to be changed directly in the script)
# ------------------------------------------------------------------------------------------------------------------

class RetentionTimeFilterApp(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Retention Time Filter")
        self.setGeometry(200, 50, 2000, 1100)

        # Default values
        self.retention_times_df = pd.DataFrame([
            [9.386, "Formaldehyde"],
            [9.520, "Methanol"],
            [9.764, "Ethanol"],
            [10.989, "Cyclohexane"],
            [12.633, "Toluene"],
            [14.977, "Cyclohexanol"],
            [16.197, "Anisole"], #[16.250, "Anisole"],
            [17.664, "Phenol"],
            [18.648, "o-Methylanisole"],
            [18.923, "m-Methylanisole"],
            [19.689, "o-Cresol"],
            [20.062, "Decalin trans"],
            [20.736, "Guaiacol"],
            [21.208, "Decalin cis"],
            [22.072, "2,5-Dimethylphenol"],
            [22.795, "2,3-Dimethylphenol"],
            [22.952, "Creosol"],
            [23.148, "3,4-Dimethylphenol"],
            [24.191, "2,3,6-Trimethylphenol"],
            [24.990, "2,4,6-Trimethylphenol"],
            [25.469, "Pentamethylbenzene"],
            [25.722, "3,4,5-Trimethylphenol"],
            [27.064, "2,3,4,6-Tetramethylphenol"],
            [29.098, "Hexamethylbenzene"],
            [30.608, "2,3,4,5,6-Pentamethylphenol"]
        ], columns=["Retention Time", "Compound"])

        # ----------------------------------------
        # Default reaction settings
        # These values act as fallback guesses and
        # can be overwritten by uploaded reference files
        # ----------------------------------------

        self.tolerance = 0.09
        self.calib_flow = 0.06
        self.experiment_flow = 0.06  # g/min

        self.slope_calibration = 1064289.148
        self.linear_co_calibration = 136790.506

        self.reactant_name = "Guaiacol"
        self.reactant_rt = 20.736
        self.reactant_formula = "C7H8O2"
        self.reactant_molar_mass = 124.139
        self.reactant_nC = 7
        self.reactant_wt = 9.992698759

        self.solvent_name = "Decalin"

        reactant_fraction = self.reactant_wt / 100
        remaining = 1 - reactant_fraction

        self.solvents = [
            {"name": "Decalin trans", "rt": 20.062, "molar_mass": 138.23, "wt_fraction": remaining * 0.40},
            {"name": "Decalin cis", "rt": 21.208, "molar_mass": 138.23, "wt_fraction": remaining * 0.60}
        ]
        self.key_products = [
            "Phenol",
            "o-Cresol",
            "2,5-Dimethylphenol",
            "Anisole",
            "Methanol",
            "Formaldehyde"
        ]

        self.groups = {
            "Other methylphenols": [
                "2,3,4,5,6-Pentamethylphenol",
                "2,3,4,6-Tetramethylphenol",
                "2,3,6-Trimethylphenol",
                "2,3-Dimethylphenol",
                "2,4,6-Trimethylphenol",
                "3,4,5-Trimethylphenol",
                "3,4-Dimethylphenol"
            ],
            "Minor compounds": [
                "Ethanol",
                "Cyclohexanol",
                "o-Methylanisole",
                "m-Methylanisole"
            ],
            "Fully deoxygenated": [
                "Toluene",
                "Pentamethylbenzene",
                "Hexamethylbenzene",
                "Cyclohexane"
            ]
        }

        self.pdf_data = pd.DataFrame()
        self.filtered_data = pd.DataFrame()
        self.csv_folder = ""
        self.plot_folder = ""

        self.init_ui()

    def init_ui(self):
        main_layout = QVBoxLayout(self)

        # Create tab widget
        self.tabs = QTabWidget()
        main_layout.addWidget(self.tabs)

        # --- Tab 1: Retention Time Filter ---
        self.main_tab = QWidget()
        self.tabs.addTab(self.main_tab, "Retention time filter")
        self.setup_main_tab()

        # --- Tab 2: Product Trend ---
        self.product_trend_tab = QWidget()
        self.tabs.addTab(self.product_trend_tab, "Product trend")
        self.setup_product_trend_tab()

        # --- Tab 3: Product Calibration ---
        self.product_calibration_tab = QWidget()
        self.tabs.addTab(self.product_calibration_tab, "Product calibration")
        self.setup_product_calibration_tab()

        # --- Tab 4: Grouped Mol% Plot ---
        self.molpercent_tab = QWidget()
        self.tabs.addTab(self.molpercent_tab, "Mol percentage over time")
        self.setup_molpercent_tab()

        # --- Tab 5: Product Bar Chart ---
        #self.hourly_product_tab = QWidget()
        #self.tabs.addTab(self.hourly_product_tab, "Products bar chart")
        #self.setup_hourly_product_tab()

        # --- Tab 6: Selectivity ---
        self.selectivity_tab = QWidget()
        self.tabs.addTab(self.selectivity_tab, "Selectivity")
        self.setup_selectivity_tab()

        # --- Tab 7: Pearson correlation ---
        self.pearson_tab = QWidget()
        self.tabs.addTab(self.pearson_tab, "Pearson correlation chart")
        self.setup_pearson_tab()

        # --- Tab 8: Carbon Balance ---
        self.carbon_balance_tab = QWidget()
        self.tabs.addTab(self.carbon_balance_tab, "Carbon balance")
        self.setup_carbon_balance_tab()

        # Set the layout on the main window
        self.setLayout(main_layout)

    # tabs #
    def setup_main_tab(self):
        layout = QVBoxLayout()

        # Load / Save retention file
        self.load_retention_button = QPushButton("Load Retention File")
        self.load_retention_button.clicked.connect(self.load_retention_table_dialog)
        layout.addWidget(self.load_retention_button)

        self.save_retention_button = QPushButton("Save New Retention Table")
        self.save_retention_button.clicked.connect(self.save_retention_table_dialog)
        layout.addWidget(self.save_retention_button)

        # Retention list
        layout.addWidget(QLabel("Retention Times (min):"))
        self.list_widget = QListWidget()
        self.list_widget.itemDoubleClicked.connect(self.remove_retention_time)
        self.refresh_retention_list()
        layout.addWidget(self.list_widget)

        # Add new retention time
        add_layout = QHBoxLayout()
        self.new_time_input = QLineEdit()
        self.new_time_input.setPlaceholderText("Add new retention time")
        self.new_name_input = QLineEdit()
        self.new_name_input.setPlaceholderText("Compound name")
        add_btn = QPushButton("Add")
        add_btn.clicked.connect(self.add_retention_time)

        self.tol_spin = QSpinBox()
        self.tol_spin.setRange(1, 1000)
        self.tol_spin.setValue(int(self.tolerance * 1000))
        self.tol_spin.setSuffix(" ms")
        self.tol_spin.valueChanged.connect(self.update_tolerance)

        # Retention time deviation
        add_layout.addWidget(self.new_time_input)
        add_layout.addWidget(self.new_name_input)
        add_layout.addWidget(add_btn)
        add_layout.addWidget(QLabel("Retention deviation:"))
        add_layout.addWidget(self.tol_spin)
        layout.addLayout(add_layout)

        # Calibration config
        flow_layout = QHBoxLayout()
        self.calib_flow_input = QDoubleSpinBox()
        self.calib_flow_input.setDecimals(3)
        self.calib_flow_input.setRange(0.001, 10.0)
        self.calib_flow_input.setSingleStep(0.001)
        self.calib_flow_input.setValue(self.calib_flow)
        self.calib_flow_input.setSuffix(" g/min")
        self.calib_flow_input.valueChanged.connect(self.update_calib_flow)

        self.experiment_flow_input = QDoubleSpinBox()
        self.experiment_flow_input.setDecimals(3)
        self.experiment_flow_input.setRange(0.001, 10.0)
        self.experiment_flow_input.setSingleStep(0.001)
        self.experiment_flow_input.setValue(self.experiment_flow)
        self.experiment_flow_input.setSuffix(" g/min")
        self.experiment_flow_input.valueChanged.connect(self.update_experiment_flow)

        self.slope_calibration_input = QLineEdit(str(self.slope_calibration))
        self.slope_calibration_input.textChanged.connect(self.update_slope_calibration)

        self.linear_co_calibration_input = QLineEdit(str(self.linear_co_calibration))
        self.linear_co_calibration_input.textChanged.connect(self.update_linear_co_calibration)

        self.reactant_rt_input = QLineEdit(str(self.reactant_rt))
        self.reactant_rt_input.textChanged.connect(self.update_reactant_rt)

        self.reactant_wt_input = QLineEdit(str(self.reactant_wt))
        self.reactant_wt_input.textChanged.connect(self.update_reactant_wt)

        flow_layout.addWidget(QLabel("Calibration Flow:"))
        flow_layout.addWidget(self.calib_flow_input)
        flow_layout.addWidget(QLabel("Experiment Flow:"))
        flow_layout.addWidget(self.experiment_flow_input)
        flow_layout.addWidget(QLabel("Slope Co.:"))
        flow_layout.addWidget(self.slope_calibration_input)
        flow_layout.addWidget(QLabel("Linear Co.:"))
        flow_layout.addWidget(self.linear_co_calibration_input)
        flow_layout.addWidget(QLabel("Reactant RT.:"))
        flow_layout.addWidget(self.reactant_rt_input)
        flow_layout.addWidget(QLabel("Reactant wt%.:"))
        flow_layout.addWidget(self.reactant_wt_input)
        layout.addLayout(flow_layout)

        # Reaction configuration
        self.load_config_button = QPushButton("Load reaction configuration")
        self.load_config_button.clicked.connect(self.load_reaction_config_dialog)
        layout.addWidget(self.load_config_button)

        # Import/Export
        btn_layout = QHBoxLayout()
        self.import_btn = QPushButton("Import PDFs")
        self.import_btn.clicked.connect(self.import_pdfs)
        self.save_btn = QPushButton("Save CSVs")
        self.save_btn.setEnabled(False)
        self.save_btn.clicked.connect(self.save_csvs)
        btn_layout.addWidget(self.import_btn)
        btn_layout.addWidget(self.save_btn)
        layout.addLayout(btn_layout)

        # Status
        self.status_label = QLabel("")
        layout.addWidget(self.status_label)

        # Graphs
        graph_layout = QVBoxLayout()

        self.figure1 = Figure(figsize=(8, 8))
        self.canvas1 = FigureCanvas(self.figure1)
        self.toolbar1 = NavigationToolbar(self.canvas1, self)
        graph_layout.addWidget(self.toolbar1)
        graph_layout.addWidget(self.canvas1)

        right_layout = QVBoxLayout()
        self.figure2 = Figure(figsize=(8, 8))
        self.canvas2 = FigureCanvas(self.figure2)
        self.toolbar2 = NavigationToolbar(self.canvas2, self)
        right_layout.addWidget(self.toolbar2)
        right_layout.addWidget(self.canvas2)
        graph_layout.addLayout(right_layout)
        layout.addLayout(graph_layout)

        self.main_tab.setLayout(layout)

    def setup_product_trend_tab(self):
        layout = QVBoxLayout()

        # Create matplotlib figure and canvas
        self.pt_figure = Figure(figsize=(8, 6))
        self.pt_canvas = FigureCanvas(self.pt_figure)
        self.toolbar3 = NavigationToolbar(self.pt_canvas, self)
        layout.addWidget(self.toolbar3)
        layout.addWidget(self.pt_canvas)

        # Create dropdown to select column
        self.area_selection_box = QComboBox()
        self.area_selection_box.addItems(["Area", "Area %"])
        layout.addWidget(self.area_selection_box)

        # Create and connect the button
        show_button = QPushButton("Show Results")
        show_button.clicked.connect(self.plot_product_trend)
        layout.addWidget(show_button)

        # Set layout to tab
        self.product_trend_tab.setLayout(layout)

    def setup_product_calibration_tab(self):
        layout = QVBoxLayout()

        title = QLabel("Calibration Setup (wt% vs A%)")
        title.setStyleSheet("font-weight: bold; font-size: 16px;")
        layout.addWidget(title)

        # Calibration table
        self.calibration_table = QTableWidget()
        self.calibration_table.setColumnCount(5)
        self.calibration_table.setHorizontalHeaderLabels(["Component Name", "Retention Time", "Molecular Formula", "Slope", "Intercept"])
        self.calibration_table.setRowCount(0)
        self.calibration_table.horizontalHeader().setStretchLastSection(True)
        self.calibration_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.calibration_table)

        # Plot area
        self.calibration_plot_canvas = FigureCanvas(Figure(figsize=(8, 6.5)))
        self.toolbar4 = NavigationToolbar(self.calibration_plot_canvas, self)
        layout.addWidget(self.toolbar4)
        layout.addWidget(self.calibration_plot_canvas)

        # Buttons

        # Load / Save retention file
        load_cali_button = QPushButton("Load Calibration File")
        load_cali_button.clicked.connect(self.load_calibration_dialog)

        calibrate_btn = QPushButton("Apply Calibration and Save CSV")
        calibrate_btn.clicked.connect(self.calibrate_and_export)

        plot_btn = QPushButton("Plot Calibrated Data")
        plot_btn.clicked.connect(self.plot_latest_calibration)

        button_layout = QHBoxLayout()
        button_layout.addWidget(load_cali_button)
        button_layout.addWidget(calibrate_btn)
        button_layout.addWidget(plot_btn)
        layout.addLayout(button_layout)

        # Status
        self.status_label_cali = QLabel("")
        layout.addWidget(self.status_label_cali)

        self.product_calibration_tab.setLayout(layout)

    def setup_molpercent_tab(self):
        layout = QVBoxLayout()

        # --- Title ---
        title = QLabel("Grouped Mol% Without Solvent")
        title.setStyleSheet("font-weight: bold; font-size: 16px;")
        layout.addWidget(title)

        # --- Matplotlib Canvas + Toolbar ---
        self.molpercent_canvas = FigureCanvas(Figure(figsize=(8, 6)))
        self.molpercent_toolbar = NavigationToolbar(self.molpercent_canvas, self)
        layout.addWidget(self.molpercent_toolbar)
        layout.addWidget(self.molpercent_canvas)

        # --- Plot Button ---
        plot_btn = QPushButton("Show Grouped Mol% Plot")
        plot_btn.clicked.connect(self.plot_grouped_molpercent)
        layout.addWidget(plot_btn)

        # --- Status Label ---
        self.molpercent_status = QLabel("")
        layout.addWidget(self.molpercent_status)

        # --- Set layout to tab ---
        self.molpercent_tab.setLayout(layout)

    def setup_selectivity_tab(self):
        layout = QVBoxLayout()

        # For selectivity calculation (exclude solvent and reactant)
        #_Reactant
        self.reactant_rt_input_2 = QLineEdit(str(self.reactant_rt))
        self.reactant_rt_input_2.textChanged.connect(self.update_reactant_rt_2)


        #Setting up the widgets
        layout.addWidget(QLabel("Reactant RT:"))
        layout.addWidget(self.reactant_rt_input_2)

        solvent_text = ", ".join(
            [f"{s.get('name', 'Solvent')} ({s.get('rt', 'n/a')})" for s in self.solvents]
        )
        layout.addWidget(QLabel(f"Loaded solvent RTs: {solvent_text}"))

        # Create matplotlib figure and canvas
        self.pt_figure_sel = Figure(figsize=(8, 6))
        self.pt_canvas_sel = FigureCanvas(self.pt_figure_sel)
        self.toolbar4 = NavigationToolbar(self.pt_canvas_sel, self)
        layout.addWidget(self.toolbar4)
        layout.addWidget(self.pt_canvas_sel)

        # Create matplotlib figure and canvas
        self.pt_figure_sel2 = Figure(figsize=(8, 6))
        self.pt_canvas_sel2 = FigureCanvas(self.pt_figure_sel2)
        self.toolbar5 = NavigationToolbar(self.pt_canvas_sel2, self)
        layout.addWidget(self.toolbar5)
        layout.addWidget(self.pt_canvas_sel2)

        # --- Create dropdown to select selectivity type ---
        self.selectivity_selection_box = QComboBox()
        self.selectivity_selection_box.addItems([
            "Selectivity_GC",
            "Selectivity_mol"
        ])
        layout.addWidget(self.selectivity_selection_box)

        # Create and connect the button
        show_button_sel = QPushButton("Show Results")
        show_button_sel.clicked.connect(self.plot_selectivity)
        layout.addWidget(show_button_sel)

        # Set layout to tab
        self.selectivity_tab.setLayout(layout)

    def setup_pearson_tab(self):
        layout = QVBoxLayout()

        # --- Button to calculate Pearson correlation ---
        btn_calc = QPushButton("Calculate Pearson Correlation")
        btn_calc.clicked.connect(self.calculate_and_plot_pearson)
        layout.addWidget(btn_calc)

        # --- Function to create diverging colormap centered at 0 ---
        def make_diverging_cmap(base_cmap_name="inferno"):
            base = plt.get_cmap(base_cmap_name)(np.linspace(0, 1, 256))
            colors = np.vstack([base[::-1], base])  # mirror for diverging effect
            return LinearSegmentedColormap.from_list(f"{base_cmap_name}_diverging", colors)

        # --- Custom colormaps (diverging versions of inferno/magma) ---
        self.custom_cmaps = {
            "inferno_diverging": make_diverging_cmap("inferno"),
            "magma_diverging": make_diverging_cmap("magma")
        }

        # --- ComboBox to select colormap ---
        self.cmap_box = QComboBox()
        self.cmap_box.addItems([
            "coolwarm","PuOr","BrBG","RdBu", "RdGy", "icefire", "twilight", "twilight_shifted","inferno_diverging", "magma_diverging"])
        self.cmap_box.currentTextChanged.connect(self.update_pearson_plot_cmap)
        layout.addWidget(self.cmap_box)

        # --- FigureCanvas for plotting ---
        self.pearson_figure = plt.Figure(figsize=(8, 6))
        self.pearson_canvas = FigureCanvas(self.pearson_figure)
        layout.addWidget(self.pearson_canvas)

        # --- Navigation Toolbar ---
        toolbar = NavigationToolbar(self.pearson_canvas, self)
        layout.addWidget(toolbar)

        self.pearson_tab.setLayout(layout)

    def setup_carbon_balance_tab(self):
        layout = QVBoxLayout()

        # Button to trigger plot
        plot_button = QPushButton("Show Carbon Recovery Plot")
        plot_button.clicked.connect(self.plot_carbon_recovery)
        layout.addWidget(plot_button)

        # --- Matplotlib canvas ---
        self.carbon_fig = plt.Figure(figsize=(8, 6))
        self.carbon_canvas = FigureCanvas(self.carbon_fig)
        layout.addWidget(self.carbon_canvas)

        # --- Navigation Toolbar ---
        toolbar = NavigationToolbar(self.carbon_canvas, self)
        layout.addWidget(toolbar)

        # Label to show save path or messages
        self.carbon_plot_label = QLabel("")
        layout.addWidget(self.carbon_plot_label)

        self.carbon_balance_tab.setLayout(layout)

    # -----------------------------
    # Config / loading functions
    # -----------------------------

    def load_reaction_config(self, filename):
        try:
            with open(filename, "r", encoding="utf-8") as f:
                cfg = json.load(f)

            self.tolerance = cfg.get("tolerance", self.tolerance)
            self.calib_flow = cfg.get("calib_flow", self.calib_flow)
            self.experiment_flow = cfg.get("experiment_flow", self.experiment_flow)

            self.slope_calibration = cfg.get("slope_calibration", self.slope_calibration)
            self.linear_co_calibration = cfg.get("linear_co_calibration", self.linear_co_calibration)

            self.reactant_name = cfg.get("reactant_name", self.reactant_name)
            self.reactant_rt = cfg.get("reactant_rt", self.reactant_rt)
            self.reactant_formula = cfg.get("reactant_formula", self.reactant_formula)
            self.reactant_molar_mass = cfg.get("reactant_molar_mass", self.reactant_molar_mass)
            self.reactant_nC = cfg.get("reactant_nC", self.reactant_nC)
            self.reactant_wt = cfg.get("reactant_wt", self.reactant_wt)

            self.solvents = cfg.get("solvents", self.solvents)

            self.key_products = cfg.get("key_products", self.key_products)
            self.groups = cfg.get("groups", self.groups)

            # update widgets if they already exist
            self.tol_spin.setValue(int(self.tolerance * 1000))
            self.calib_flow_input.setValue(self.calib_flow)
            self.experiment_flow_input.setValue(self.experiment_flow)
            self.slope_calibration_input.setText(str(self.slope_calibration))
            self.linear_co_calibration_input.setText(str(self.linear_co_calibration))
            self.reactant_rt_input.setText(str(self.reactant_rt))
            self.reactant_wt_input.setText(str(self.reactant_wt))

            if hasattr(self, "reactant_rt_input_2"):
                self.reactant_rt_input_2.setText(str(self.reactant_rt))

            self.status_label.setText("Reaction config loaded successfully.")

        except Exception as e:
            self.status_label.setText(f"Error loading reaction config: {e}")
            print(f"[DEBUG] Error loading reaction config: {e}")

    def load_reaction_config_dialog(self):
        file_name, _ = QFileDialog.getOpenFileName(
            self,
            "Select Reaction Config File",
            "",
            "JSON Files (*.json)"
        )
        if file_name:
            self.load_reaction_config(file_name)

    def get_solvent_rts(self):
        return [
            s["rt"] for s in self.solvents
            if isinstance(s, dict) and s.get("rt") is not None
        ]

    def compute_total_solvent_mol_flow(self, total_feed_rate_g_min):
        total_solvent_mol = 0.0

        for solvent in self.solvents:
            try:
                wt_fraction = float(solvent.get("wt_fraction", 0))
                mw = float(solvent.get("molar_mass", 0))
                if wt_fraction > 0 and mw > 0:
                    total_solvent_mol += (total_feed_rate_g_min * wt_fraction) / mw
            except Exception:
                continue

        return total_solvent_mol
    ## Functions tab_1 ##

    # Updates
    def update_calib_flow(self, value):
        self.calib_flow = value

    def update_experiment_flow(self, value):
        self.experiment_flow = value

    def update_slope_calibration(self, value):
        try:
            self.slope_calibration = float(value)
        except ValueError:
            pass

    def update_linear_co_calibration(self, value):
        try:
            self.linear_co_calibration = float(value)
        except ValueError:
            pass

    def update_tolerance(self):
        self.tolerance = self.tol_spin.value() / 1000.0
        self.status_label.setText(f"New retention deviation: {self.tolerance:.3f} min")

    def update_reactant_rt (self, value):
        try:
            self.reactant_rt = float(value)
        except ValueError:
            pass

    def update_reactant_wt(self, value):
        try:
            self.reactant_wt = float(value)
        except ValueError:
            pass



                    # Retention time

    def refresh_retention_list(self):
        self.list_widget.clear()
        for _, row in self.retention_times_df.sort_values("Retention Time").iterrows():
            self.list_widget.addItem(f"{row['Retention Time']:.5f} - {row['Compound']}")

    def load_retention_table_dialog(self):
        file_name, _ = QFileDialog.getOpenFileName(
            self, "Select Retention Table File", "", "CSV or TXT Files (*.csv *.txt)"
        )
        if file_name:
            self.load_retention_table(file_name)

    def load_retention_table(self, filename):
        try:
            if filename.endswith('.csv'):
                df = pd.read_csv(filename)
            elif filename.endswith('.txt'):
                df = pd.read_csv(filename, delimiter='\t', header=None, names=["Retention Time", "Compound"])
            else:
                raise ValueError("Unsupported format. Use .csv or .txt")

            if "Retention Time" not in df.columns or "Compound" not in df.columns:
                raise ValueError("File must contain 'Retention Time' and 'Compound' columns")

            self.retention_times_df = df
            self.refresh_retention_list()
            self.status_label.setText("Retention table loaded.")
        except Exception as e:
            self.status_label.setText(f"Error loading table: {e}")
        print(self.retention_times_df)

    def save_retention_table_dialog(self):
        file_name, selected_filter = QFileDialog.getSaveFileName(
            self, "Save Retention Table", "", "CSV Files (*.csv);;Text Files (*.txt)"
        )
        if file_name:
            try:
                if selected_filter == "CSV Files (*.csv)" or file_name.endswith('.csv'):
                    if not file_name.endswith('.csv'):
                        file_name += '.csv'
                    self.retention_times_df.to_csv(file_name, index=False)
                elif selected_filter == "Text Files (*.txt)" or file_name.endswith('.txt'):
                    if not file_name.endswith('.txt'):
                        file_name += '.txt'
                    self.retention_times_df.to_csv(file_name, index=False, sep='\t')
                else:
                    # Default to csv if extension/filter not recognized
                    if not file_name.endswith('.csv'):
                        file_name += '.csv'
                    self.retention_times_df.to_csv(file_name, index=False)

                self.status_label.setText(f"Retention table saved to {file_name}")
            except Exception as e:
                self.status_label.setText(f"Error saving table: {e}")

    def add_retention_time(self):
        try:
            time = float(self.new_time_input.text())
            name = self.new_name_input.text().strip()
            if not name:
                self.status_label.setText("Please enter a compound name.")
                return
            if not any(abs(time - t) < 1e-5 for t in self.retention_times_df["Retention Time"]):
                self.retention_times_df.loc[len(self.retention_times_df)] = [time, name]
                self.refresh_retention_list()
                self.new_time_input.clear()
                self.new_name_input.clear()
                self.status_label.setText(f"Time {time:.5f} ({name}) added.")
            else:
                self.status_label.setText("Time already in the list.")
        except ValueError:
            self.status_label.setText("Please enter a valid number.")

    def remove_retention_time(self, item):
        try:
            time_str = item.text().split(" - ")[0]
            val = float(time_str)
            self.retention_times_df = self.retention_times_df[
                (self.retention_times_df["Retention Time"].sub(val).abs() > 1e-5)
            ]
            self.refresh_retention_list()
            self.status_label.setText(f"Time {val:.5f} removed.")
        except Exception:
            self.status_label.setText("Error removing time.")


    # PDF processing
    def extract_front_signal_results(self, pdf_path):
        try:
            with fitz.open(pdf_path) as doc:
                all_text = "".join([page.get_text() for page in doc])
        except Exception as e:
            print(f"Error opening {pdf_path}: {e}")
            return pd.DataFrame()

        match = re.search(r'Front Signal\s+Results(.*?)(Back Signal Results|Aux Det|Totals|Acquired)', all_text, re.DOTALL)
        if not match:
            print(f"'Front Signal Results' section not found in {pdf_path}")
            return pd.DataFrame()

        section_text = re.sub(r'^.*Page.*\n?', '', match.group(1), flags=re.MULTILINE)
        numbers = re.findall(r'\d+[\.,]?\d*', section_text)
        numbers = [float(num.replace(',', '.')) for num in numbers]
        if len(numbers) % 5 != 0:
            print(f"Unexpected format in {pdf_path}, number of values: {len(numbers)}")

        rows = [numbers[i:i + 5] for i in range(0, len(numbers) - 4, 5)]
        return pd.DataFrame(rows, columns=["Retention Time", "Area", "Area %", "Height", "Height %"])

    def filter_peaks(self, df):
        return df[df["Retention Time"].apply(
            lambda x: any(abs(x - rt) <= self.tolerance for rt in self.retention_times_df["Retention Time"])
        )]

    def process_reactant_conversion(self, df):
        reactant_rt = self.reactant_rt
        wt_initial = self.reactant_wt

        correction_factor = self.calib_flow / self.experiment_flow

        df = df.copy()
        df["Sample Name"] = df.get("File", "").str.replace(".pdf", "", regex=False)
        df[f"wt% {self.reactant_name}"] = None
        df[f"{self.reactant_name} Conversion (%)"] = None

        for name, group in df.groupby("Sample Name"):
            match = group.loc[group["Retention Time"].sub(reactant_rt).abs() <= self.tolerance]

            if not match.empty:
                best_peak = match.loc[match["Area"].idxmax()]
                corrected_slope = self.slope_calibration * correction_factor
                wt_percent = (best_peak["Area"] - self.linear_co_calibration) / corrected_slope
                conversion = (wt_initial - wt_percent) / wt_initial * 100

                idx = best_peak.name
                df.at[idx, f"wt% {self.reactant_name}"] = round(wt_percent, 6)
                df.at[idx, f"{self.reactant_name} Conversion (%)"] = round(conversion, 4)

        df[f"wt% {self.reactant_name}"] = pd.to_numeric(
            df[f"wt% {self.reactant_name}"], errors="coerce"
        )
        df[f"{self.reactant_name} Conversion (%)"] = pd.to_numeric(
            df[f"{self.reactant_name} Conversion (%)"], errors="coerce"
        )

        return df

    def import_pdfs(self):
        try:
            file_paths, _ = QFileDialog.getOpenFileNames(self, "Select PDF files", filter="PDF Files (*.pdf)")
            if not file_paths:
                self.status_label.setText("Import cancelled.")
                return

            # Create CSV folder for results output
            base_folder = os.path.dirname(file_paths[0])

            csv_folder = os.path.join(base_folder, "CSV_Files")
            os.makedirs(csv_folder, exist_ok=True)

            # Create Plot summary folder
            plot_folder = os.path.join(base_folder, "Plot_Summary")
            os.makedirs(plot_folder, exist_ok=True)

            self.csv_folder = csv_folder
            self.plot_folder = plot_folder

            def extract_numeric_part(filename):
                base = os.path.basename(filename).replace(".pdf", "")
                numeric_part = re.sub(r'[^\d.]', '', base)
                try:
                    return float(numeric_part)
                except ValueError:
                    return float('inf')

            all_data = []
            for path in sorted(file_paths, key=extract_numeric_part):
                df = self.extract_front_signal_results(path)
                if not df.empty:
                    df["File"] = os.path.basename(path)
                    df["Sample Name"] = df["File"].replace(".pdf", "")
                    all_data.append(df)

            if all_data:
                self.pdf_data = pd.concat(all_data, ignore_index=True)
                self.status_label.setText(f"{len(all_data)} PDFs successfully imported.")
                self.save_btn.setEnabled(True)
                self.plot_conversion()
            else:
                self.status_label.setText("No data extracted.")
        except Exception as e:
            print(f"General error: {e}")

    def save_csvs(self):
        if self.pdf_data.empty or not self.csv_folder:
            self.status_label.setText("Nothing to save or unknown source folder.")
            return

        all_path = os.path.join(self.csv_folder, "complete.csv")
        self.pdf_data.to_csv(all_path, index=False)

        self.filtered_data = self.filter_peaks(self.pdf_data.copy())
        self.filtered_data["Sample Name"] = self.filtered_data["File"].str.replace(".pdf", "")
        self.filtered_data = self.process_reactant_conversion(self.filtered_data)

        filtered_path = os.path.join(self.csv_folder, "filtered.csv")
        self.filtered_data.to_csv(filtered_path, index=False)
        self.status_label.setText("CSVs saved to the PDFs folder.")
        self.plot_conversion()

    def plot_conversion(self):
        try:
            if self.pdf_data.empty:
                return

            df = self.filtered_data if not self.filtered_data.empty else self.pdf_data
            conv_col = f"{self.reactant_name} Conversion (%)"
            df = df[df[conv_col].notna()].copy()
            if df.empty:
                return

            df["Sample Name"] = df["Sample Name"].astype(float)
            grouped = df.groupby("Sample Name")[conv_col].mean().sort_index()
            grouped_save = grouped

            grouped.index.name = "Time on stream (h)"
            print(grouped_save)

            conversion_path = os.path.join(self.csv_folder, "conversion.csv")
            grouped_save.to_csv(conversion_path, index=True)
            self.status_label.setText("Conversion.csv saved to the PDF folder.")

            # Plot 1
            self.figure1.clear()
            ax1 = self.figure1.add_subplot(111)
            ax1.plot(grouped.index, grouped.values, marker='o', linestyle='-', color='#9C2F45')
            ax1.set_title(f"{self.reactant_name} Conversion (0 to 100%)", fontname='Arial', fontsize=12)
            ax1.set_ylabel("Conversion (%)",fontname='Arial', fontsize=12)
            ax1.set_xlabel("TOS (h)", fontname='Arial', fontsize=12)
            ax1.set_ylim(0, 100)
            ax1.grid(True)

            self.canvas1.draw()

            # Save graph 1 with temporary size
            original_size1 = self.figure1.get_size_inches()
            self.figure1.set_size_inches(8, 6)
            photo_path_1 = os.path.join(self.plot_folder, "conversion_plot.png")
            self.figure1.savefig(photo_path_1, dpi=300, bbox_inches='tight')
            self.figure1.set_size_inches(*original_size1)  # Restore
            self.canvas1.draw()


            # Plot 2
            self.figure2.clear()
            ax2 = self.figure2.add_subplot(111)
            ax2.plot(grouped.index, grouped.values, marker='o', linestyle='-', color='#4954B0')
            ax2.set_title(f"{self.reactant_name} Conversion (zoom)", fontname='Arial', fontsize=12)
            ax2.set_ylabel("Conversion (%)",fontname='Arial', fontsize=12)
            ax2.set_xlabel("TOS (h)", fontname='Arial', fontsize=12)
            ax2.grid(True)

            self.canvas2.draw()

            # Save graph 2 with temporary size
            original_size2 = self.figure2.get_size_inches()
            self.figure2.set_size_inches(8, 6)
            photo_path_2 = os.path.join(self.plot_folder, "conversion_zoom_plot.png")
            self.figure2.savefig(photo_path_2, dpi=300, bbox_inches='tight')
            self.figure2.set_size_inches(*original_size2)
            self.canvas2.draw()

        except Exception as e:
            print(f"Error plotting: {e}")

    ## Functions tab_2 ##

    def plot_product_trend(self):
        self.variable_a = self.area_selection_box.currentText()

        if self.filtered_data.empty:
            QMessageBox.warning(self, "Warning", "Empty dataset")
            return

        df = self.filtered_data.copy()
        df["Sample Name"] = df["Sample Name"].astype(float)

        trend_data = {}
        for rt, comp in self.retention_times_df.itertuples(index=False):
            # Select rows within tolerance of the target retention time
            tmp = df[df["Retention Time"].sub(rt).abs() <= self.tolerance]
            if tmp.empty:
                continue

            # Group by Sample Name and calculate mean Area% (or the column you have)
            grouped = tmp.groupby("Sample Name")[self.variable_a].mean().sort_index()

            # Align by subtracting the first value (optional, you may want absolute values)
            first = grouped.iloc[0]
            aligned = grouped #- first (minus to make it relative to the first sample)

            # Store the aligned area percentage trend for this compound
            trend_data[comp] = aligned

        self.pt_figure.clear()
        ax = self.pt_figure.add_subplot(111)

        for comp, series in sorted(trend_data.items()):
            ax.plot(series.index, series.values, marker='o', linestyle='-', label=comp)

        ax.set_xlabel("Time on stream (h)")
        ax.set_ylabel(self.variable_a+" (a.u.)")
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
        self.pt_figure.tight_layout()
        ax.grid(True)

        self.pt_canvas.draw()


    ## Functions tab_3 ##

    def load_calibration_dialog(self):
        file_name, _ = QFileDialog.getOpenFileName(
            self,
            "Select Calibration Table File",
            "",
            "CSV or TXT Files (*.csv *.txt)"
        )
        if file_name:
            self.load_calibration_table(file_name)

    def load_calibration_table(self, filename):
        try:
            # --- Read file ---
            if filename.endswith(".csv"):
                df = pd.read_csv(filename)
            elif filename.endswith(".txt"):
                # Assume tab-separated and no header
                df = pd.read_csv(
                    filename,
                    sep="\t",
                    header=None,
                    names=["RT", "Compound", "Formula", "Slope", "Intercept"]
                )
            else:
                QMessageBox.warning(self, "Warning", "Unsupported format. Use .csv or .txt")
                return

            # --- Validate columns ---
            required_cols = ["RT", "Compound", "Formula", "Slope", "Intercept"]
            if not all(col in df.columns for col in required_cols):
                QMessageBox.warning(
                    self,
                    "Warning",
                    f"File must contain columns: {', '.join(required_cols)}"
                )
                return

            # --- Save and update status ---
            self.uploaded_calibration_df = df
            self.status_label_cali.setText("Calibration table loaded successfully.")
            print(self.uploaded_calibration_df)

            # --- Fill the QTableWidget ---
            self.calibration_table.setRowCount(len(df))

            for row, record in enumerate(df.itertuples(index=False)):
                name_item = QTableWidgetItem(str(record.Compound))
                name_item.setFlags(name_item.flags() & ~Qt.ItemFlag.ItemIsEditable) # flag to avoid edition

                rt_item = QTableWidgetItem(str(record.RT))
                rt_item.setFlags(rt_item.flags() & ~Qt.ItemFlag.ItemIsEditable) # flag to avoid edition

                formula_item = QTableWidgetItem(str(record.Formula))

                slope_item = QTableWidgetItem(str(record.Slope))

                intercept_item = QTableWidgetItem(str(record.Intercept))

                self.calibration_table.setItem(row, 0, name_item)
                self.calibration_table.setItem(row, 1, rt_item)
                self.calibration_table.setItem(row, 2, formula_item)
                self.calibration_table.setItem(row, 3, slope_item)
                self.calibration_table.setItem(row, 4, intercept_item)

            print("[DEBUG] Calibration table populated successfully.")

        except Exception as e:
            self.status_label_cali.setText(f"Error loading calibration table: {e}")
            print("[DEBUG] Error loading calibration table:", e)

    def calibrate_and_export(self):
        if self.filtered_data.empty or self.retention_times_df.empty:
            QMessageBox.warning(self, "Warning", "Filtered data or retention times data is empty.")
            return

        # --- Step 1: Read calibration table into dictionary ---
        calibration_data = {}
        for row in range(self.calibration_table.rowCount()):
            rt_item = self.calibration_table.item(row, 1)
            formula_item = self.calibration_table.item(row, 2)
            slope_item = self.calibration_table.item(row, 3)
            intercept_item = self.calibration_table.item(row, 4)

            if not rt_item or not slope_item or not intercept_item:
                continue

            try:
                rt_key = float(rt_item.text().replace(',', '.'))
                formula = formula_item.text()
                slope = float(slope_item.text().replace(',', '.'))
                intercept = float(intercept_item.text().replace(',', '.'))
                calibration_data[rt_key] = (formula, slope, intercept)
            except ValueError:
                continue

        if not calibration_data:
            QMessageBox.warning(self, "Warning", "No valid calibration data found.")
            return

        print(f"[DEBUG] Loaded calibration data: {calibration_data}")
        print(f"[DEBUG] Tolerance set to: {self.tolerance}")

        # --- Step 2: Work on a copy of the filtered data ---
        df = self.filtered_data.copy().reset_index(drop=True)

        # Ensure 'wt%' column exists
        if "wt%" not in df.columns:
            df["wt%"] = ""

        # Ensure 'wt%' column exists
        if "Molecular formula" not in df.columns:
            df["Molecular formula"] = ""

        # --- Step 3: Fill calibrated values only when match exists ---
        for index, row in df.iterrows():
            try:
                rt = float(row["Retention Time"])
                area = float(row["Area"])
            except Exception:
                continue  # leave blank if parsing fails

            # Initialize as blank
            wt_percent = ""
            molecular_formula = ""

            # Match closest RT within tolerance
            for rt_ref, (formula, slope, intercept) in calibration_data.items():
                if np.isclose(rt, rt_ref, atol=self.tolerance):
                    if slope != 0:
                        wt_percent = (area - intercept) / slope
                        molecular_formula = formula
                    break

            # Assign only if matched; otherwise remains blank
            df.at[index, "wt%"] = wt_percent
            df.at[index, "Molecular formula"] = molecular_formula

        # --- Step 4: Compute g/min, MW, mol/min, mol%, and gC/min per group with solvent reference ---
        total_feed_rate_g_min = self.experiment_flow  # feed rate in g/min

        # Solvent information
        # self.solvents must be a list of dicts with:
        # name, rt, molar_mass, wt_fraction
        # Example:
        # self.solvents = [
        #     {"name": "Decalin trans", "rt": 20.062, "molar_mass": 138.23, "wt_fraction": 0.45},
        #     {"name": "Decalin cis", "rt": 21.208, "molar_mass": 138.23, "wt_fraction":  0.45}
        # ]

        reactant_fraction_wt = self.reactant_wt / 100.0
        total_solvent_fraction_wt = sum(
            float(s.get("wt_fraction", 0)) for s in self.solvents
        )

        # Optional consistency check
        total_fraction = reactant_fraction_wt + total_solvent_fraction_wt
        if not np.isclose(total_fraction, 1.0, atol=1e-3):
            QMessageBox.warning(
                self,
                "Warning",
                f"Reactant + solvent weight fractions do not sum to 1.0 (current total = {total_fraction:.4f})."
            )
            return

        # Ensure wt% column exists and is numeric
        df["wt%"] = pd.to_numeric(df.get("wt%", 0), errors="coerce").fillna(0).clip(lower=0)

        # --- 4a: Compute g/min for measured compounds (non-solvent) ---
        df["g/min"] = (df["wt%"] / 100) * total_feed_rate_g_min

        # --- 4b: Compute molecular weight ---
        def compute_mw(formula):
            if not formula or formula.strip() == "":
                return np.nan
            try:
                pattern = r'([A-Z][a-z]*)(\d*)'
                matches = re.findall(pattern, formula)
                mw = 0.0
                for elem, count in matches:
                    count = int(count) if count else 1
                    mw += getattr(pt, elem).mass * count
                return mw
            except Exception as e:
                print(f"[ERROR] Could not parse formula '{formula}': {e}")
                return np.nan

        df["MW"] = df["Molecular formula"].apply(compute_mw)

        # --- 4c: Compute mol/min ---
        df["mol/min"] = df["g/min"] / df["MW"]
        df["mol/min"] = df["mol/min"].fillna(0).clip(lower=0)

        # --- helper: total solvent molar flow ---
        def compute_total_solvent_mol_flow():
            total_solvent_mol = 0.0

            for solvent in self.solvents:
                try:
                    wt_fraction = float(solvent.get("wt_fraction", 0))
                    mw = float(solvent.get("molar_mass", 0))

                    if wt_fraction > 0 and mw > 0:
                        solvent_g_min = total_feed_rate_g_min * wt_fraction
                        total_solvent_mol += solvent_g_min / mw
                except Exception:
                    continue

            return total_solvent_mol

        total_solvent_mol = compute_total_solvent_mol_flow()

        # --- 4d: Compute mol% per Sample Name group including solvent implicitly ---
        def compute_group_molperc(group):
            total_mol = group["mol/min"].sum() + total_solvent_mol
            if total_mol > 0:
                group["mol%"] = group["mol/min"] / total_mol * 100
            else:
                group["mol%"] = 0
            return group

        df = df.groupby("Sample Name", group_keys=False).apply(compute_group_molperc)

        # --- 4e: Compute mol% per Sample Name group without solvent implicitly but with reactant---
        def compute_group_molperc_ws(group):
            total_mol = group["mol/min"].sum()
            if total_mol > 0:
                group["mol%_without_solv"] = group["mol/min"] / total_mol * 100
            else:
                group["mol%_without_solv"] = 0
            return group

        df = df.groupby("Sample Name", group_keys=False).apply(compute_group_molperc_ws)

        # --- 4f: Compute mol% normalized(products only, sum = 100% per Sample Name, exclude reactant) --> the same as selectivity mol%
        def compute_molperc_normalized_products(group):
            # Exclude the reactant
            products = group[group["Molecular formula"] != self.reactant_formula]
            total_mol_products = products["mol/min"].sum()

            if total_mol_products > 0:
                # Assign normalized mol% only for products
                group["mol%_normalized"] = (group["mol/min"] / total_mol_products) * 100
                # Set reactant row to NaN for normalized mol%
                group.loc[group["Molecular formula"] == self.reactant_formula, "mol%_normalized"] = np.nan
            else:
                group["mol%_normalized"] = np.nan
            return group

        df = df.groupby("Sample Name", group_keys=False).apply(compute_molperc_normalized_products)

        # --- 4g: Compute gC/min ---
        def gmin_carbon(row):
            formula = row["Molecular formula"]
            mol_per_min = row["mol/min"]
            if mol_per_min <= 0 or not formula:
                return 0.0
            try:
                matches = re.findall(r'(C)(\d*)', formula)
                nC = sum(int(count) if count else 1 for elem, count in matches)
                return mol_per_min * nC * pt.C.mass
            except:
                return 0.0

        df["gC/min"] = df.apply(gmin_carbon, axis=1)

        # --- 4h: Carbon balance per Sample Name ---
        # Sum total output carbon flow per sample (reactant + products)

        carbon_summary = (
            df.groupby("Sample Name", as_index=False)["gC/min"]
              .sum()
              .rename(columns={"gC/min": "Total gC_out/min"})
        )

        # Compute dynamic reactant carbon input (based on any specified reactant)
        reactant_flow_g_min = total_feed_rate_g_min * (self.reactant_wt / 100.0)

        # Compute molecular weight of reactant consistently
        reactant_mw = compute_mw(self.reactant_formula) if hasattr(self, "reactant_formula") else getattr(self, "reactant_molar_mass", np.nan)

        # Compute number of carbons in reactant formula (if not explicitly set)
        if hasattr(self, "reactant_nC") and self.reactant_nC:
            reactant_nC = self.reactant_nC
        else:
            try:
                matches = re.findall(r'(C)(\d*)', self.reactant_formula)
                reactant_nC = sum(int(count) if count else 1 for _, count in matches)
            except:
                reactant_nC = 0

        # Compute molar and carbon inlet flow
        reactant_mol_min = reactant_flow_g_min / reactant_mw if reactant_mw and reactant_mw > 0 else 0.0
        gC_in_min = reactant_mol_min * reactant_nC * pt.C.mass

        # Attach inlet and recovery columns
        carbon_summary["gC_in/min"] = gC_in_min
        carbon_summary["Carbon recovery (%)"] = (
            carbon_summary["Total gC_out/min"] / gC_in_min * 100
        ).replace([np.inf, -np.inf], 0).fillna(0).clip(lower=0)

        # Merge summary back into df
        df = df.merge(carbon_summary, on="Sample Name", how="left")

        # Drop the intermediate columns if you only want Carbon recovery (%)
        df = df.drop(columns=["Total gC_out/min", "gC_in/min"], errors="ignore")

        # --- Keep Carbon recovery only in the first row of each Sample Name group ---
        for col in ["Carbon recovery (%)"]:
            df[col] = df.groupby("Sample Name")[col].transform(
                lambda x: x.mask(~x.index.isin([x.index[0]]))
            )

        # --- Optional: reorder columns so MW is next to Molecular formula ---
        cols = df.columns.tolist()
        if "MW" in cols and "Molecular formula" in cols:
            mf_idx = cols.index("Molecular formula")
            cols.insert(mf_idx + 1, cols.pop(cols.index("MW")))
            df = df[cols]

        # --- Step 5: Save calibrated CSV, preserving all original columns ---
        output_path = os.path.join(self.csv_folder, "calibrated.csv")
        df.to_csv(output_path, index=False)

        # --- Step 6: Store for plotting ---
        self.latest_calibrated_df = df

        print(f"[DEBUG] Saved calibrated file to: {output_path}")
        print("[DEBUG] Preview of calibrated data:")
        print(df.head(10))

    def plot_latest_calibration(self):
        if not hasattr(self, 'latest_calibrated_df') or self.latest_calibrated_df.empty:
            QMessageBox.warning(self, "Warning", "No calibrated data available to plot.")
            return

        df = self.latest_calibrated_df.copy()
        required_cols = ["Sample Name", "mol%_without_solv", "Retention Time"]
        if not all(col in df.columns for col in required_cols):
            QMessageBox.warning(self, "Warning", f"Required columns {required_cols} not found in calibrated data.")
            return

        # --- Read calibration table (RT, slope, intercept, compound) ---
        calibration_data = {}
        calibration_labels = {}
        for row in range(self.calibration_table.rowCount()):
            compound_item = self.calibration_table.item(row, 0)
            rt_item = self.calibration_table.item(row, 1)
            slope_item = self.calibration_table.item(row, 3)
            intercept_item = self.calibration_table.item(row, 4)

            try:
                compound = compound_item.text() if compound_item else "Unknown"
                rt = float(rt_item.text().replace(',', '.'))
                slope = float(slope_item.text().replace(',', '.')) if slope_item else None
                intercept = float(intercept_item.text().replace(',', '.')) if intercept_item else None
                calibration_data[rt] = (slope, intercept)
                calibration_labels[rt] = compound
            except Exception:
                continue

        # Prepare plotting canvas
        self.calibration_plot_canvas.figure.clear()
        ax = self.calibration_plot_canvas.figure.add_subplot(111)

        # Ensure numeric
        df["Sample Name"] = pd.to_numeric(df["Sample Name"], errors="coerce")

        # Sort by RT for consistency
        df = df.sort_values("Retention Time").reset_index(drop=True)

        # --- NEW TOLERANCE-BASED GROUPING (robust, order-independent) ---
        used = set()
        groups = []
        rt_values = df["Retention Time"].values

        for i, rt in enumerate(rt_values):
            if i in used:
                continue

            # Find all rows within tolerance of this RT
            idx_group = df.index[df["Retention Time"].sub(rt).abs() <= self.tolerance]
            used.update(idx_group)
            groups.append(df.loc[idx_group].copy())

        # --- Plot each matched/calibrated group ---
        for group in groups:
            rt_mean = group["Retention Time"].mean()

            # Find matching calibration peak (same logic as calibration routine)
            matched_rt = None
            for rt_cal in calibration_data.keys():
                if abs(rt_mean - rt_cal) <= self.tolerance:
                    matched_rt = rt_cal
                    break

            if matched_rt is None:
                continue  # No calibrated equivalent → skip

            slope, intercept = calibration_data[matched_rt]
            compound = calibration_labels.get(matched_rt, "Unknown")

            # Only plot components with valid calibration
            if slope is None or intercept is None:
                continue

            # Ensure mol% is numeric
            group_plot = group[pd.to_numeric(group["mol%_without_solv"], errors="coerce").notna()].copy()
            if group_plot.empty:
                continue

            group_plot = group_plot.sort_values("Sample Name")

            ax.plot(
                group_plot["Sample Name"],
                group_plot["mol%_without_solv"],
                marker="o",
                linestyle="-",
                label=compound,
            )

        # Labels & styling
        ax.set_title("Calibrated mol% vs time on stream")
        ax.set_xlabel("Time on stream (h)")
        ax.set_ylabel("mol%")
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
        self.calibration_plot_canvas.figure.subplots_adjust(right=0.75)
        ax.grid(True)

        # Redraw the canvas
        self.calibration_plot_canvas.draw()

    ## Functions tab_4 ##

    def plot_grouped_molpercent(self):
        # --- Check for data ---
        if not hasattr(self, "latest_calibrated_df") or self.latest_calibrated_df.empty:
            QMessageBox.warning(self, "Warning", "No calibrated data available to plot.")
            return

        df = self.latest_calibrated_df.copy()
        df["Sample Name"] = pd.to_numeric(df["Sample Name"], errors="coerce")
        df["mol%_without_solv"] = pd.to_numeric(df["mol%_without_solv"], errors="coerce").fillna(0)

        tol = getattr(self, "tolerance", 0.05)  # fallback tolerance
        exclude_rts = [getattr(self, "reactant_rt", None),
                       getattr(self, "solvent_1_rt", None),
                       getattr(self, "solvent_2_rt", None)]
        exclude_rts = [rt for rt in exclude_rts if rt is not None]

        trend_data = []

        # --- Iterate over retention times and compounds ---
        for rt, comp_name in self.retention_times_df.itertuples(index=False):
            # skip reactants/solvents
            if any(abs(rt - ex_rt) <= tol for ex_rt in exclude_rts):
                continue

            # select rows within tolerance
            tmp = df[df["Retention Time"].sub(rt).abs() <= tol]
            if tmp.empty:
                continue

            grouped = tmp.groupby("Sample Name")[["mol%_without_solv"]].sum().reset_index()
            grouped["Component"] = comp_name  # standardize column name
            trend_data.append(grouped)

        if not trend_data:
            print("No product peaks found within tolerance.")
            return

        # --- Concatenate all compounds ---
        result = pd.concat(trend_data, ignore_index=True)
        result["mol_percentage"] = result["mol%_without_solv"]

        # --- Assign groups ---
        result["Group"] = None

        for group_name, compounds in self.groups.items():
            result.loc[result["Component"].isin(compounds), "Group"] = group_name

        for comp in self.key_products:
            result.loc[result["Component"] == comp, "Group"] = comp

        result["Group"] = result["Group"].fillna("Other")

        # --- Aggregate by Sample Name and Group ---
        grouped_result = (
            result.groupby(["Sample Name", "Group"])[["mol_percentage"]]
            .sum()
            .reset_index()
        )

        # --- Compute averages for reference ---
        avg_df = grouped_result.groupby("Group")[["mol_percentage"]].mean().reset_index()
        avg_df["Sample Name"] = "Average"

        # --- Setup plot style ---
        import matplotlib.pyplot as plt
        import matplotlib as mpl

        mpl.rcParams.update({
            "font.family": "Arial",
            "font.size": 12,
            "axes.titlesize": 12,
            "axes.labelsize": 12,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.fontsize": 12,
        })

        # --- generate colormap colors ---
        unique_groups = grouped_result["Group"].unique()
        cmap = plt.get_cmap("icefire", len(unique_groups))
        colors = {grp: cmap(i / (len(unique_groups) - 1)) for i, grp in enumerate(unique_groups)}

        # --- target color in RGB (0-1) ---
        target_rgb = np.array([0.200, 0.125, 0.137])

        # --- replace if close ---
        for grp, col in colors.items():
            rgb = np.array(col[:3])  # ignore alpha
            if np.allclose(rgb, target_rgb, atol=0.02):  # allow tiny difference
                colors[grp] = "#C0C0C0"

        # --- Trend plot ---
        self.molpercent_canvas.figure.clear()
        ax = self.molpercent_canvas.figure.add_subplot(111)
        for grp, grp_data in grouped_result.groupby("Group"):
            ax.plot(
                grp_data["Sample Name"],
                grp_data["mol_percentage"],
                marker="o",
                linestyle="-",
                label=grp,
                color=colors[grp]
            )
        # Title with Arial 12 explicitly
        ax.set_title("Mol percentage over time on stream", fontsize=12, fontname='Arial')
        ax.set_xlabel("Time on stream (h)")
        ax.set_ylabel("Mol (%)")
        # Rotate x-axis labels
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=12, fontname='Arial')
        # Defined y-axis labels
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=12, fontname='Arial')
        ax.legend()
        ax.grid(True)
        self.molpercent_canvas.draw()

        # --- Save figure ---
        fig = self.molpercent_canvas.figure
        original_size = fig.get_size_inches()
        fig.set_size_inches(8, 6)
        photo_path = os.path.join(self.plot_folder, "mol_percentage.png")
        fig.savefig(photo_path, dpi=300, bbox_inches='tight')
        fig.set_size_inches(*original_size)
        self.molpercent_canvas.draw()

    ## Functions tab_5 ## removed ##

    ## Functions tab_6 ##

    def update_reactant_rt_2 (self, value):
        try:
            self.reactant_rt = float(value)
        except ValueError:
            pass

    def plot_selectivity(self):
        selected_option = self.selectivity_selection_box.currentText()

        if self.latest_calibrated_df.empty:
            QMessageBox.warning(self, "Warning", "Empty dataset")
            return

        df = self.latest_calibrated_df.copy()
        df["Sample Name"] = df["Sample Name"].astype(float)

        tol = self.tolerance
        exclude_rts = [self.reactant_rt] + self.get_solvent_rts()
        exclude_rts = [rt for rt in exclude_rts if rt is not None]

        trend_data = []

        # --- iterate through compounds ---
        for rt, comp_name in self.retention_times_df.itertuples(index=False):
            if any(abs(rt - ex_rt) <= tol for ex_rt in exclude_rts):
                continue

            tmp = df[df["Retention Time"].sub(rt).abs() <= tol]
            if tmp.empty:
                continue

            grouped = tmp.groupby("Sample Name")[["Area %", "wt%", "mol%", "mol/min"]].sum().reset_index()
            grouped["Component"] = comp_name
            trend_data.append(grouped)

        if not trend_data:
            print("No product peaks found")
            return

        result = pd.concat(trend_data)

        # --- vectorized selectivity calculation ---
        totals = result.groupby("Sample Name")[["Area %", "wt%", "mol%"]].transform("sum")
        result["Selectivity_GC"] = result["Area %"] / totals["Area %"] * 100
        result["Selectivity_wt%"] = result["wt%"] / totals["wt%"] * 100
        result["Selectivity_mol"] = result["mol%"] / totals["mol%"] * 100

        # --- conversion-based selectivity using mol/min ---
        n_reactant_feed = (self.experiment_flow * self.reactant_wt / 100) / self.reactant_molar_mass  # mol/min

        # reactant mol/min per sample
        reactant_data = df[df["Retention Time"].sub(self.reactant_rt).abs() <= tol]
        if reactant_data.empty:
            print("⚠️ Reactant peak not found — conversion-based selectivity will be NaN.")
            result["Selectivity_mol (conversion-based)"] = np.nan
        else:
            reactant_grouped = (
                reactant_data.groupby("Sample Name")[["mol/min"]]
                .sum()
                .reset_index()
                .rename(columns={"mol/min": "mol_min_reactant"})
            )

            # Merge reactant mol/min into product table
            result = result.merge(reactant_grouped, on="Sample Name", how="left")

            # Conversion-based selectivity: mol/min product / (n_reactant_feed - mol/min reactant)
            result["Selectivity_mol (conversion-based)"] = result["mol/min"] / (
                    n_reactant_feed - result["mol_min_reactant"]
            ) * 100

        # ======================================================
        # === GROUPING STEP ====================================
        # ======================================================
        result["Group"] = None

        # assign each component to a group
        for group_name, compounds in self.groups.items():
            result.loc[result["Component"].isin(compounds), "Group"] = group_name

        # assign key products as their own group
        for comp in self.key_products:
            result.loc[result["Component"] == comp, "Group"] = comp

        # everything else as “Other”
        result["Group"] = result["Group"].fillna("Other")

        # --- compute selectivity by group ---
        grouped_result = result.groupby(["Sample Name", "Group"])[
            ["Selectivity_GC", "Selectivity_wt%", "Selectivity_mol", "Selectivity_mol (conversion-based)"]
        ].sum().reset_index()

        # --- compute average over time ---
        avg_df = grouped_result.groupby("Group")[
            ["Selectivity_GC", "Selectivity_wt%", "Selectivity_mol", "Selectivity_mol (conversion-based)"]
        ].mean().reset_index()
        avg_df["Sample Name"] = "Average"

        # ======================================================
        # === PLOTS ============================================
        # ======================================================
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import numpy as np
        import os

        # set global style to Arial 12
        mpl.rcParams.update({
            "font.family": "Arial",
            "font.size": 12,
            "axes.titlesize": 12,
            "axes.labelsize": 12,
            "xtick.labelsize": 12,
            "ytick.labelsize": 12,
            "legend.fontsize": 12,
        })

        # --- generate colormap colors ---
        unique_groups = grouped_result["Group"].unique()
        cmap = plt.get_cmap("icefire", len(unique_groups))
        colors = {grp: cmap(i / (len(unique_groups) - 1)) for i, grp in enumerate(unique_groups)}

        # --- replace target color if close ---
        target_rgb = np.array([0.200, 0.125, 0.137])
        for grp, col in colors.items():
            rgb = np.array(col[:3])  # ignore alpha
            if np.allclose(rgb, target_rgb, atol=0.02):
                colors[grp] = "#C0C0C0"

        # ======================================================
        # === TREND PLOT =======================================
        # ======================================================
        self.pt_figure_sel.clear()
        ax = self.pt_figure_sel.add_subplot(111)
        for grp, grp_data in grouped_result.groupby("Group"):
            ax.plot(
                grp_data["Sample Name"],
                grp_data[selected_option],
                marker="o",
                linestyle="-",
                label=grp,
                color=colors[grp]
            )
        ax.set_xlabel("Time on stream (h)")
        ax.set_ylabel("Selectivity (%)")
        ax.set_title("Selectivity over time on stream", fontsize=12, fontname='Arial')
        ax.set_xticklabels(ax.get_xticklabels(), fontsize=12, fontname='Arial')
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=12, fontname='Arial')
        ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5))
        self.pt_figure_sel.tight_layout()
        ax.grid(True)
        self.pt_canvas_sel.draw()

        # --- save trend figure ---
        fig1 = self.pt_canvas_sel.figure
        original_size1 = fig1.get_size_inches()
        fig1.set_size_inches(8, 6)
        photo_path1 = os.path.join(self.plot_folder, "trend_selectivity_group_plot.png")
        fig1.savefig(photo_path1, dpi=300, bbox_inches='tight')
        fig1.set_size_inches(*original_size1)
        self.pt_canvas_sel.draw()

        # ======================================================
        # === AVERAGE BAR PLOT =================================
        # ======================================================
        self.pt_figure_sel2.clear()
        ax2 = self.pt_figure_sel2.add_subplot(111)
        bars = ax2.bar(
            avg_df["Group"],
            avg_df[selected_option],
            color=[colors[g] for g in avg_df["Group"]],
            edgecolor="black"
        )
        ax2.set_ylabel("Averaged selectivity (%)")
        ax2.set_xticklabels(avg_df["Group"], rotation=45, ha="right", fontsize=12, fontname='Arial')
        ax2.set_title("Averaged selectivity by group", fontsize=12, fontname='Arial')
        ax2.set_yticklabels(ax2.get_yticklabels(), fontsize=12, fontname='Arial')
        self.pt_canvas_sel2.draw()

        # --- save average figure ---
        fig = self.pt_canvas_sel2.figure
        original_size = fig.get_size_inches()
        fig.set_size_inches(8, 6)
        photo_path = os.path.join(self.plot_folder, "avg_selectivity_group_plot.png")
        fig.savefig(photo_path, dpi=300, bbox_inches='tight')
        fig.set_size_inches(*original_size)
        self.pt_canvas_sel2.draw()

        # --- save CSV results ---
        csv_paths = {
            "full": os.path.join(self.csv_folder, "selectivity_group_results.csv"),
            "avg": os.path.join(self.csv_folder, "selectivity_group_avg_results.csv")
        }

        grouped_result.to_csv(csv_paths["full"], index=False)
        avg_df.to_csv(csv_paths["avg"], index=False)
        print(f"Selectivity by group results saved to {csv_paths['full']}")

    ## Functions 7 ##

    def update_pearson_plot_cmap(self):
        if not hasattr(self, "pearson_matrix"):
            return  # correlation not calculated yet

        mask = np.triu(np.ones_like(self.pearson_matrix, dtype=bool), k=1)

        # Select colormap
        selected_cmap = self.cmap_box.currentText()
        if hasattr(self, "custom_cmaps") and selected_cmap in self.custom_cmaps:
            cmap = self.custom_cmaps[selected_cmap]
            norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1) if selected_cmap in ["inferno_diverging",
                                                                                 "magma_diverging"] else None
        else:
            cmap = selected_cmap
            norm = None

        # --- Plot ---
        # Set default font globally (affects title, labels, tick labels, etc.)
        plt.rcParams.update({'font.family': 'Arial', 'font.size': 12})
        self.pearson_figure.clear()
        ax = self.pearson_figure.add_subplot(111)
        sns.heatmap(
            self.pearson_matrix,
            mask=mask,
            annot=True,
            fmt=".2f",
            cmap=cmap,
            norm=norm,
            ax=ax,
            cbar_kws={"shrink": 0.8}
        )

        # Title with Arial 12 explicitly
        ax.set_title("Pearson correlation map", fontsize=12, fontname='Arial')
        # Rotate x-axis labels
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=12, fontname='Arial')
        # Defined y-axis labels
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=12, fontname='Arial')

        plt.tight_layout()
        self.pearson_canvas.draw()

        # --- Automatically save the plot ---
        if hasattr(self, "plot_folder"):
            original_size = self.pearson_figure.get_size_inches()
            self.pearson_figure.set_size_inches(14, 8)
            pearson_photo_path = os.path.join(self.plot_folder, "pearson_plot.png")
            self.pearson_figure.savefig(pearson_photo_path, dpi=300, bbox_inches='tight')
            self.pearson_figure.set_size_inches(*original_size)

    def calculate_and_plot_pearson(self):

        df = self.latest_calibrated_df.copy()

        if self.latest_calibrated_df.empty:
            QMessageBox.warning(self, "Warning", "Empty dataset")
            return

        df["Sample Name"] = df["Sample Name"].astype(float)

        trend_data = {}

        # Loop through all retention times from self.retention_times_df
        for rt, comp in self.retention_times_df.itertuples(index=False):
            # Skip reactant and solvent peaks
            skip_rts = [self.reactant_rt] + [s["rt"] for s in self.solvents if s.get("rt") is not None]
            if any(abs(rt - srt) <= self.tolerance for srt in skip_rts):
                continue

            # Select rows within tolerance of this retention time
            tmp = df[df["Retention Time"].sub(rt).abs() <= self.tolerance]
            if tmp.empty:
                continue

            # Take mol%_normalized, indexed by Sample Name
            series = tmp.groupby("Sample Name")["mol%_normalized"].mean().sort_index()

            # Skip series that are all zeros or constant
            if (series != 0).sum() < 2 or series.nunique() == 1:
                continue

            trend_data[comp] = series

        if not trend_data:
            QMessageBox.warning(self, "Warning", "No valid compounds found for correlation.")
            return

        # Create DataFrame: columns=compounds, index=Sample Name
        df_trend = pd.DataFrame(trend_data)

        # Compute Pearson correlation
        self.pearson_matrix = df_trend.corr(method='pearson')  # store for cmap updates

        # Mask strictly upper triangle, keep diagonal
        mask = np.triu(np.ones_like(self.pearson_matrix, dtype=bool), k=1)

        # --- Select colormap (handle custom ones) ---
        selected_cmap = self.cmap_box.currentText()
        if hasattr(self, "custom_cmaps") and selected_cmap in self.custom_cmaps:
            cmap = self.custom_cmaps[selected_cmap]
            # Use TwoSlopeNorm for diverging maps
            if selected_cmap in ["inferno_diverging", "magma_diverging"]:
                norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
            else:
                norm = None
        else:
            cmap = selected_cmap
            norm = None


        # --- Plot ---
        # Set default font globally (affects title, labels, tick labels, etc.)
        plt.rcParams.update({'font.family': 'Arial', 'font.size': 12})

        self.pearson_figure.clear()
        ax = self.pearson_figure.add_subplot(111)
        sns.heatmap(
            self.pearson_matrix,
            mask=mask,
            annot=True,
            fmt=".2f",
            cmap=cmap,
            norm=norm,  # ensure diverging colormap is centered
            ax=ax,
            cbar_kws={"shrink": 0.8}
        )

        # Title with Arial 12 explicitly
        ax.set_title("Pearson correlation map", fontsize=12, fontname='Arial')
        # Rotate x-axis labels
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=12, fontname='Arial')
        # Defined y-axis labels
        ax.set_yticklabels(ax.get_yticklabels(), fontsize=12, fontname='Arial')

        plt.tight_layout()

        # Save graph 2 with temporary size
        original_size = self.pearson_figure.get_size_inches()
        self.pearson_figure.set_size_inches(14, 8)
        pearson_photo_path = os.path.join(self.plot_folder, "pearson_plot.png")
        self.pearson_figure.savefig(pearson_photo_path, dpi=300, bbox_inches='tight')
        self.pearson_figure.set_size_inches(*original_size)
        self.pearson_canvas.draw()

    ## Functions 8 ##

    def plot_carbon_recovery(self):
        # Copy dataframe
        df = self.latest_calibrated_df.copy()

        if df.empty:
            QMessageBox.warning(self, "Warning", "Empty dataset")
            return

        # Ensure Sample Name is numeric for proper sorting
        df['Sample Name'] = pd.to_numeric(df['Sample Name'], errors='coerce')

        # Summarize by Sample Name (remove repetitions) and sort
        carbon_df = df[['Sample Name', 'Carbon recovery (%)']].drop_duplicates(subset='Sample Name')
        carbon_df = carbon_df.sort_values('Sample Name')

        # Save summarized CSV
        csv_path = os.path.join(self.csv_folder, "carbon_recovery.csv")
        carbon_df.to_csv(csv_path, index=False)

        # --- Plotting ---
        self.carbon_fig.clear()
        ax = self.carbon_fig.add_subplot(111)

        # Plot line in black with fully purple dots
        ax.plot(
            carbon_df['Sample Name'],
            carbon_df['Carbon recovery (%)'],
            color='purple',
            linestyle='--',
            marker='o',
            markerfacecolor='purple',  # fill
            markeredgecolor='purple',  # edge
            markeredgewidth=1.5,
            markersize=6,
            label='Carbon Recovery'
        )

        # Set labels, title, grid, and legend
        ax.set_xlabel('Time on stream (h)', fontname='Arial', fontsize=12)
        ax.set_ylabel('Carbon recovery (%)', fontname='Arial', fontsize=12)
        ax.set_title('Carbon recovery over time', fontname='Arial', fontsize=12)
        ax.grid(True)
        ax.legend(prop={'family': 'Arial', 'size': 12})

        # Draw canvas
        self.carbon_canvas.draw()

        # --- Save figure with temporary size ---
        original_size = self.carbon_fig.get_size_inches()
        self.carbon_fig.set_size_inches(8, 6)
        fig_path = os.path.join(self.plot_folder, "carbon_plot.png")
        self.carbon_fig.savefig(fig_path, dpi=300, bbox_inches='tight')  # tight ensures labels aren't cut
        self.carbon_fig.set_size_inches(*original_size)

        # Update label
        self.carbon_plot_label.setText(f"Plot saved to: {fig_path}\nCSV saved to: {csv_path}")


def main():
    app = QApplication(sys.argv)
    window = RetentionTimeFilterApp()  # This class now includes a QTabWidget
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()