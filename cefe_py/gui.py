import sys
import os
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QPushButton, QLabel, QRadioButton, QGroupBox, QButtonGroup, 
    QFrame, QSizePolicy
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QPalette, QColor

import cefe_py as ce
from .render_engine import CefeRenderEngine

class CefeMainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("CEFE High-Fidelity Scientific Simulator")
        self.setGeometry(100, 100, 1200, 800)
        
        # Apply dark aesthetic
        self.apply_dark_theme()
        
        # Main layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QHBoxLayout(central_widget)
        main_layout.setContentsMargins(10, 10, 10, 10)
        
        # Sidebar for controls
        sidebar = QFrame()
        sidebar.setFixedWidth(300)
        sidebar.setStyleSheet("""
            QFrame {
                background-color: #1A1A1D;
                border-right: 1px solid #333;
            }
            QLabel {
                font-weight: 800;
                font-size: 18px;
                color: #FFFFFF;
                letter-spacing: 1px;
            }
            QGroupBox {
                font-weight: 700;
                color: #00E5FF;
                border: 1px solid #333;
                border-radius: 8px;
                margin-top: 15px;
                padding-top: 20px;
                font-size: 14px;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                subcontrol-position: top left;
                padding: 0 10px;
                top: -8px;
            }
            QRadioButton {
                color: #CCCCCC;
                font-size: 14px;
                padding: 5px;
            }
            QRadioButton::indicator {
                width: 16px;
                height: 16px;
                border-radius: 8px;
                border: 2px solid #555;
            }
            QRadioButton::indicator:checked {
                background-color: #00E5FF;
                border: 2px solid #00E5FF;
            }
            QPushButton {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #0052D4, stop:0.5 #4364F7, stop:1 #6FB1FC);
                color: white;
                border-radius: 8px;
                padding: 12px;
                font-weight: 800;
                font-size: 15px;
                border: none;
            }
            QPushButton:hover {
                background: qlineargradient(x1:0, y1:0, x2:1, y2:0, stop:0 #4364F7, stop:1 #6FB1FC);
            }
            QPushButton:disabled {
                background: #333;
                color: #666;
            }
        """)
        sidebar_layout = QVBoxLayout(sidebar)
        
        # Title
        title = QLabel("CEFE Control Panel")
        title.setAlignment(Qt.AlignmentFlag.AlignCenter)
        sidebar_layout.addWidget(title)
        
        # Mode Group
        mode_group = QGroupBox("Engine Mode")
        mode_layout = QVBoxLayout()
        self.mode_qm = QRadioButton("QM Double-Slit (Schrödinger)")
        self.mode_qm.setChecked(True)
        self.mode_qft = QRadioButton("Interacting QFT (Non-Linear)")
        self.mode_qft.setEnabled(True)
        self.mode_heg = QRadioButton("Holographic Gravity")
        
        self.mode_btn_group = QButtonGroup()
        self.mode_btn_group.addButton(self.mode_qm, 1)
        self.mode_btn_group.addButton(self.mode_qft, 2)
        self.mode_btn_group.addButton(self.mode_heg, 3)
        
        mode_layout.addWidget(self.mode_qm)
        mode_layout.addWidget(self.mode_qft)
        mode_layout.addWidget(self.mode_heg)
        mode_group.setLayout(mode_layout)
        sidebar_layout.addWidget(mode_group)
        
        # Action Buttons
        self.start_btn = QPushButton("Initialize Engine")
        self.start_btn.clicked.connect(self.init_engine)
        sidebar_layout.addWidget(self.start_btn)
        
        self.play_btn = QPushButton("▶ Run Simulation")
        self.play_btn.clicked.connect(self.toggle_play)
        self.play_btn.setEnabled(False)
        sidebar_layout.addWidget(self.play_btn)
        
        sidebar_layout.addStretch()
        
        # Simulation Viewport
        self.viewport = QFrame()
        self.viewport.setStyleSheet("background-color: #0F0F12; border-radius: 8px; border: 1px solid #333;")
        self.viewport_layout = QVBoxLayout(self.viewport)
        self.viewport_layout.setContentsMargins(5, 5, 5, 5)
        
        # Add a placeholder label
        self.placeholder = QLabel("Select a mode and initialize engine to view.")
        self.placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.placeholder.setStyleSheet("color: #777; font-size: 18px;")
        self.viewport_layout.addWidget(self.placeholder)
        
        # Assemble
        main_layout.addWidget(sidebar)
        main_layout.addWidget(self.viewport)
        
        self.renderer = None

    def apply_dark_theme(self):
        palette = QPalette()
        palette.setColor(QPalette.ColorRole.Window, QColor(15, 15, 18))
        palette.setColor(QPalette.ColorRole.WindowText, Qt.GlobalColor.white)
        palette.setColor(QPalette.ColorRole.Base, QColor(25, 25, 30))
        palette.setColor(QPalette.ColorRole.AlternateBase, QColor(15, 15, 18))
        palette.setColor(QPalette.ColorRole.ToolTipBase, Qt.GlobalColor.white)
        palette.setColor(QPalette.ColorRole.ToolTipText, Qt.GlobalColor.white)
        palette.setColor(QPalette.ColorRole.Text, Qt.GlobalColor.white)
        palette.setColor(QPalette.ColorRole.Button, QColor(40, 40, 45))
        palette.setColor(QPalette.ColorRole.ButtonText, Qt.GlobalColor.white)
        palette.setColor(QPalette.ColorRole.BrightText, Qt.GlobalColor.red)
        palette.setColor(QPalette.ColorRole.Link, QColor(0, 229, 255))
        palette.setColor(QPalette.ColorRole.Highlight, QColor(0, 229, 255))
        palette.setColor(QPalette.ColorRole.HighlightedText, Qt.GlobalColor.black)
        self.setPalette(palette)

    def init_engine(self):
        # Clear viewport
        for i in reversed(range(self.viewport_layout.count())): 
            self.viewport_layout.itemAt(i).widget().setParent(None)
            
        mode_id = self.mode_btn_group.checkedId()
        
        if mode_id == 1:
            print("Initializing QM C++ Engine...")
            grid = ce.CausalDiamondGrid(radius=2.0, spacing=0.1)
            engine = ce.QuantumWaveEngine(grid)
            engine.set_mass(1.0)
            engine.set_hbar(1.0)
            engine.initialize_state(target_t=0.0)
            engine.set_gaussian_packet(x0=-1.0, y0=0.0, sigma=0.2, px=5.0, py=0.0)
            engine.set_double_slit_potential(slit_width=0.4, slit_separation=1.0, barrier_thickness=0.2, barrier_x=0.0)
            engine.build_hamiltonian(target_t=0.0)
            
            self.renderer = CefeRenderEngine(engine, mode="qm")
            self.viewport_layout.addWidget(self.renderer)
            
            self.play_btn.setEnabled(True)
            self.play_btn.setText("▶ Run Simulation")
            
        elif mode_id == 2:
            print("Initializing QFT Interacting C++ Engine...")
            grid = ce.CausalDiamondGrid(radius=2.0, spacing=0.1)
            engine = ce.LatticeQFTEngine(grid)
            engine.set_mass(1.0)
            engine.set_coupling(10.0) # non-linear coupling
            engine.initialize_state()
            engine.set_opposing_packets(1.0, 0.2, 5.0, 5.0)
            engine.build_laplacian()
            
            self.renderer = CefeRenderEngine(engine, mode="qft")
            self.viewport_layout.addWidget(self.renderer)
            
            self.play_btn.setEnabled(True)
            self.play_btn.setText("▶ Run Simulation")
            
        elif mode_id == 3:
            print("Initializing Holographic Emergent Gravity Engine...")
            radius = 3.0
            spacing = 0.5
            ads_metric = ce.AntiDeSitterMetric(L=10.0)
            grid = ce.CausalDiamondGrid(radius=radius, spacing=spacing, metric=ads_metric)
            
            engine = ce.EntropicFieldEngine(grid)
            engine.set_mass(1.0)
            engine.initialize_field_state(target_t=0.0, temperature=5.0)
            
            self.renderer = CefeRenderEngine(engine, mode="heg")
            self.viewport_layout.addWidget(self.renderer)
            
            self.play_btn.setEnabled(True)
            self.play_btn.setText("▶ Run Simulation")

    def toggle_play(self):
        if self.renderer:
            if self.renderer.is_running:
                self.renderer.stop_animation()
                self.play_btn.setText("▶ Resume Simulation")
            else:
                self.renderer.start_animation()
                self.play_btn.setText("⏸ Pause Simulation")

def main():
    app = QApplication(sys.argv)
    window = CefeMainWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()
