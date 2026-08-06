import numpy as np
from PyQt6.QtWidgets import QWidget, QVBoxLayout
from PyQt6.QtCore import QThread, pyqtSignal
import pyqtgraph.opengl as gl
import cefe_py as ce
import time

class EngineWorker(QThread):
    # Signal emitted when a new frame is ready to render
    frame_ready = pyqtSignal(object)

    def __init__(self, engine, mode, dt):
        super().__init__()
        self.engine = engine
        self.mode = mode
        self.dt = dt
        self.is_running = False

    def run(self):
        self.is_running = True
        while self.is_running:
            # The heavy C++ matrix multiplication happens here!
            # Since C++ releases the GIL, this won't block the Python GUI thread.
            
            if self.mode == "qm":
                self.engine.step_forward()
                prob = self.engine.get_probability_density()
            elif self.mode == "qft":
                self.engine.step_forward(self.dt)
                prob = self.engine.get_field_amplitude()
                prob = [p*p for p in prob] # Square for pseudo-density visibility
            elif self.mode == "heg" or self.mode == "holographic":
                self.engine.step_forward(self.dt)
                prob = self.engine.get_field_energy()
                
            self.frame_ready.emit(prob)
            
            # Rate limiting to avoid overwhelming the GUI thread queue
            time.sleep(0.016) # ~60 FPS cap

    def stop(self):
        self.is_running = False
        self.wait()

class CefeRenderEngine(QWidget):
    @property
    def is_running(self):
        return self.worker.is_running if hasattr(self, 'worker') else False

    def __init__(self, cefe_engine, mode="qm", parent=None):
        super().__init__(parent)
        self.engine = cefe_engine
        self.mode = mode
        
        # Initialize Layout
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 0, 0, 0)
        
        # Setup Pyqtgraph OpenGL View
        self.view = gl.GLViewWidget()
        self.layout.addWidget(self.view)
        
        # Configure Camera
        self.view.setCameraPosition(distance=40, elevation=30, azimuth=45)
        
        # Add Grid
        g = gl.GLGridItem()
        g.scale(2, 2, 2)
        self.view.addItem(g)
        
        # Initialize 3D Surface
        self.size = 50
        self.L = 2.0
        self.x_vals = np.linspace(-self.L, self.L, self.size)
        self.y_vals = np.linspace(-self.L, self.L, self.size)
        
        self.z_vals = np.zeros((self.size, self.size))
        
        self.surface = gl.GLSurfacePlotItem(
            x=self.x_vals, y=self.y_vals, z=self.z_vals, 
            shader='heightColor', computeNormals=True, smooth=True
        )
        
        self.surface.shader()['colorMap'] = np.array([
            0.0, 0.0, 0.0, 1.0,
            0.5, 0.0, 0.5, 1.0,
            1.0, 0.5, 0.0, 1.0,
            1.0, 1.0, 0.0, 1.0,
            1.0, 1.0, 1.0, 1.0
        ])
        
        self.view.addItem(self.surface)
        
        # Animation State
        self.dt = 0.01
        if self.mode == "qm":
            self.engine.prepare_crank_nicolson(self.dt)
            
        self.worker = EngineWorker(self.engine, self.mode, self.dt)
        self.worker.frame_ready.connect(self.update_frame)

    def start_animation(self):
        if not self.worker.is_running:
            self.worker.start()
            
    def stop_animation(self):
        if self.worker.is_running:
            self.worker.stop()

    def update_frame(self, prob):
        if self.mode == "heg" or self.mode == "holographic":
            energy = prob
            current_time = time.time() * 3.0
            X, Y = np.meshgrid(self.x_vals, self.y_vals)
            R = np.sqrt(X**2 + Y**2)
            Z = (np.log10(abs(energy) + 1.0)) * np.sin(5.0 * R - current_time) * np.exp(-0.3 * R) * 5.0
            self.surface.setData(x=self.x_vals, y=self.y_vals, z=Z)
            return

        total_points = len(prob)
        actual_size = int(np.sqrt(total_points))
        
        if actual_size > 0:
            Z_flat = np.array(prob[:actual_size*actual_size])
            Z = Z_flat.reshape((actual_size, actual_size))
            
            if self.mode == "qm":
                Z_scaled = Z * 15.0
            elif self.mode == "qft":
                Z_scaled = Z * 0.5
            else:
                Z_scaled = Z
            
            if actual_size != self.size:
                self.size = actual_size
                self.x_vals = np.linspace(-self.L, self.L, self.size)
                self.y_vals = np.linspace(-self.L, self.L, self.size)
                
            self.surface.setData(x=self.x_vals, y=self.y_vals, z=Z_scaled)

    def render_3d_wave_animation(self, steps=100, dt=0.01):
        self.dt = dt
        self.show()
        self.start_animation()
        
    def render_2d_heatmap(self, steps=100, dt=0.01):
        print("2D Heatmap rendering mapped to standard 3D viewport.")
        self.render_3d_wave_animation(steps, dt)
