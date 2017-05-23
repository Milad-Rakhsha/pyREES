import sys
import os
from PyQt5.QtWidgets import QWidget
from PyQt5.QtWidgets import QToolTip
from PyQt5.QtWidgets import QPushButton
from PyQt5.QtWidgets import QCheckBox
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtWidgets import QDesktopWidget
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QAction
from PyQt5.QtWidgets import qApp
from PyQt5.QtWidgets import QTextEdit
from PyQt5.QtWidgets import QTreeView
from PyQt5.QtWidgets import QVBoxLayout
from PyQt5.QtWidgets import QDockWidget
from PyQt5.QtWidgets import QMenu
from PyQt5.QtWidgets import QFileDialog
from PyQt5.QtGui import QFont, QIcon, QStandardItem, QStandardItemModel, QPixmap
from PyQt5.QtCore import QCoreApplication, QSize
from PyQt5.Qt import  Qt
from REESVisualization.widget import initialize_opengl, RenderWidget


#my_data = [
#    ("Alice", [
#        ("Keys", []),
#        ("Purse", [
#            ("Cellphone", [])
#        ])
#    ]),
#    ("Bob", [
#        ("Wallet", [
#            ("Credit card", []),
#            ("Money", [])
#        ])
#    ])
#]


class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()

        self.initialize_user_interface()

    def initialize_user_interface(self):
        self.setAcceptDrops(True)
        self.create_render_widget()
        self.create_icons()
        self.create_actions()
        self.create_main_menu()
        self.create_toolbar()
        self.create_layout()
        self.create_statusbar()
        self.show()

    def create_render_widget(self):
        self.render_widget = RenderWidget()
        self.render_widget.setMinimumSize(1280, 720)  # We go for a HD 720p60 format

    def create_statusbar(self):
        self.statusBar().showMessage('Ready')

    def create_icons(self):
        self.exit_icon = QIcon('resources/images/exit_icon.png')
        self.profiling_icon = QIcon('resources/images/profile_icon.png')
        self.record_motion_icon = QIcon('resources/images/record_motion_icon.png')
        self.record_movie_icon = QIcon('resources/images/record_movie_icon.png')
        self.run_simulation_icon = QIcon('resources/images/run_simulation_icon.png')
        self.open_file_icon = QIcon('resources/images/open_file_icon.png')
        self.save_file_icon = QIcon('resources/images/save_file_icon.png')
        self.keep_visualization_icon = QIcon('resources/images/on_icon.png')

    def create_actions(self):
        self.open_file_action = QAction(self.open_file_icon, 'Open File', self)
        self.open_file_action.setShortcut('Ctrl+F')
        self.open_file_action.setStatusTip('Open File')
        self.open_file_action.triggered.connect(
            lambda: self.open_file_dialog()
        )

        self.save_file_action = QAction(self.save_file_icon, 'Save File', self)
        self.save_file_action.setShortcut('Ctrl+S')
        self.save_file_action.setStatusTip('Save File')
        self.save_file_action.triggered.connect(
            lambda: self.save_file_dialog()
        )

        self.exit_action = QAction(self.exit_icon, 'Exit', self)
        self.exit_action.setShortcut('Ctrl+Q')
        self.exit_action.setStatusTip('Exit application')
        self.exit_action.triggered.connect(self.close)

        self.run_simulation_action = QAction(self.run_simulation_icon, 'Run Simulation', self)
        self.run_simulation_action.setShortcut('Ctrl+P')
        self.run_simulation_action.setStatusTip('Run Simulation')
        self.run_simulation_action.setCheckable(True)
        self.run_simulation_action.setChecked(False)
        self.run_simulation_action.triggered.connect(
            lambda: self.render_widget.engine.solver_params.set_state(self.run_simulation_action.isChecked())
        )

        self.profiling_action = QAction(self.profiling_icon, 'Profiling', self)
        self.profiling_action.setShortcut('Ctrl+T')
        self.profiling_action.setStatusTip('Profiling')
        self.profiling_action.setCheckable(True)
        self.profiling_action.setChecked(False)
        self.profiling_action.triggered.connect(
            lambda: self.render_widget.engine.profiler.set_state(self.profiling_action.isChecked())
        )

        self.record_motion_action = QAction(self.record_motion_icon, 'Record Motion', self)
        self.record_motion_action.setShortcut('Ctrl+M')
        self.record_motion_action.setStatusTip('Record Motion')
        self.record_motion_action.setCheckable(True)
        self.record_motion_action.setChecked(False)
        self.record_motion_action.triggered.connect(
            lambda: self.render_widget.engine.motion_recorder.set_state(self.record_motion_action.isChecked())
        )

        self.record_movie_action = QAction(self.record_movie_icon, 'Record Movie', self)
        self.record_movie_action.setShortcut('Ctrl+M')
        self.record_movie_action.setStatusTip('Record Movie')
        self.record_movie_action.setCheckable(True)
        self.record_movie_action.setChecked(False)
        self.record_movie_action.triggered.connect(
            lambda: self.render_widget.movie_recorder.set_state(self.record_movie_action.isChecked())
        )

        self.keep_visualization_action = QAction(self.keep_visualization_icon, 'Keep Visualization', self)
        self.keep_visualization_action.setShortcut('Ctrl+V')
        self.keep_visualization_action.setStatusTip('Keep Visualization')
        self.keep_visualization_action.setCheckable(True)
        self.keep_visualization_action.setChecked(self.render_widget.keep_visualization)
        self.keep_visualization_action.triggered.connect(
            lambda: self.render_widget.set_keep_visualization(self.keep_visualization_action.isChecked())
        )

    def create_main_menu(self):
        self.menuBar().setNativeMenuBar(False)

        self.file_menu = self.menuBar().addMenu('File')
        self.file_menu.addAction(self.open_file_action)
        self.file_menu.addAction(self.save_file_action)
        self.file_menu.addAction(self.exit_action)

        self.actions_menu = self.menuBar().addMenu('Action')
        self.actions_menu.addAction(self.run_simulation_action)
        self.actions_menu.addAction(self.profiling_action)
        self.actions_menu.addAction(self.record_motion_action)
        self.actions_menu.addAction(self.record_movie_action)

        self.view_menu = self.menuBar().addMenu('View')
        self.view_menu.addAction(self.keep_visualization_action)

    def create_toolbar(self):
        self.toolbar = self.addToolBar('Main')
        self.toolbar.setIconSize(QSize(16, 16))

        self.toolbar.addAction(self.open_file_action)
        self.toolbar.addAction(self.save_file_action)
        self.toolbar.addSeparator()
        self.toolbar.addAction(self.run_simulation_action)
        self.toolbar.addAction(self.keep_visualization_action)
        self.toolbar.addSeparator()
        self.toolbar.addAction(self.profiling_action)
        self.toolbar.addAction(self.record_motion_action)
        self.toolbar.addAction(self.record_movie_action)
        self.toolbar.addSeparator()
        self.toolbar.addAction(self.exit_action)

        self.toolbar.setOrientation(Qt.Vertical)
        self.addToolBar(Qt.LeftToolBarArea, self.toolbar)

    def create_layout(self):
        dock_render_widget = QDockWidget()
        dock_render_widget.setWidget(self.render_widget)
        self.addDockWidget(Qt.LeftDockWidgetArea, dock_render_widget)

        dock_render_widget.toggleViewAction().setStatusTip('Toggle GL Viewer')
        dock_render_widget.toggleViewAction().setText('GL Viewer')
        dock_render_widget.toggleViewAction().setShortcut('Ctrl+G')
        self.view_menu.addAction(dock_render_widget.toggleViewAction())

        self.setGeometry(300, 300, 640, 480)
        self.setWindowTitle('Python Simulation')

    def dragEnterEvent(self, e):
        if e.mimeData().hasText() and e.mimeData().hasUrls:
            e.accept()
        else:
            e.ignore()

    def dropEvent(self, e):
        for url in e.mimeData().urls():
            path = str(url.toLocalFile())
            self.open_file(path)

    def open_file_dialog(self):
        fileinfo = QFileDialog.getOpenFileName(self, 'Open File')
        path = fileinfo[0]
        self.open_file(path)

    def open_file(self, path):
        ext = os.path.splitext(path)[-1].lower()
        if ext == '.xml' and os.path.isfile(path):
            self.render_widget.open_file(path)

        self.record_movie_action.setChecked(self.render_widget.movie_recorder.on)
        self.profiling_action.setChecked(self.render_widget.engine.profiler.on)
        self.record_motion_action.setChecked(self.render_widget.engine.motion_recorder.on)
        self.run_simulation_action.setChecked(self.render_widget.engine.solver_params.on)

    def save_file_dialog(self):
        fileinfo = QFileDialog.getSaveFileName(self, 'Save File')
        path = fileinfo[0]
        ext = os.path.splitext(path)[-1].lower()
        if ext == '.xml':
            self.render_widget.save_file(path)
