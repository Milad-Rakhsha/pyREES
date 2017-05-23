from OpenGL.GL import *
import weakref
import numpy as np
import xml.etree.ElementTree as ET
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QOpenGLWidget
from PyQt5.QtCore import QTimer
from PyQt5.QtGui import QSurfaceFormat
from PyQt5.QtGui import QOpenGLVersionProfile
from PyQt5.QtCore import Qt
import REESSimulation.api as API
import REESSimulation.xml as XML
import REESUtility.util as UTIL
import REESVisualization.xml as VISXML
from REESVisualization.types import *
import REESMesh.factory as MESH_FACTORY


OpenGL.ERROR_CHECKING = True
OpenGL.CHECK_CONTEXT = True


class RenderWidget(QOpenGLWidget):

    def __init__(self):
        super(QOpenGLWidget, self).__init__()
        self.resize(400, 400)
        self.m_update_timer = QTimer()
        self.m_update_timer.timeout.connect(self.simulation_step)

        self.gl = 0

        self.camera = None
        self.light_setup = None
        self.grid = None
        self.clear_color = None
        self.materials = None
        self.scene_graphs = []
        self.movie_recorder = None
        self.contact_scale = 0.5

        self.keep_visualization = False  # Controls if visualization is loaded whening opening an XML file.

        self.engine = None

        self.trackball = Trackball()   # A track ball used to convert mouse move events into rotations
        self.dolly_mode = False        # If on then up-down mouse moves make one moves back and forth towards target
        self.pan_mode = False          # If on then one translates in the screen space
        self.trackball_mode = False    # If on then one rotates around the camera center
        self.fpv_mode = False          # If on then one rotates around the camera position
        self.selection_mode = False    # If on then one is making a selection of an object
        self.dolly_sensitivity = 0.025
        self.pan_sensitivity = 0.025
        self.anchor_x = None            # Original x-position when doing a mouse operation
        self.anchor_y = None            # Original t-position when doing a mouse operation
        self.anchor_eye = None
        self.anchor_center = None
        self.anchor_up = None
        self.height = None              # Used to convert from window space into screen space when clicking inside the window


    def compute_normalized_device_coordinates(self, sx, sy):
        '''

        :param sx:     X-coordinate of screen position (x-pixel value)
        :param sy:     Y-coordinate of screen position (y-pixel value)
        :return:
        '''
        viewport = glGetFloatv(GL_VIEWPORT)
        ratio = self.devicePixelRatio()
        nx = (2.0 * ratio * sx) / viewport[2] - 1.0
        ny = (2.0 * ratio * sy) / viewport[3] - 1.0
        nx = max(-1.0, min(1.0, nx))
        ny = max(-1.0, min(1.0, ny))
        return nx, ny

    def initializeGL(self):
        version_profile = QOpenGLVersionProfile()
        version_profile.setVersion(4, 1)
        version_profile.setProfile(QSurfaceFormat.CoreProfile)
        self.gl = self.context().versionFunctions(version_profile)
        if not self.gl:
            raise RuntimeError("unable to apply OpenGL version profile")
        self.gl.initializeOpenGLFunctions()

        camera, light_setup, materials, clear_color, grid, mov_rec, contact_scale = VISXML.load('resources/scene.xml')

        self.camera = camera
        self.light_setup = light_setup
        self.materials = materials
        self.clear_color = clear_color
        self.grid = grid
        self.movie_recorder = mov_rec
        self.contact_scale = contact_scale

        self.scene_graphs.clear()
        self.create_grid_scene_graph()

        glClearColor(self.clear_color[0], self.clear_color[1], self.clear_color[2], 0.0)
        glEnable(GL_DEPTH_TEST)
        glDepthFunc(GL_LESS)
        glEnable(GL_CULL_FACE)

    def create_rigid_bodies_graph(self):

        scene_graph = SceneGraph()

        scene_graph.camera = self.camera
        scene_graph.light_setup = self.light_setup

        vertex_shader = Shader('resources/shaders/vertex.glsl', GL_VERTEX_SHADER)
        geometry_shader = Shader('resources/shaders/geometry.glsl', GL_GEOMETRY_SHADER)
        fragment_shader = Shader('resources/shaders/fragment.glsl', GL_FRAGMENT_SHADER)

        scene_graph.program = ShaderProgram([vertex_shader, geometry_shader, fragment_shader])

        for shape in self.engine.shapes.values():
            vao = VAO()
            vertex_data, index_data = create_mesh_array_data(shape.mesh)
            vbo = VBO(vertex_data, index_data)
            shape_node = ShapeNode(scene_graph.program, vao, vbo)
            scene_graph.shape_nodes[shape.name] = shape_node

        for body in self.engine.rigid_bodies.values():
            shape_node = scene_graph.shape_nodes[body.shape.name]
            instance = InstanceNode(shape_node, self.materials[body.visual_material], body)
            instance.update_transform(body.r, body.q)
            scene_graph.instance_nodes[body.name] = instance

        self.scene_graphs.append(scene_graph)

    def create_grid_scene_graph(self):
        scene_graph = SceneGraph()

        scene_graph.camera = self.camera
        scene_graph.light_setup = self.light_setup

        vertex_shader = Shader('resources/shaders/debug_vertex.glsl', GL_VERTEX_SHADER)
        fragment_shader = Shader('resources/shaders/debug_fragment.glsl', GL_FRAGMENT_SHADER)

        scene_graph.program = ShaderProgram([vertex_shader, fragment_shader])

        vao = VAO()
        vertex_data, index_data = create_wire_grid_data(self.grid)
        vbo = VBO(vertex_data, index_data, make_triangles=False)
        shape_node = ShapeNode(scene_graph.program, vao, vbo)
        instance = InstanceNode(shape_node, self.materials['grid'])
        scene_graph.instance_nodes['grid'] = instance

        self.scene_graphs.append(scene_graph)

    def create_contact_points_graph(self):
        scene_graph = SceneGraph()

        scene_graph.camera = self.camera
        scene_graph.light_setup = self.light_setup

        vertex_shader = Shader('resources/shaders/debug_vertex.glsl', GL_VERTEX_SHADER)
        fragment_shader = Shader('resources/shaders/debug_fragment.glsl', GL_FRAGMENT_SHADER)

        scene_graph.program = ShaderProgram([vertex_shader, fragment_shader])

        vao = VAO()
        slices = 12
        segments = 12

        H = 1.0/(1.0 - 1.0/8.0)

        base = MESH_FACTORY.make_sphere(H/8.0,slices, segments)
        shaft = MESH_FACTORY.make_cylinder(H/16.0, 5.0*H/8.0, slices)
        head = MESH_FACTORY.make_cone(H/8.0, 2.0*H/8.0, slices)
        MESH.translate(shaft, V3.make(0.0, 5.0*H/16.0, 0.0))
        MESH.translate(head, V3.make(0.0, 5.0*H/8.0, 0.0))

        mesh = MESH.join([base, shaft, head])

        vertex_data, index_data = create_mesh_array_data(mesh)
        vbo = VBO(vertex_data, index_data)
        shape_node = ShapeNode(scene_graph.program, vao, vbo)
        instance = ContactPointsInstanceNode(shape_node, self.materials['contact'], self.contact_scale, self.engine)
        scene_graph.instance_nodes['contact'] = instance

        self.scene_graphs.append(scene_graph)

    def create_joint_connectors_graph(self):
        scene_graph = SceneGraph()

        scene_graph.camera = self.camera
        scene_graph.light_setup = self.light_setup

        vertex_shader = Shader('resources/shaders/debug_vertex.glsl', GL_VERTEX_SHADER)
        fragment_shader = Shader('resources/shaders/debug_fragment.glsl', GL_FRAGMENT_SHADER)

        scene_graph.program = ShaderProgram([vertex_shader, fragment_shader])

        vao = VAO()
        slices = 12
        segments = 12

        H = 1.0/(1.0 - 1.0/8.0)

        base = MESH_FACTORY.make_sphere(H/8.0,slices, segments)
        shaft = MESH_FACTORY.make_cylinder(H/16.0, 5.0*H/8.0, slices)
        head = MESH_FACTORY.make_cone(H/8.0, 2.0*H/8.0, slices)
        MESH.translate(shaft, V3.make(0.0, 5.0*H/16.0, 0.0))
        MESH.translate(head, V3.make(0.0, 5.0*H/8.0, 0.0))

        mesh = MESH.join([base, shaft, head])

        vertex_data, index_data = create_mesh_array_data(mesh)
        vbo = VBO(vertex_data, index_data)
        shape_node = ShapeNode(scene_graph.program, vao, vbo)

        instance = JointConnectorsInstanceNode(shape_node,
                                               self.materials['socket'],
                                               self.materials['plug'],
                                               self.contact_scale,
                                               self.engine
                                               )
        scene_graph.instance_nodes['connectors'] = instance

        self.scene_graphs.append(scene_graph)

    def paintGL(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        for scene_graph in self.scene_graphs:
            scene_graph.render()

    def resizeGL(self, width, height):
        glViewport(0, 0, width, height)
        self.height = height

    def mouseMoveEvent(self, e):
        x = e.x()
        y = self.height - e.y()

        nx, ny = self.compute_normalized_device_coordinates(x, y)

        left = (int(e.buttons()) & Qt.LeftButton) != 0
        middle = (int(e.buttons()) & Qt.MidButton) != 0
        right = (int(e.buttons()) & Qt.RightButton) != 0
        ctrl = (int(QApplication.keyboardModifiers()) & Qt.ControlModifier) != 0
        shift = (int(QApplication.keyboardModifiers()) & Qt.ShiftModifier) != 0
        alt = (int(QApplication.keyboardModifiers()) & Qt.AltModifier) != 0

        self.camera.update(self.anchor_eye, self.anchor_center, self.anchor_up)

        if self.dolly_mode:
            distance = self.dolly_sensitivity * (y - self.anchor_y)
            self.camera.dolly(-distance)

        if self.pan_mode:
            x_distance = self.pan_sensitivity * (x - self.anchor_x)
            y_distance = self.pan_sensitivity * (y - self.anchor_y)
            self.camera.pan(-x_distance, -y_distance)

        if self.trackball_mode:
            self.trackball.move_to(nx, ny)
#            self.trackball.move_to(x, y)
            self.camera.orbit(self.trackball.rotation_matrix.transpose())

        if self.fpv_mode:
            self.trackball.move_to(nx, ny)
            self.camera.rotate(self.trackball.rotation_matrix)

        if self.selection_mode:
            p, r = self.camera.get_ray(nx, ny)
            #self.select_tool.move_selection(p, r, self.camera.dof, self.engine)
            #update_render_manager(self.render_manager, self.engine)

        if not self.m_update_timer.isActive():
            self.update()

    def mousePressEvent(self, e):
        x = e.x()
        y = self.height - e.y()

        nx, ny = self.compute_normalized_device_coordinates(x, y)

        left = (int(e.buttons()) & Qt.LeftButton) != 0
        middle = (int(e.buttons()) & Qt.MidButton) != 0
        right = (int(e.buttons()) & Qt.RightButton) != 0
        ctrl = (int(QApplication.keyboardModifiers()) & Qt.ControlModifier) != 0
        shift = (int(QApplication.keyboardModifiers()) & Qt.ShiftModifier) != 0
        alt = (int(QApplication.keyboardModifiers()) & Qt.AltModifier) != 0

        if alt and left:
            self.dolly_mode = True
        elif shift and left:
            self.pan_mode = True
            self.camera.center_locked = False
        elif ctrl and left:
            self.selection_mode = True
        elif left:
            self.trackball_mode = True
        elif right:
            self.fpv_mode = True
            self.camera.center_locked = False

        self.trackball.reset()
        self.anchor_x = x
        self.anchor_y = y
        self.anchor_eye = np.copy(self.camera.eye)
        self.anchor_center = np.copy(self.camera.center)
        self.anchor_up = np.copy(self.camera.up)

        if self.trackball_mode:
            self.trackball.click_at(nx, ny)

        if self.fpv_mode:
            self.trackball.click_at(nx, ny)

        if self.selection_mode:
            p,r = self.camera.get_ray(nx, ny)
            #self.select_tool.select(p, r, self.engine )
            self.selection_mode = True

    def mouseReleaseEvent(self, e):
        #if self.selection_mode:
        #   self.select_tool.deselect()
        self.dolly_mode = False
        self.pan_mode = False
        self.trackball_mode = False
        self.fpv_mode = False
        self.camera.center_locked = True
        self.selection_mode = False

    def simulation_step(self):
        did_step = API.run(self.engine)

        if did_step and self.movie_recorder is not None:
            self.movie_recorder.record(self)

        self.update()

    def open_file(self, filename):
        xml = ET.parse(filename)
        root = xml.getroot()

        if not self.keep_visualization:
            camera, light_setup, materials, clear_color, grid, mov_rec, contact_scale = VISXML.load_from_elementtree(root)

            self.camera = camera
            self.light_setup = light_setup
            self.materials = materials
            self.clear_color = clear_color
            self.grid = grid
            self.movie_recorder = mov_rec
            self.contact_scale = contact_scale

        self.scene_graphs.clear()
        self.create_grid_scene_graph()

        self.engine = XML.load_from_elementtree(root)
        self.create_rigid_bodies_graph()
        self.create_contact_points_graph()
        self.create_joint_connectors_graph()

        glClearColor(self.clear_color[0], self.clear_color[1], self.clear_color[2], 0.0)

        self.m_update_timer.start(1000/self.engine.solver_params.fps)

        self.update()

    def save_file(self, filename):
        root = ET.Element('scene')

        VISXML.save_to_elementtree(self.camera,
                                   self.light_setup,
                                   self.materials,
                                   self.clear_color,
                                   self.grid,
                                   self.movie_recorder,
                                   self.contact_scale,
                                   root
                                   )

        XML.save_to_elementtree(self.engine, root)

        UTIL.xml_pretty_indent(root)
        tree = ET.ElementTree(root)
        tree.write(filename)

    def set_keep_visualization(self, keep_it):
        self.keep_visualization = keep_it
        if keep_it:
            print('RenderWidget: Keeping visualization on file open')
        else:
            print('RenderWidget: Loading visualization on file open')


def initialize_opengl():
    format = QSurfaceFormat()
    format.setProfile(QSurfaceFormat.CoreProfile)
    format.setVersion(4, 5)
    format.setSamples(4)
    QSurfaceFormat.setDefaultFormat(format)

