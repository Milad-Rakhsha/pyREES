import networkx as nx
import REESMath.vector3 as V3
import REESMath.quaternion as Q
import REESMath.coordsys as C
import REESMath.matrix3 as M3


class MaterialBehaviour:
    def __init__(self):
        self.mu = V3.ones()   # Coefficients of Friction
        self.epsilon = 0.0    # Coefficient of restitution


class MaterialLibrary:
    def __init__(self):
        self.storage = dict()
        self.storage[('default', 'default')] = MaterialBehaviour()

    def get_behaviour(self, A, B):
        key = (A, B) if A < B else (B, A)
        if key in self.storage:
            return self.storage[key]
        return self.storage[('default', 'default')]

    def exist_behaviour(self, A, B):
        key = (A, B) if A < B else (B, A)
        if key in self.storage:
            return True
        return False

    def exist_material(self, name):
        for key in self.storage:
            if name in key:
                return True
        return False


class KeyframeMotion:

    def __init__(self):
        self.keyframes = []

    def create_keyframe(self, time, r, q):
        self.keyframes.append([time, r, q])

    def clear(self):
        self.keyframes.clear()


class ForceCalculator:
    def __init__(self, force_type, name):
        self.force_type = force_type
        self.name = name


class Gravity(ForceCalculator):
    def __init__(self, name):
        super().__init__('Gravity', name)
        self.g = 9.81              # Acceleration of gravity
        self.up = V3.j()  # Up direction

    def compute(self, body, r, q, v, w):
        F = - body.mass * self.g * self.up
        T = V3.zero()
        return F, T


class Damping(ForceCalculator):
    def __init__(self, name):
        super().__init__('Damping', name)
        self.alpha = 0.001              # Linear damping
        self.beta = 0.001               # Angular damping

    def compute(self, body, r, q, v, w):
        F = - v * self.alpha
        T = - w * self.beta
        return F, T


class Shape:

    def __init__(self, name):
        self.name = name
        self.mesh = None             # Polygonal mesh assumed to be in body frame coordinates
        self.mass = 0.0              # Total mass of shape assuming unit-mass-density
        self.inertia = V3.zero()     # Body frame inertia tensor assumping unit-mass-density
        self.r = V3.zero()           # Translation from body frame to model frame
        self.q = Q.identity()        # Rotation from body frame to model frame


class RigidBody:

    def __init__(self, name):
        self.name = name
        self.q = Q.identity()       # Orientation stored as a quaternion
        self.r = V3.zero()          # Center of mass position
        self.v = V3.zero()          # Linear velocity
        self.w = V3.zero()          # Angular velocity
        self.mass = 0.0             # Total mass
        self.inertia = V3.zero()    # Body frame inertia tensor
        self.is_fixed = False
        self.is_scripted = False
        self.is_free = True
        self.is_active = True
        self.use_finite_update = True               # Toggle between infinitesimal and finite updates of orientation
        self.finite_update_rotation_axis = None     # If a vector3 is given then finite updates are done wrt. this axis
        self.forces = []                            # External forces (like gravity and damping) acting on this body
        self.material = 'default'                   # The material this rigid body is made up of
        self.shape = None                           # Geometry/Shape of rigid body
        self.kdop = None                            # World space kdop bvh
        self.joints = []                            # Joints connected to this body
        self.scripted_motion = None
        self.visual_material = None                 # string value that makes a reference to the name of the visual material to use when rendering this body.


class Constraint:

    def __init__(self, txt, bodyA, bodyB):
        self.type = txt
        self.bodyA = bodyA
        self.bodyB = bodyB


class ContactPoint(Constraint):

    def __init__(self, bodyA, bodyB, position=V3.zero(), normal=V3.k(), gap=0.0):
        super().__init__('ContactPoint', bodyA, bodyB)

        if abs(1.0 - V3.norm(normal)) > 0.1:
            raise RuntimeError('ContactPoint.init() was called with non-unit size normal')

        self.p = position
        self.n = normal
        self.g = gap


class JointConnector:

    def __init__(self, body):
        self.body = body
        self.transform = C.CoordSys()

    def get_local_anchor(self):
        return self.transform.r

    def get_local_axis_1(self):
        return Q.rotate(self.transform.q, V3.i())

    def get_local_axis_2(self):
        return Q.rotate(self.transform.q, V3.j())

    def get_local_axis_3(self):
        return Q.rotate(self.transform.q, V3.k())

    def get_world_anchor(self):
        return Q.rotate(self.body.q, self.transform.r) + self.body.r

    def get_world_axis_1(self):
        return Q.rotate(self.transform.q, self.get_local_axis_1())

    def get_world_axis_2(self):
        return Q.rotate(self.transform.q, self.get_local_axis_2())

    def get_world_axis_3(self):
        return Q.rotate(self.transform.q, self.get_local_axis_3())

    def get_world_arm(self):
        return Q.rotate(self.body.q, self.transform.r)


class Joint(Constraint):

    def __init__(self, type, name, socket_body, plug_body):
        super().__init__(type, socket_body, plug_body)
        self.name = name
        self.socket = JointConnector(socket_body)
        self.plug = JointConnector(plug_body)
        self.error_reduction = 0.1


class BallJoint(Joint):

    def __init__(self, name, socket_body, plug_body):
        super().__init__('ball', name, socket_body, plug_body)

    def compute_jacobians(self):
        JA_v = M3.identity()
        JA_w = - M3.star(self.socket.get_world_arm())
        JB_v = -M3.identity()
        JB_w = M3.star(self.plug.get_world_arm())
        return JA_v, JA_w, JB_v, JB_w

    def compute_error_terms(self, fps):

        if fps < 0.0:
            raise RuntimeError('compute_error_terms() Illegal fps value')

        if fps > 200.0:
            raise RuntimeWarning('compute_error_terms() Unlikely large fps value')

        k_correction = fps * self.error_reduction
        b = k_correction * (self.plug.get_world_anchor() - self.socket.get_world_anchor())
        return b


class GraphEdge:
    def __init__(self):
        self.contacts = []


class MotionRecorder:

    def __init__(self):
        self.on = False
        self.load = False                 # Boolean flag telling if channel data should be loaded when opening xml file
        self.save = False                 # Boolean flag telling if channel data should be saved when savj g xml file
        self.path = './'
        self.filename = 'motion.xml'
        self.storage = dict()             # Maps body to keyframe motion

    def record(self, time, body):
        if body is None:
            return

        if body.name not in self.storage:
            self.storage[body.name] = KeyframeMotion()

        self.storage[body.name].create_keyframe(time, body.r, body.q)

    def clear(self):
        self.storage.clear()

    def set_state(self, is_on):
        self.on = is_on
        if is_on:
            print('Motion recorder is on')
        else:
            print('Motion recorder is off')


class Profiler:

    def __init__(self):
        self.on = False
        self.path = './'
        self.filename = 'profiling.xml'

    def set_state(self, is_on):
        self.on = is_on
        if is_on:
            print('Profiler is on')
        else:
            print('Profiler is off')


class SolverParameters:

    def __init__(self):
        self.total_time = 10.0        # The total allowed simulation time.
        self.current_time = 0.0       # The current simulation time.
        self.time_step = 0.001        # The time step size to use when taking one simulation solver step. Observe that there could be many steps for a single frame update.
        self.fps = 60.0               # Number of frames to be made per second. This parameters helps control movie recording.
        self.on = False               # Boolean flag that turns simulation on and off.
        self.mode = 'simulate'        # The current mode of the solver, can be play (of recorded motion) or simulate.
        self.load_restart = False     # If true then body states will be initialized with restart file
        self.save_restart = False     # If true then simulator will save restart data after each simulation step.
        self.restart_path = './'
        self.restart_filename = 'restart.xml'

    def set_state(self, is_on):
        self.on = is_on
        if is_on:
            print('Simulation solver is on')
        else:
            print('Simulation solver is off')


class Engine:

    def __init__(self):
        self.graph = nx.Graph()

        self.rigid_bodies = dict()
        self.forces = dict()
        self.shapes = dict()

        self.contact_points = []
        self.joints = dict()
        self.joint_limits = []
        self.joint_motors = []

        self.material_library = MaterialLibrary()
        self.solver_params = SolverParameters()
        self.profiler = Profiler()
        self.motion_recorder = MotionRecorder()

