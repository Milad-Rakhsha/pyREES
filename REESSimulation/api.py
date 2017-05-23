from REESSimulation.solver import *
from REESSimulation.types import *


def run(engine):
    if engine is None:
        return

    if engine.solver_params.mode == 'simulate':
        return simulate(engine)

    if engine.solver_params.mode == 'play':
        return play(engine)

    return False


def create_rigid_body(engine, body_name):
    if body_name in engine.rigid_bodies:
        raise RuntimeError('connect() rigid body already exist with that name')
    engine.rigid_bodies[body_name] = RigidBody(body_name)


def create_shape_from_obj(engine, shape_name, filename):
    if shape_name in engine.shapes:
        raise RuntimeError('create_shape_from_obj(): shape with that name already exist')
    shape = Shape(shape_name)
    shape.mesh = MESH.read_obj(filename)
    transform_shape_into_body_frame(shape)
    engine.shapes[shape_name] = shape


def create_cylinder_shape(engine, shape_name, radius, height):
    if shape_name in engine.shapes:
        raise RuntimeError('create_cylinder_shape(): shape with that name already exist')
    shape = Shape(shape_name)
    shape.mesh = MESH_FACTORY.make_cylinder(radius, height, slices=12)
    transform_shape_into_body_frame(shape)
    engine.shapes[shape_name] = shape


def create_convex_shape(engine, shape_name, points):
    if shape_name in engine.shapes:
        raise RuntimeError('create_convex_shape(): shape with that name already exist')
    shape = Shape(shape_name)
    shape.mesh = MESH_FACTORY.make_convex_hull(points)
    transform_shape_into_body_frame(shape)
    engine.shapes[shape_name] = shape


def create_capsule_shape(engine, shape_name, radius, height):
    if shape_name in engine.shapes:
        raise RuntimeError('create_capsule_shape(): shape with that name already exist')
    shape = Shape(shape_name)
    shape.mesh = MESH_FACTORY.make_capsule(radius, height, slices=12, segments=12)
    transform_shape_into_body_frame(shape)
    engine.shapes[shape_name] = shape


def create_cone_shape(engine, shape_name, radius, height):
    if shape_name in engine.shapes:
        raise RuntimeError('create_cone_shape(): shape with that name already exist')
    shape = Shape(shape_name)
    shape.mesh = MESH_FACTORY.make_cone(radius, height, slices=12)
    transform_shape_into_body_frame(shape)
    engine.shapes[shape_name] = shape


def create_conical_shape(engine, shape_name, bottom_radius, top_radius, height):
    if shape_name in engine.shapes:
        raise RuntimeError('create_conical_shape(): shape with that name already exist')
    shape = Shape(shape_name)
    shape.mesh = MESH_FACTORY.make_conical(bottom_radius, top_radius, height, slices=12)
    transform_shape_into_body_frame(shape)
    engine.shapes[shape_name] = shape


def create_sphere_shape(engine, shape_name, radius):
    if shape_name in engine.shapes:
        raise RuntimeError('create_sphere_shape(): shape with that name already exist')
    shape = Shape(shape_name)
    shape.mesh = MESH_FACTORY.make_sphere(radius, slices=12, segments=12)
    transform_shape_into_body_frame(shape)
    engine.shapes[shape_name] = shape


def create_ellipsoid_shape(engine, shape_name, a, b, c):
    if shape_name in engine.shapes:
        raise RuntimeError('create_ellipsoid_shape(): shape with that name already exist')
    shape = Shape(shape_name)
    shape.mesh = MESH_FACTORY.make_ellipsoid(a, b, c, slices=12, segments=12)
    transform_shape_into_body_frame(shape)
    engine.shapes[shape_name] = shape


def create_tetrahedron_shape(engine, shape_name, p0, p1, p2, p3):
    if shape_name in engine.shapes:
        raise RuntimeError('create_ellipsoid_shape(): shape with that name already exist')
    shape = Shape(shape_name)
    shape.mesh = MESH_FACTORY.make_tetrahedron(p0, p1, p2, p3)
    transform_shape_into_body_frame(shape)
    engine.shapes[shape_name] = shape


def create_cuboid_shape(engine, shape_name, p0, p1, p2, p3, p4, p5, p6, p7):
    if shape_name in engine.shapes:
        raise RuntimeError('create_cuboid_shape(): shape with that name already exist')
    shape = Shape(shape_name)
    shape.mesh = MESH_FACTORY.make_cuboid(p0, p1, p2, p3, p4, p5, p6, p7)
    transform_shape_into_body_frame(shape)
    engine.shapes[shape_name] = shape


def create_box_shape(engine, shape_name, width, height, depth):
    if shape_name in engine.shapes:
        raise RuntimeError('create_ellipsoid_shape(): shape with that name already exist')
    shape = Shape(shape_name)
    shape.mesh = MESH_FACTORY.make_box(width, height, depth)
    transform_shape_into_body_frame(shape)
    engine.shapes[shape_name] = shape


def connect_shape(engine, body_name, shape_name):
    if body_name in engine.rigid_bodies:
        body = engine.rigid_bodies[body_name]
    else:
        raise RuntimeError('connect() no such rigid body exist with that name')
    if shape_name in engine.shapes:
        shape = engine.shapes[shape_name]
    else:
        raise RuntimeError('connect() no such shape exist with that name')
    body.shape = shape


def set_position(engine, body_name, r, use_model_frame=False):
    if body_name in engine.rigid_bodies:
        body = engine.rigid_bodies[body_name]
    else:
        raise RuntimeError('set_position() no such rigid body exist with that name')
    if use_model_frame:

        r_bf2wcs = body.r
        q_bf2wcs = body.q

        r_bf2mf = body.shape.r
        q_bf2mf = body.shape.q
        #                                                                      |q r| |x|
        # By definition we have the rigid body transformations T x = q*x + r = |0 q| |1|
        #
        #  T_bf2wcs =   T_mf2wcs T_bf2mf
        #
        #   |q_bf2wcs  r_bf2wcs|   |q_mf2wcs  r_mf2wcs|  |q_bf2mf   r_bf2mf |
        #   |0            1    | = |0            1    |  |0            1    |
        #
        # So
        #
        #   q_bf2wcs = q_mf2wcs q_bf2mf
        #   r_bf2wcs = q_mf2wcs r_bf2mf + r_mf2wcs      (*)
        #
        # From this we can solve for the current model to world coordinate transformation
        q_mf2wcs = Q.prod(q_bf2wcs, Q.conjugate(q_bf2mf))
        r_mf2wcs = r_bf2wcs - Q.rotate(q_mf2wcs, r_bf2mf)
        #
        # Now we wish to change the position of model frame origin to (x,y,z) wrt the world frame
        #
        # so
        #
        #   r_mf2wcs = [x, y, z]
        #
        # The orientation of the model frame wrt. world frame is unchanged, hence we want to
        # know what r_bf2wcs should be. For this we use (*) to compute it
        #
        r_mf2wcs = r
        r_bf2wcs = Q.rotate(q_mf2wcs, r_bf2mf) + r_mf2wcs

        # Finally we can update the origin of the body model frame origin to reflect the desired position of the model
        #  frame origin in the world coordinate system
        body.r = r_bf2wcs
    else:
        body.r = r


def set_orientation(engine, body_name, q, use_model_frame=False):
    if body_name in engine.rigid_bodies:
        body = engine.rigid_bodies[body_name]
    else:
        raise RuntimeError('set_position() no such rigid body exist with that name')

    if use_model_frame:

        r_bf2wcs = body.r
        q_bf2wcs = body.q

        r_bf2mf = body.shape.r
        q_bf2mf = body.shape.q
        #                                                                      |q r| |x|
        # By definition we have the rigid body transformations T x = q*x + r = |0 q| |1|
        #
        #  T_bf2wcs =   T_mf2wcs T_bf2mf
        #
        #   |q_bf2wcs  r_bf2wcs|   |q_mf2wcs  r_mf2wcs|  |q_bf2mf   r_bf2mf |
        #   |0            1    | = |0            1    |  |0            1    |
        #
        # So
        #
        #   q_bf2wcs = q_mf2wcs q_bf2mf
        #   r_bf2wcs = q_mf2wcs r_bf2mf + r_mf2wcs
        #
        # From this we can solve for the current model to world coordinate transformation
        #
        q_mf2wcs = Q.prod(q_bf2wcs, Q.conjugate(q_bf2mf))
        r_mf2wcs = r_bf2wcs - Q.rotate(q_mf2wcs, r_bf2mf)
        #
        # Now we wish to change the orientation of the model frame origin  wrt the world frame
        #
        # so we must now have
        #
        #   q_mf2wcs = [qs, qx, qy, qz]
        #
        q_mf2wcs = Q.unit(q)
        #
        # Change the orientation of the model frame means that both the orientation and
        #  position of the body frame will change wrt. the world coordinate system.
        q_bf2wcs = Q.prod(q_mf2wcs, q_bf2mf)
        r_bf2wcs = Q.rotate(q_mf2wcs, r_bf2mf) + r_mf2wcs

        body.q = q_bf2wcs
        body.r = r_bf2wcs
    else:
        body.q = Q.unit(q)


def set_velocity(engine, body_name, v):
    if body_name in engine.rigid_bodies:
        body = engine.rigid_bodies[body_name]
    else:
        raise RuntimeError('set_velocity() no such rigid body exist with that name')
    body.v = v


def set_spin(engine, body_name, w):
    if body_name in engine.rigid_bodies:
        body = engine.rigid_bodies[body_name]
    else:
        raise RuntimeError('set_spin() no such rigid body exist with that name')
    body.w = w


def set_mass_properties_from_shape(engine, body_name, density):
    if body_name in engine.rigid_bodies:
        body = engine.rigid_bodies[body_name]
    else:
        raise RuntimeError('set_mass_properties_from_shape() no such rigid body exist with that name')
    if body.shape is None:
        raise RuntimeError('set_mass_properties_from_shape() rigid body did not have a shape')
    body.mass = body.shape.mass*density
    body.inertia = body.shape.inertia*density


def set_mass(engine, body_name, mass):
    if body_name in engine.rigid_bodies:
        body = engine.rigid_bodies[body_name]
    else:
        raise RuntimeError('set_mass() no such rigid body exist with that name')
    if mass <= 0.0:
        raise RuntimeError('set_mass() illegal mass value')
    body.mass = mass


def set_inertia(engine, body_name, inertia):
    if body_name in engine.rigid_bodies:
        body = engine.rigid_bodies[body_name]
    else:
        raise RuntimeError('set_inertia() no such rigid body exist with that name')

    if inertia[0] <= 0.0:
        raise RuntimeError('set_inertia() Illegal Ixx value')

    if inertia[1] <= 0.0:
        raise RuntimeError('set_inertia() Illegal Iyy value')

    if inertia[2] <= 0.0:
        raise RuntimeError('set_inertia() Illegal Izz value')

    body.inertia = inertia


def set_visual_material(engine, body_name, visual_material_name):
    if body_name in engine.rigid_bodies:
        body = engine.rigid_bodies[body_name]
    else:
        raise RuntimeError('set_visual_material() no such rigid body exist with that name')

    if visual_material_name is None:
        raise RuntimeError('set_visual_material() must give a visual material name')

    body.visual_material = visual_material_name


def create_gravity_force(engine, force_name, g, up):
    if force_name in engine.forces:
        raise RuntimeError('create_gravity(): Force already exist with that name')
    if g <= 0.0:
        raise RuntimeError('create_gravity(): Illegal value for gravitational acceleration')

    gravity = Gravity(force_name)
    gravity.up = V3.unit(up)
    gravity.g = g

    engine.forces[force_name] = gravity


def create_damping_force(engine, force_name, alpha, beta):
    if force_name in engine.forces:
        raise RuntimeError('create_damping(): Force already exist with that name')
    if alpha <= 0:
        raise RuntimeError('create_damping(): Illegal value for alpha')
    if beta <= 0:
        raise RuntimeError('create_damping(): Illegal value for beta')

    damping = Damping(force_name)
    damping.alpha = alpha
    damping.beta = beta

    engine.forces[force_name] = damping


def connect_force(engine, body_name, force_name):
    if body_name in engine.rigid_bodies:
        body = engine.rigid_bodies[body_name]
    else:
        raise RuntimeError('connect_force() no such rigid body exist with that name')

    if force_name in engine.forces:
        force = engine.forces[force_name]
    else:
        raise RuntimeError('connect_force() no such force exist with that name')

    if force in body.forces:
        raise RuntimeError('connect_force() force was already connected to body')

    body.forces.append(force)


def set_body_type(engine, body_name, body_type):
    if body_name not in engine.rigid_bodies:
        raise RuntimeError('set_body_type() no such rigid body exist with that name')

    body = engine.rigid_bodies[body_name]

    body.is_free = False
    body.is_fixed = False
    body.is_scripted = False

    if body_type in ['free', 'Free', 'FREE']:
        body.is_free = True
    elif body_type in ['fixed', 'Fixed', 'FIXED']:
        body.is_fixed = True
    elif body_type in ['scripted', 'Scripted', 'SCRIPTED']:
        body.is_scripted = True
    else:
        raise RuntimeError('set_body_type(): Unsupported body type found')


def is_fixed_body(engine, body_name):
    if body_name not in engine.rigid_bodies:
        raise RuntimeError('is_fixed_body() no such rigid body exist with that name')
    body = engine.rigid_bodies[body_name]
    return body.is_fixed


def is_free_body(engine, body_name):
    if body_name not in engine.rigid_bodies:
        raise RuntimeError('is_free_body() no such rigid body exist with that name')
    body = engine.rigid_bodies[body_name]
    return body.is_free


def is_scripted_body(engine, body_name):
    if body_name not in engine.rigid_bodies:
        raise RuntimeError('is_scripted_body() no such rigid body exist with that name')
    body = engine.rigid_bodies[body_name]
    return body.is_scripted


def set_body_active(engine, body_name, value):
    if body_name not in engine.rigid_bodies:
        raise RuntimeError('set_body_type() no such rigid body exist with that name')

    body = engine.rigid_bodies[body_name]
    body.is_active = value


def is_active_body(engine, body_name):
    if body_name not in engine.rigid_bodies:
        raise RuntimeError('is_active_body() no such rigid body exist with that name')
    body = engine.rigid_bodies[body_name]
    return body.is_active


def set_body_material(engine, body_name, material_name):
    if body_name not in engine.rigid_bodies:
        raise RuntimeError('set_body_material() no such rigid body exist with that name')

    body = engine.rigid_bodies[body_name]

    if not engine.material_library.exist_material(material_name):
        raise RuntimeError('set_body_material() no such material exist')

    body.material = material_name


def create_material_behavior(engine, A, B, epsilon, mu):
    if engine.material_library.exist_behaviour(A, B):
        raise RuntimeError('create_material_behavior() behaviour already exist')

    if epsilon < 0.0:
        raise RuntimeError('create_material_behavior() illegal epsilon value')

    if mu[0] < 0.0:
        raise RuntimeError('create_material_behavior() illegal mu_x value')

    if mu[1] < 0.0:
        raise RuntimeError('create_material_behavior() illegal mu_y value')

    if mu[2] < 0.0:
        raise RuntimeError('create_material_behavior() illegal mu_z value')

    tmp = [A, B]
    tmp.sort()
    key = tuple(tmp)

    behaviour = MaterialBehaviour()
    behaviour.epsilon = epsilon
    behaviour.mu = mu
    engine.material_library.storage[key] = behaviour


def set_finite_update(engine, body_name, value):
    if body_name not in engine.rigid_bodies:
        raise RuntimeError('set_finite_update_rotation_axis() no such rigid body exist with that name')
    body = engine.rigid_bodies[body_name]
    body.use_finite_update = value


def set_finite_update_rotation_axis(engine, body_name, axis):
    if body_name not in engine.rigid_bodies:
        raise RuntimeError('set_finite_update_rotation_axis() no such rigid body exist with that name')
    body = engine.rigid_bodies[body_name]
    body.finite_update_rotation_axis = V3.unit(axis)


def create_keyframe(engine, body_name, time, r, q):
    if body_name not in engine.rigid_bodies:
        raise RuntimeError('create_keyframe() no such rigid body exist with that name')

    body = engine.rigid_bodies[body_name]

    if not body.is_scripted:
        raise RuntimeError('create_keyframe() body was not scripted')

    if body.scripted_motion is None:
        body.scripted_motion = KeyframeMotion()

    body.scripted_motion.create_keyframe(time, r, Q.unit(q))

    body.scripted_motion.keyframes.sort()

    body.r = body.scripted_motion.keyframes[0][1]
    body.q = body.scripted_motion.keyframes[0][2]


def generate_unique_name(engine, name):
    import datetime
    import random
    n = random.random()
    unique_name = name + '_' + str(n) + '_' + str(datetime.datetime.now())
    return unique_name


def compute_shape_extends(engine, shape_name, use_model_frame=True):
    if shape_name in engine.shapes:
        shape = engine.shapes[shape_name]
    else:
        raise RuntimeError('get_shape_extends() no such shape exist with that name')
    if use_model_frame:
        model_mesh = MESH.join([shape.mesh])
        MESH.rotate(model_mesh, shape.q)
        MESH.translate(model_mesh, shape.r)
        (l, u) = MESH.aabb(model_mesh)
    else:
        (l, u) = MESH.aabb(shape.mesh)
    return u-l


def rotate_shape(engine, shape_name, q):
    if shape_name in engine.shapes:
        shape = engine.shapes[shape_name]
    else:
        raise RuntimeError('resize_shape() no such shape exist with that name')
    # Make sure mesh are in model frame
    MESH.rotate(shape.mesh, shape.q)
    MESH.translate(shape.mesh, shape.r)
    # Rotate mesh in its model frame
    MESH.rotate(shape.mesh, q)
    # Compute new body frame mesh
    transform_shape_into_body_frame(shape)


def scale_shape(engine, shape_name, sx, sy, sz):
    if shape_name in engine.shapes:
        shape = engine.shapes[shape_name]
    else:
        raise RuntimeError('resize_shape() no such shape exist with that name')
    # Make sure mesh are in model frame
    MESH.rotate(shape.mesh, shape.q)
    MESH.translate(shape.mesh, shape.r)
    # Compute the scaling in model frame
    MESH.scale(shape.mesh, sx, sy, sz)
    # Compute new body frame mesh
    transform_shape_into_body_frame(shape)


def center_shape(engine, shape_name):
    if shape_name in engine.shapes:
        shape = engine.shapes[shape_name]
    else:
        raise RuntimeError('center_shape() no such shape exist with that name')
    # Make sure mesh are in model frame
    MESH.rotate(shape.mesh, shape.q)
    MESH.translate(shape.mesh, shape.r)
    # Compute the new center in model frame
    (l, u) = MESH.aabb(shape.mesh)
    center = (l+u)/2.0
    MESH.translate(shape.mesh, -center)
    # Compute new body frame mesh
    transform_shape_into_body_frame(shape)


def resize_shape_to_unit_box(engine, shape_name, keep_model_center=False):
    if shape_name in engine.shapes:
        shape = engine.shapes[shape_name]
    else:
        raise RuntimeError('resize_shape() no such shape exist with that name')
    # Make sure mesh are in model frame
    MESH.rotate(shape.mesh, shape.q)
    MESH.translate(shape.mesh, shape.r)
    # Compute the scaling in model frame
    (l, u) = MESH.aabb(shape.mesh)
    center = (l+u)/2.0
    MESH.translate(shape.mesh, -center)
    s = 1.0 / np.max(u-l)
    MESH.scale(shape.mesh, s, s, s)
    if keep_model_center:
        MESH.translate(shape.mesh, center)
    # Compute new body frame mesh
    transform_shape_into_body_frame(shape)


def create_ball_joint(engine, joint_name, socket_body_name, plug_body_name):
    if joint_name in engine.joints:
        raise RuntimeError('create_ball_joint() joint already exist with that name')

    if socket_body_name not in engine.rigid_bodies:
        raise RuntimeError('create_ball_joint() no such rigid body exist with name ' + socket_body_name)

    socket_body = engine.rigid_bodies[socket_body_name]

    if plug_body_name not in engine.rigid_bodies:
        raise RuntimeError('create_ball_joint() no such rigid body exist with name ' + plug_body_name)

    plug_body = engine.rigid_bodies[plug_body_name]

    engine.joints[joint_name] = BallJoint(joint_name, socket_body, plug_body)


def set_socket_connector(engine, joint_name, r, q):
    if joint_name not in engine.joints:
        raise RuntimeError('set_socket_connector() No joint with that name')

    joint = engine.joints[joint_name]
    # TODO Need to carefully design control over rigging frame
    joint.socket.transform.r = r
    joint.socket.transform.q = q


def set_plug_connector(engine, joint_name, r, q):
    if joint_name not in engine.joints:
        raise RuntimeError('set_plug_connector() No joint with that name')

    joint = engine.joints[joint_name]
    # TODO Need to carefully design control over rigging frame
    joint.plug.transform.r = r
    joint.plug.transform.q = q


def set_joint_error_reduction(engine, joint_name, error_reduction):
    if joint_name not in engine.joints:
        raise RuntimeError('set_joint_parameters() No joint with that name')

    joint = engine.joints[joint_name]

    if error_reduction < 0.0:
        raise RuntimeError('set_joint_parameters() Illegal error_reduction value')

    if error_reduction > 1.0:
        raise RuntimeError('set_joint_parameters() Illegal error_reduction value')

    joint.error_reduction = error_reduction





