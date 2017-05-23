import sys

sys.path.append('/usr/local/lib/python/')


import openmesh as OM
import numpy as np
import scipy.sparse as sparse
import REESMath.vector3 as V3
import REESMath.quaternion as Q
import REESMath.matrix3 as M3
import REESMass.mass as MASS
import REESMesh.mesh as MESH
import REESMesh.factory as MESH_FACTORY
import REESUtility.util as UTIL
from REESSimulation.types import *
from REESMath.functions import sinc
from math import cos
import xml.etree.ElementTree as ET
import os.path


def transform_shape_into_body_frame(shape):
    prop = MASS.compute_mass_properties(shape.mesh, 1.0)
    (shape.r, shape.q, shape.mass, shape.inertia) = MASS.xform_model_2_body_space(prop)
    #
    # shape.r and shape.q gives the rigid body transform from body to model space
    # We need to do the inverse transform here
    #
    MESH.translate(shape.mesh, -shape.r)
    MESH.rotate(shape.mesh, Q.conjugate(shape.q))


def get_position_vector(engine):
    x = np.zeros(len(engine.rigid_bodies)*7, dtype=np.float64)
    k = 0
    for body in engine.rigid_bodies.values():
        offset = 7*k
        x[offset:offset+3] = body.r
        x[offset+3:offset+7] = body.q
        k += 1
    return x


def set_position_vector(engine, x):
    k = 0
    for body in engine.rigid_bodies.values():
        offset = 7*k
        body.r = x[offset:offset+3]
        body.q = x[offset+3:offset+7]
        k += 1


def get_velocity_vector(engine):
    u = np.zeros((len(engine.rigid_bodies)*6,), dtype=np.float64)
    k = 0
    for body in engine.rigid_bodies.values():
        offset = 6*k
        u[offset:offset+3] = body.v
        u[offset+3:offset+6] = body.w
        k += 1
    return u


def set_velocity_vector(engine, u):
    k = 0
    for body in engine.rigid_bodies.values():
        offset = 6*k
        body.v = u[offset:offset+3]
        body.w = u[offset+3:offset+6]
        k += 1


def compute_keyframe_motion(keyframes, time):
    if time < keyframes[0][0]:
        r = keyframes[0][1]
        q = keyframes[0][2]
        return r, q, V3.zero(), V3.zero()

    N = len(keyframes)

    if time >= keyframes[N-1][0]:
        r = keyframes[N-1][1]
        q = keyframes[N-1][2]
        return r, q, V3.zero(), V3.zero()

    for n in range(N):
        keyframe0 = keyframes[n]
        keyframe1 = keyframes[(n+1) % N]

        time0 = keyframe0[0]
        r0 = keyframe0[1]
        q0 = keyframe0[2]

        time1 = keyframe1[0]
        r1 = keyframe1[1]
        q1 = keyframe1[2]

        if time0 <= time < time1:
            dT = time1 - time0
            t = (time - time0) / dT
            r = r0 + (r1-r0)*t
            q = Q.slerp(q0, q1, t)
            v = (r1-r0)/dT
            dQ = Q.prod(q1, Q.conjugate(q0))
            (theta, n) = Q.to_angle_axis(dQ)
            w = n*(theta/dT)
            return r, q, v, w
    raise RuntimeError('compute_keyframe_motion(): Internal error, could not find keyframes to interpolate from')


def record_motion(engine):

    if not engine.motion_recorder.on:
        return

    for body in engine.rigid_bodies.values():
        engine.motion_recorder.record(engine.solver_params.current_time, body)


def save_restart_data(engine):

    if not engine.solver_params.save_restart:
        return

    filename = os.path.join(engine.solver_params.restart_path, engine.solver_params.restart_filename)

    ext = os.path.splitext(filename)[-1].lower()

    if ext != '.xml':
        raise RuntimeError('save_restart_data(): save: file was not a xml file')

    root = ET.Element('restart')

    root.attrib['time'] = str(engine.solver_params.current_time)

    for body in engine.rigid_bodies.values():
        body_tag = ET.SubElement(root, 'body')
        body_tag.attrib['name'] = body.name
        body_tag.attrib['r'] = UTIL.array2string(body.r)
        body_tag.attrib['q'] = UTIL.array2string(body.q)
        body_tag.attrib['v'] = UTIL.array2string(body.v)
        body_tag.attrib['w'] = UTIL.array2string(body.w)

    UTIL.xml_pretty_indent(root)
    tree = ET.ElementTree(root)
    tree.write(filename)


def position_update(engine, x, u, dt):

    new_time = engine.solver_params.current_time + dt

    k = 0

    for body in engine.rigid_bodies.values():
        x_offset = 7*k
        u_offset = 6*k

        r = x[x_offset:x_offset+3]
        q = x[x_offset+3:x_offset+7]
        v = u[u_offset:u_offset+3]
        w = u[u_offset+3:u_offset+6]

        if body.is_free and not body.use_finite_update:

            r += v * dt

            q += Q.prod(
                Q.from_vector3(w),
                q
            ) * dt * 0.5

        elif body.is_free and body.use_finite_update and body.finite_update_rotation_axis is not None:

            r += v * dt

            # Split angular velocity into parallel (finite update) and orthogonal components (infinitesimal update)
            w_finite = np.dot(w, body.finite_update_rotation_axis) * body.finite_update_rotation_axis
            w_infinitesimal = w - w_finite

            # Apply finite rotation around rotation axis
            theta = np.linalg.norm(w_finite)*dt*0.5
            vec = w_finite * sinc(theta) * dt * 0.5
            q_finite = Q.from_array([cos(theta), vec[0], vec[1], vec[2]])
            q = Q.prod(q_finite, q)

            # Apply infinitesimal rotation update
            q += Q.prod(
                Q.from_vector3(w_infinitesimal),
                q
            ) * dt * 0.5

        elif body.is_free and body.use_finite_update and body.finite_update_rotation_axis is None:

            r += v * dt

            # Apply finite rotation around axis n = w / |w| with angle theta = |w|*dt
            theta = np.linalg.norm(w)*dt*0.5
            vec = body.w * sinc(theta) * dt * 0.5
            q_finite = Q.from_array([cos(theta), vec[0], vec[1], vec[2]])
            q = Q.prod(q_finite, q)

        elif body.is_scripted:

            if body.scripted_motion is None:
                raise RuntimeError('position_update(): scripted body did not have any scripted motion?')
            else:
                r, q, v, w = compute_keyframe_motion(body.scripted_motion.keyframes, new_time)

            u[u_offset:u_offset+3] = v
            u[u_offset+3:u_offset+6] = w

        x[x_offset:x_offset+3] = r
        x[x_offset+3:x_offset+7] = Q.unit(q)

        k += 1


def compute_total_external_forces(engine, x, u):
    f_ext = np.zeros((len(engine.rigid_bodies)*6,), dtype=np.float64)
    k = 0
    for body in engine.rigid_bodies.values():
        x_offset = 7*k
        u_offset = 6*k

        F = V3.zero()
        T = V3.zero()

        if body.is_free:
            r = x[x_offset:x_offset+3]
            q = x[x_offset+3:x_offset+7]
            v = u[u_offset:u_offset+3]
            w = u[u_offset+3:u_offset+6]

            for force_type in body.forces:
                (Fi, Ti) = force_type.compute(body, r, q, v, w)
                F += Fi
                T += Ti
            R = Q.to_matrix(q)
            I = MASS.update_inertia_tensor(R, body.inertia)
            T -= np.cross(w, np.dot(I, w), axis=0)

        f_ext[u_offset:u_offset+3] = F
        f_ext[u_offset+3:u_offset+6] = T

        k += 1
    return f_ext


def compute_inverse_mass_matrix(engine, x):
    N = len(engine.rigid_bodies)
    blocks = [np.zeros((3, 3), dtype=np.float64) for _ in range(N*2)]
    k = 0
    for body in engine.rigid_bodies.values():
        x_offset = 7*k  # Position offset into x-array

        if body.is_free:
            q = x[x_offset+3:x_offset+7]  # Extract rotation part
            R = Q.to_matrix(q)
            I = MASS.update_inertia_tensor(R, 1.0 / body.inertia)
            m = 1.0 / body.mass
            blocks[2*k] = M3.diag(m, m, m)
            blocks[2*k+1] = I
        k += 1
    D = sparse.block_diag(blocks)  # TODO Kenny 2017-02-12: Verify this is the most efficient way of doing this!
    W = D.tobsr(blocksize=(3, 3))
    return W


def simulation_stepper(engine, dt):
    x = get_position_vector(engine)
    u = get_velocity_vector(engine)

    f_ext = compute_total_external_forces(engine, x, u)
    W = compute_inverse_mass_matrix(engine, x)

    # TODO Do collision detection and solve for constraint forces here

    u += W.dot(f_ext*dt)

    position_update(engine, x, u, dt)

    set_position_vector(engine, x)
    set_velocity_vector(engine, u)

    engine.solver_params.current_time += dt


def motion_player(engine, dt):

    new_time = engine.solver_params.current_time + dt

    if len(engine.motion_recorder.storage) <= 0:
        raise RuntimeError('player() No recorded motion data exist')

    for body in engine.rigid_bodies.values():

        if body.name not in engine.motion_recorder.storage:
            raise RuntimeError('player() Body was not in recorded motion storage')

        keyframe_motion = engine.motion_recorder.storage[body.name]

        if keyframe_motion is None:
            raise RuntimeError('player() could not find motion data for this body')

        r, q, v, w = compute_keyframe_motion(keyframe_motion.keyframes, new_time)

        body.r = r
        body.q = Q.unit(q)
        body.v = v
        body.w = w


def simulate(engine):
    if engine is None:
        return False

    if not engine.solver_params.on:
        return False

    if not engine.solver_params.mode == 'simulate':
        return False

    if engine.solver_params.current_time >= engine.solver_params.total_time:
        return False

    if engine.solver_params.current_time == 0.0:
        record_motion(engine)

    T = engine.solver_params.total_time - engine.solver_params.current_time
    T_safe = min(T, 1.0 / engine.solver_params.fps)
    dt = engine.solver_params.time_step

    while T_safe > 0:
        dt_safe = min(dt, T_safe)
        simulation_stepper(engine, dt_safe)
        T_safe = 0 if dt_safe < dt else T_safe - dt_safe

    record_motion(engine)

    save_restart_data(engine)

    # Just for debugging visualization
    if len(engine.contact_points) == 0:
        engine.contact_points = [ContactPoint(V3.rand(-1.0, 1.0), V3.rand(-1.0, 1.0), ) for _ in range(100)]

    return True


def play(engine):
    if engine is None:
        return False

    if not engine.solver_params.on:
        return False

    if not engine.solver_params.mode == 'play':
        return False

    if engine.solver_params.current_time >= engine.solver_params.total_time:
        return False

    T = engine.solver_params.total_time - engine.solver_params.current_time
    T_safe = min(T, 1.0 / engine.solver_params.fps)
    dt = engine.solver_params.time_step

    while T_safe > 0:
        dt_safe = min(dt, T_safe)
        motion_player(engine, dt_safe)
        T_safe = 0 if dt_safe < dt else T_safe - dt_safe

    return True


