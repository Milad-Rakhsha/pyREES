import sys

sys.path.append('/usr/local/lib/python/')


import openmesh as OM
import xml.etree.ElementTree as ET
import textwrap
import os.path
import REESMass.mass as MASS
import REESMesh.mesh as MESH
import REESUtility.util as UTIL
from REESSimulation.types import *
from REESSimulation.api import *
from REESSimulation.procedural import *
import REESMath.vector3 as V3
import REESMath.quaternion as Q
import REESMath.coordsys as C
import numpy as np
from math import pi, cos, sin, sqrt


def __load_solver_parameters(engine, engine_tag):
    params = SolverParameters()
    engine.solver_params = params
    solver_param_tag = engine_tag.find('solver_parameters')
    if solver_param_tag is None:
        return
    params.set_state(UTIL.bool_from_xml(solver_param_tag, 'on', params.on))
    params.fps = UTIL.float_from_xml(solver_param_tag, 'fps', params.fps)
    params.time_step = UTIL.float_from_xml(solver_param_tag, 'time_step', params.time_step)
    params.current_time = UTIL.float_from_xml(solver_param_tag, 'current_time', params.current_time)
    params.total_time = UTIL.float_from_xml(solver_param_tag, 'total_time', params.total_time)
    params.mode = UTIL.string_from_xml(solver_param_tag, 'mode', params.mode)
    params.load_restart = UTIL.bool_from_xml(solver_param_tag, 'load_restart', params.load_restart)
    params.save_restart = UTIL.bool_from_xml(solver_param_tag, 'save_restart', params.save_restart)
    params.restart_path = UTIL.string_from_xml(solver_param_tag, 'restart_path', params.restart_path)
    params.restart_filename = UTIL.string_from_xml(solver_param_tag, 'restart_filename', params.restart_filename)


def __save_solver_parameters(engine, engine_tag):
    solver_param_tag = ET.SubElement(engine_tag, 'solver_parameters')
    solver_param_tag.attrib['on'] = str(engine.solver_params.on)
    solver_param_tag.attrib['fps'] = str(engine.solver_params.fps)
    solver_param_tag.attrib['time_step'] = str(engine.solver_params.time_step)
    solver_param_tag.attrib['current_time'] = str(engine.solver_params.current_time)
    solver_param_tag.attrib['total_time'] = str(engine.solver_params.total_time)
    solver_param_tag.attrib['mode'] = str(engine.solver_params.mode)
    solver_param_tag.attrib['load_restart'] = str(engine.solver_params.load_restart)
    solver_param_tag.attrib['save_restart'] = str(engine.solver_params.save_restart)
    solver_param_tag.attrib['restart_path'] = engine.solver_params.restart_path
    solver_param_tag.attrib['restart_filename'] = engine.solver_params.restart_filename


def __load_motion_recorder(engine, engine_tag):
    motion_recorder = MotionRecorder()
    engine.motion_recorder = motion_recorder
    motion_recorder_tag = engine_tag.find('motion_recorder')
    if motion_recorder_tag is None:
        return
    motion_recorder.path = UTIL.string_from_xml(motion_recorder_tag, 'path', motion_recorder.path)
    motion_recorder.filename = UTIL.string_from_xml(motion_recorder_tag, 'filename', motion_recorder.filename)
    motion_recorder.set_state(UTIL.bool_from_xml(motion_recorder_tag, 'on', motion_recorder.on))
    motion_recorder.load = UTIL.bool_from_xml(motion_recorder_tag, 'load', motion_recorder.load)
    motion_recorder.save = UTIL.bool_from_xml(motion_recorder_tag, 'save', motion_recorder.save)

    if motion_recorder.load:

        if not os.path.exists(engine.motion_recorder.path):
            raise RuntimeError('motion path did not exist')

        filename = os.path.join(engine.motion_recorder.path, engine.motion_recorder.filename)

        ext = os.path.splitext(filename)[-1].lower()

        if ext != '.xml':
            raise RuntimeError('__load_motion_recorder: motion file=' + filename + ' was not a xml file')
        if not os.path.isfile(filename):
            raise RuntimeError('__load_motion_recorder: motion filename=' + filename + ' was not a file')
        if not os.path.exists(filename):
            raise RuntimeError('__load_motion_recorder: motion file' + filename + ' did not exist')

        motion_xml = ET.parse(filename)
        motion_root = motion_xml.getroot()

        for channel_tag in motion_root.iter('channel'):

            body_name = UTIL.string_from_xml(channel_tag, 'body', None)
            if body_name is None:
                raise RuntimeError('__load_motion_recorder(): body is missing on channel element')

            if body_name in engine.motion_recorder.storage:
                raise RuntimeError('__load_motion_recorder(): body motion data have already been loaded')

            recorded_motion = KeyframeMotion()
            for keyframe_tag in channel_tag.iter('keyframe'):
                time = UTIL.float_from_xml(keyframe_tag, 'time', None)
                r = UTIL.vector3_from_xml(keyframe_tag, 'r', None)
                q = UTIL.quaternion_from_xml(keyframe_tag, 'q', None)
                if time is None:
                    raise RuntimeError('__load_motion_recorder(): Keyframe missing time attribute')
                if r is None:
                    raise RuntimeError('__load_motion_recorder(): Keyframe missing r attribute')
                if q is None:
                    raise RuntimeError('__load_motion_recorder(): Keyframe missing q attribute')
                recorded_motion.create_keyframe(time, r, q)

            recorded_motion.keyframes.sort()

            motion_recorder.storage[body_name] = recorded_motion


def __save_motion_recorder(engine, engine_tag):
    motion_recorder_tag = ET.SubElement(engine_tag, 'motion_recorder')
    motion_recorder_tag.attrib['path'] = engine.motion_recorder.path
    motion_recorder_tag.attrib['filename'] = engine.motion_recorder.filename
    motion_recorder_tag.attrib['on'] = str(engine.motion_recorder.on)
    motion_recorder_tag.attrib['load'] = str(engine.motion_recorder.load)
    motion_recorder_tag.attrib['save'] = str(engine.motion_recorder.save)
    if engine.motion_recorder.save:

        filename = os.path.join(engine.motion_recorder.path, engine.motion_recorder.filename)

        ext = os.path.splitext(filename)[-1].lower()

        if ext != '.xml':
            raise RuntimeError('__save_motion_recorder: motion file was not a xml file')

        motion_root = ET.Element('motion')

        for body_name, motion in engine.motion_recorder.storage.items():

            channel_tag = ET.SubElement(motion_root, 'channel')
            channel_tag.attrib['body'] = body_name

            channel = engine.motion_recorder.storage[body_name]

            for keyframe in channel.keyframes:
                key_tag = ET.SubElement(channel_tag, 'keyframe')
                key_tag.attrib['time'] = str(keyframe[0])
                key_tag.attrib['r'] = UTIL.array2string(keyframe[1])
                key_tag.attrib['q'] = UTIL.array2string(keyframe[2])

        UTIL.xml_pretty_indent(motion_root)
        tree = ET.ElementTree(motion_root)
        tree.write(filename)


def __load_profiler(engine, engine_tag):
    engine.profiler = Profiler()
    profiler_tag = engine_tag.find('profiler')
    if profiler_tag is None:
        return
    engine.profiler.path = UTIL.string_from_xml(profiler_tag, 'path', engine.profiler.path)
    engine.profiler.filename = UTIL.string_from_xml(profiler_tag, 'filename', engine.profiler.filename)
    engine.profiler.set_state(UTIL.bool_from_xml(profiler_tag, 'on', engine.profiler.on))


def __save_profiler(engine, engine_tag):
    profiler_tag = ET.SubElement(engine_tag, 'profiler')
    profiler_tag.attrib['path'] = engine.profiler.path
    profiler_tag.attrib['filename'] = engine.profiler.filename
    profiler_tag.attrib['on'] = str(engine.profiler.on)


def __load_material_library(engine, engine_tag):
    engine.material_library = MaterialLibrary()
    material_library_tag = engine_tag.find('material_library')
    if material_library_tag is None:
        return
    for behavior_tag in material_library_tag.iter('behaviour'):
        behavior = MaterialBehaviour()
        key_string_list = UTIL.string_from_xml(behavior_tag, 'materials', None)
        if key_string_list is None:
            raise RuntimeError('__load_material_library(): Missing materials attribute')
        key = UTIL.string_list_2_sorted_tuple(key_string_list)
        behavior.mu = UTIL.vector3_from_xml(behavior_tag, 'friction', behavior.mu)
        behavior.epsilon = UTIL.float_from_xml(behavior_tag, 'restitution', behavior.epsilon)
        engine.material_library.storage[key] = behavior


def __save_material_library(engine, engine_tag):
    library_tag = ET.SubElement(engine_tag, 'material_library')
    if engine.material_library is not None:
        for key, behavior in engine.material_library.storage.items():
            behavior_tag = ET.SubElement(library_tag, 'behaviour')
            behavior_tag.attrib['materials'] = str(key)
            behavior_tag.attrib['friction'] = UTIL.array2string(behavior.mu)
            behavior_tag.attrib['restitution'] = str(behavior.epsilon)


def __load_gravity_force(engine, gravity_tag):
    name = UTIL.string_from_xml(gravity_tag, 'name', None)

    if name is None:
        raise RuntimeError('__load_gravity_force() missing name attribute')

    if name in engine.forces:
        raise RuntimeError('__load_gravity_force(): force with name =' + name + ' already exist')

    gravity = Gravity(name)
    gravity.up = UTIL.vector3_from_xml(gravity_tag, 'up', gravity.up)
    gravity.g = UTIL.float_from_xml(gravity_tag, 'g', gravity.g)

    engine.forces[name] = gravity


def __save_gravity_force(gravity, engine_tag):
    gravity_tag = ET.SubElement(engine_tag, 'gravity')
    gravity_tag.attrib['name'] = gravity.name
    gravity_tag.attrib['g'] = str(gravity.g)
    gravity_tag.attrib['up'] = UTIL.array2string(gravity.up)


def __load_damping_force(engine, damping_tag):
    name = UTIL.string_from_xml(damping_tag, 'name', None)

    if name is None:
        raise RuntimeError('__load_damping_force() missing name')

    if name in engine.forces:
        raise RuntimeError('__load_damping_force(): force with name =' + name + ' already exist')

    damping = Damping(name)
    damping.alpha = UTIL.float_from_xml(damping_tag, 'alpha', damping.alpha)
    damping.beta = UTIL.float_from_xml(damping_tag, 'beta', damping.beta)

    engine.forces[name] = damping


def __save_damping_force(damping, engine_tag):
    damping_tag = ET.SubElement(engine_tag, 'damping')
    damping_tag.attrib['name'] = damping.name
    damping_tag.attrib['alpha'] = str(damping.alpha)
    damping_tag.attrib['beta'] = str(damping.beta)


def __load_shape(engine, shape_tag):
    name = UTIL.string_from_xml(shape_tag, 'name', None)

    if name is None:
        raise RuntimeError('__load_shape() missing name')

    if name in engine.shapes:
        raise RuntimeError('__load_shape(): shape with name=' + name + ' already exist')

    shape = Shape(name)

    mesh_file_tag = shape_tag.find('mesh_file')
    mesh_tag = shape_tag.find('mesh')
    cylinder_tag = shape_tag.find('cylinder')
    capsule_tag = shape_tag.find('capsule')
    cone_tag = shape_tag.find('cone')
    conical_tag = shape_tag.find('conical')
    sphere_tag = shape_tag.find('sphere')
    ellipsoid_tag = shape_tag.find('ellipsoid')
    tetrahedron_tag = shape_tag.find('tetrahedron')
    box_tag = shape_tag.find('box')
    cuboid_tag = shape_tag.find('cuboid')
    convex_tag = shape_tag.find('convex')

    if mesh_file_tag is not None:
        filename = UTIL.string_from_xml(mesh_file_tag, 'filename', None)
        if filename is None:
            raise RuntimeError('__load_shape() missing filename')
        shape.mesh = MESH.read_obj(filename)
    elif mesh_tag is not None:
        lut = dict()

        shape.mesh = OM.TriMesh()
        for vertex_tag in mesh_tag.iter('vertex'):
            idx = UTIL.int_from_xml(vertex_tag, 'idx', None)

            if idx is None:
                raise RuntimeError('__load_shape() missing idx attribute on vertex tag')

            position = UTIL.vector3_from_xml(vertex_tag, 'r', None)

            if position is None:
                raise RuntimeError('__load_shape() missing position attribute on vertex tag')

            vh = shape.mesh.add_vertex(OM.TriMesh.Point(position[0], position[1], position[2]))
            lut[idx] = vh

        for triangle_tag in mesh_tag.iter('triangle'):
            i = UTIL.int_from_xml(triangle_tag, 'i', None)
            j = UTIL.int_from_xml(triangle_tag, 'j', None)
            k = UTIL.int_from_xml(triangle_tag, 'k', None)

            if i is None:
                raise RuntimeError('__load_shape() missing i attribute on triangle tag')
            if j is None:
                raise RuntimeError('__load_shape() missing j attribute on triangle tag')
            if k is None:
                raise RuntimeError('__load_shape() missing k attribute on triangle tag')

            shape.mesh.add_face(lut[i], lut[j], lut[k])

        shape.mesh.request_face_normals()
        shape.mesh.update_face_normals()
    elif cylinder_tag is not None:
        radius = UTIL.float_from_xml(cylinder_tag,'radius', 1.0)
        height = UTIL.float_from_xml(cylinder_tag,'height', 2.0)
        slices = UTIL.int_from_xml(cylinder_tag,'slices', 12)
        shape.mesh = MESH_FACTORY.make_cylinder(radius, height, slices)
    elif capsule_tag is not None:
        radius = UTIL.float_from_xml(capsule_tag, 'radius', 1.0)
        height = UTIL.float_from_xml(capsule_tag, 'height', 2.0)
        slices = UTIL.int_from_xml(capsule_tag, 'slices', 12)
        segments = UTIL.int_from_xml(capsule_tag, 'segments', 12)

        shape.mesh = MESH_FACTORY.make_capsule(radius, height, slices, segments)
    elif cone_tag is not None:
        radius = UTIL.float_from_xml(cone_tag, 'radius', 1.0)
        height = UTIL.float_from_xml(cone_tag, 'height', 2.0)
        slices = UTIL.int_from_xml(cone_tag, 'slices', 12)

        shape.mesh = MESH_FACTORY.make_cone(radius, height, slices)
    elif conical_tag is not None:
        bottom_radius = UTIL.float_from_xml(conical_tag, 'bottom_radius', 1.0)
        top_radius = UTIL.float_from_xml(conical_tag, 'top_radius', 2.0)
        height = UTIL.float_from_xml(conical_tag, 'height', 2.0)
        slices = UTIL.int_from_xml(conical_tag, 'slices', 12)

        shape.mesh = MESH_FACTORY.make_conical(bottom_radius, top_radius, height, slices)
    elif sphere_tag is not None:
        radius = UTIL.float_from_xml(sphere_tag, 'bottom_radius', 1.0)
        slices = UTIL.int_from_xml(sphere_tag, 'slices', 12)
        segments = UTIL.int_from_xml(sphere_tag, 'segments', 12)

        shape.mesh = MESH_FACTORY.make_sphere(radius, slices, segments)
    elif ellipsoid_tag is not None:
        a = UTIL.float_from_xml(ellipsoid_tag, 'a', 1.0)
        b = UTIL.float_from_xml(ellipsoid_tag, 'b', 1.0)
        c = UTIL.float_from_xml(ellipsoid_tag, 'c', 1.0)
        slices = UTIL.int_from_xml(ellipsoid_tag, 'slices', 12)
        segments = UTIL.int_from_xml(ellipsoid_tag, 'segments', 12)

        shape.mesh = MESH_FACTORY.make_ellipsoid(a, b, c, slices, segments)
    elif tetrahedron_tag is not None:
        p0 = UTIL.vector3_from_xml(tetrahedron_tag, 'p0', V3.make(0.0, 0.0, 0.0))
        p1 = UTIL.vector3_from_xml(tetrahedron_tag, 'p1', V3.make(0.0, 0.0, 1.0))
        p2 = UTIL.vector3_from_xml(tetrahedron_tag, 'p2', V3.make(1.0, 0.0, 0.0))
        p3 = UTIL.vector3_from_xml(tetrahedron_tag, 'p3', V3.make(0.0, 1.0, 0.0))

        shape.mesh = MESH_FACTORY.make_tetrahedron(p0, p1, p2, p3)
    elif box_tag is not None:

        width = UTIL.float_from_xml(box_tag, 'width', 1.0)
        height = UTIL.float_from_xml(box_tag, 'height', 1.0)
        depth = UTIL.float_from_xml(box_tag, 'depth', 1.0)

        shape.mesh = MESH_FACTORY.make_box(width, height, depth)
    elif cuboid_tag is not None:
        p0 = UTIL.vector3_from_xml(cuboid_tag, 'p0', None)
        p1 = UTIL.vector3_from_xml(cuboid_tag, 'p1', None)
        p2 = UTIL.vector3_from_xml(cuboid_tag, 'p2', None)
        p3 = UTIL.vector3_from_xml(cuboid_tag, 'p3', None)
        p4 = UTIL.vector3_from_xml(cuboid_tag, 'p4', None)
        p5 = UTIL.vector3_from_xml(cuboid_tag, 'p5', None)
        p6 = UTIL.vector3_from_xml(cuboid_tag, 'p6', None)
        p7 = UTIL.vector3_from_xml(cuboid_tag, 'p7', None)
        if p0 is None:
            raise RuntimeError('__load_shape(): missing p0 on cuboid element')
        if p1 is None:
            raise RuntimeError('__load_shape(): missing p1 on cuboid element')
        if p2 is None:
            raise RuntimeError('__load_shape(): missing p2 on cuboid element')
        if p3 is None:
            raise RuntimeError('__load_shape(): missing p3 on cuboid element')
        if p4 is None:
            raise RuntimeError('__load_shape(): missing p4 on cuboid element')
        if p5 is None:
            raise RuntimeError('__load_shape(): missing p5 on cuboid element')
        if p6 is None:
            raise RuntimeError('__load_shape(): missing p6 on cuboid element')
        if p7 is None:
            raise RuntimeError('__load_shape(): missing p7 on cuboid element')
        shape.mesh = MESH_FACTORY.make_cuboid(p0, p1, p2, p3, p4, p5, p6, p7)
    elif convex_tag is not None:
        points = []
        for point_tag in convex_tag.iter('point'):
            p = UTIL.vector3_from_xml(point_tag, 'p', None)
            if p is None:
                raise RuntimeError('__load_shape(): missing p on point element')
            points.append(p)
        shape.mesh = MESH_FACTORY.make_convex_hull(points)
    else:
        raise RuntimeError('__load_shape(): Unrecognized shape data')

    transform_shape_into_body_frame(shape)
    engine.shapes[name] = shape


def __save_shape(shape, engine_tag):
    shape_tag = ET.SubElement(engine_tag, 'shape')
    shape_tag.attrib['name'] = shape.name

    mesh_tag = ET.SubElement(shape_tag, 'mesh')
    vertices_tag = ET.SubElement(mesh_tag, 'vertices')
    triangles_tag = ET.SubElement(mesh_tag, 'triangles')

    for vh in shape.mesh.vertices():
        r = MESH.get_vertex_coords(shape.mesh, vh)
        vertex_tag = ET.SubElement(vertices_tag, 'vertex')
        vertex_tag.attrib['idx'] = str(vh.idx())
        vertex_tag.attrib['r'] = UTIL.array2string(r)

    for fh in shape.mesh.faces():
        indices = []
        for vh in shape.mesh.fv(fh):
            indices.append(vh.idx())
        if len(indices) != 3:
            raise RuntimeError('save(): Non triangle face encountered')
        triangle_tag = ET.SubElement(triangles_tag, 'triangle')
        triangle_tag.attrib['i'] = str(indices[0])
        triangle_tag.attrib['j'] = str(indices[1])
        triangle_tag.attrib['k'] = str(indices[2])


def __load_body(engine, body_tag):
    name = UTIL.string_from_xml(body_tag, 'name', None)

    if name is None:
        raise RuntimeError('__load_body(): Missing name on body')

    if name in engine.rigid_bodies:
        raise RuntimeError('__load_body(): body with name=' + name + ' already exist')

    body = RigidBody(name)

    body_type = UTIL.string_from_xml(body_tag, 'type', 'free')

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
        raise RuntimeError('__load_body(): Unsupported body type found')

    body.is_active = UTIL.bool_from_xml(body_tag, 'active', body.is_active)
    body.use_finite_update = UTIL.bool_from_xml(body_tag, 'use_finite_update', body.use_finite_update)
    if body.finite_update_rotation_axis is not None:
        body.finite_update_rotation_axis = UTIL.vector3_from_xml(body_tag, 'finite_update_rotation_axis', body.finite_update_rotation_axis)

    body.v = UTIL.vector3_from_xml(body_tag, 'v', body.v)
    body.w = UTIL.vector3_from_xml(body_tag, 'w', body.w)

    body.visual_material = UTIL.string_from_xml(body_tag, 'visual_material', body.visual_material)
    if body.visual_material is None:
        raise RuntimeError('__load_body(): visual material is missing')

    body.material = UTIL.string_from_xml(body_tag, 'material', body.material)
    if body.material is None:
        raise RuntimeWarning('__load_body(): material is missing using default instead')

    shape_ref = UTIL.string_from_xml(body_tag, 'shape', None)

    if shape_ref is None:
        raise RuntimeError('__load_body(): shape reference is missing')

    if shape_ref not in engine.shapes:
        raise RuntimeError('__load_body(): no such shape= ' + shape_ref + ' exist')

    body.shape = engine.shapes[shape_ref]

    for force_tag in body_tag.iter('force'):

        force_ref = UTIL.string_from_xml(force_tag, 'ref', None)

        if force_ref is None:
            raise RuntimeError('__load_body(): force reference attribute is missing')

        if force_ref not in engine.forces:
            raise RuntimeError('__load_body() no such force ' + force_ref + ' exist')

        force = engine.forces[force_ref]

        if force in body.forces:
            raise RuntimeError('__load_body() force' + force_ref + ' already added to body')

        body.forces.append(force)

    if body.is_scripted:

        if body.scripted_motion is not None:
            raise RuntimeError('__load_body(): Body' + name + ' already had a scripted motion')

        motion = KeyframeMotion()

        for keyframe_tag in body_tag.iter('keyframe'):
            time = UTIL.float_from_xml(keyframe_tag, 'time', None)
            r = UTIL.vector3_from_xml(keyframe_tag, 'r', None)
            q = UTIL.quaternion_from_xml(keyframe_tag, 'q', None)

            if time is None or r is None or q is None:
                raise RuntimeError('__load_body(): Keyframe missing time, r or q attribute on body = ' + name)

            motion.create_keyframe(time, r, q)

        body.scripted_motion = motion

        motion.keyframes.sort()

        body.r = motion.keyframes[0][1]
        body.q = motion.keyframes[0][2]

    engine.rigid_bodies[body.name] = body

    r = UTIL.vector3_from_xml(body_tag, 'r', body.r)
    q = UTIL.quaternion_from_xml(body_tag, 'q', body.q)
    use_model_frame = UTIL.bool_from_xml(body_tag, 'use_model_frame', True)
    set_position(engine, name, r, use_model_frame)
    set_orientation(engine, name, q, use_model_frame)

    # Note if density is used it will override mass and inertia
    #  settings. This means that mass properties are set based on shape
    #  information.
    mass = UTIL.float_from_xml(body_tag, 'mass', body.mass)
    inertia = UTIL.vector3_from_xml(body_tag, 'inertia', body.inertia)
    density = UTIL.float_from_xml(body_tag, 'density', None)

    if density is None:
        if mass > 0.0:
            set_mass(engine, name, mass)
        else:
            raise RuntimeError('__load_body(): Either missing mass or density attributes or mass attribute was '
                               'non-positive on body = ' + name)

        if inertia[0] > 0.0 and inertia[1] > 0.0 and inertia[2] > 0.0:
            set_inertia(engine, name, inertia)
        else:
            raise RuntimeError(
                '__load_body(): Either missing inertia or density attributes or inertia attribute has non-positive '
                'entries on body = ' + name)
    else:
        if density > 0.0:
            set_mass_properties_from_shape(engine, name, density)
        else:
            raise RuntimeError('__load_body(): Density attributes was '
                               'non-positive on body = ' + name)


def __save_body(body, engine_tag):
    body_tag = ET.SubElement(engine_tag, 'body')

    body_tag.attrib['name'] = body.name
    body_tag.attrib['active'] = str(body.is_active)
    body_tag.attrib['material'] = body.material

    type_value = ''
    if body.is_free:
        type_value = 'free'
    if body.is_scripted:
        type_value = 'scripted'
    if body.is_fixed:
        type_value = 'fixed'

    body_tag.attrib['type'] = type_value
    body_tag.attrib['use_finite_update'] = str(body.use_finite_update)

    if body.finite_update_rotation_axis is not None:
        body_tag.attrib['finite_update_rotation_axis'] = UTIL.array2string(body.finite_update_rotation_axis)

    body_tag.attrib['r'] = UTIL.array2string(body.r)
    body_tag.attrib['q'] = UTIL.array2string(body.q)
    body_tag.attrib['v'] = UTIL.array2string(body.v)
    body_tag.attrib['w'] = UTIL.array2string(body.w)
    body_tag.attrib['mass'] = str(body.mass)
    body_tag.attrib['inertia'] = UTIL.array2string(body.inertia)

    if body.shape is not None:
        body_tag.attrib['shape'] = body.shape.name

    body_tag.attrib['visual_material'] = body.visual_material

    for force in body.forces:
        force_tag = ET.SubElement(body_tag, 'force')
        force_tag.attrib['ref'] = force.name

    if body.scripted_motion is not None:
        for keyframe in body.scripted_motion.keyframes:
            key_tag = ET.SubElement(body_tag, 'keyframe')
            key_tag.attrib['time'] = str(keyframe[0])
            key_tag.attrib['r'] = UTIL.array2string(keyframe[1])
            key_tag.attrib['q'] = UTIL.array2string(keyframe[2])


def __load_restart_data(engine):
    if not engine.solver_params.load_restart:
        return

    filename = os.path.join(engine.solver_params.restart_path, engine.solver_params.restart_filename)

    ext = os.path.splitext(filename)[-1].lower()

    if ext != '.xml':
        raise RuntimeError('__load_restart_data(): file was not a xml file'+filename)

    if not os.path.isfile(filename):
        return

    if not os.path.exists(filename):
        return

    xml = ET.parse(filename)
    root = xml.getroot()

    engine.solver_params.current_time = UTIL.float_from_xml(root,'time', engine.solver_params.current_time)

    for body_tag in root.iter('body'):
        name = UTIL.string_from_xml(body_tag, 'name', None)
        if name is None:
            raise RuntimeError('__load_restart_data(): missing name')
        if name not in engine.rigid_bodies:
            raise RuntimeError('__load_restart_data(): body with that name did not exist')
        body = engine.rigid_bodies[name]
        if body is None:
            raise RuntimeError('__load_restart_data(): body was None')

        body.r = UTIL.vector3_from_xml(body_tag, 'r', body.r)
        body.q = UTIL.quaternion_from_xml(body_tag, 'q', body.q)
        body.v = UTIL.vector3_from_xml(body_tag, 'v', body.v)
        body.w = UTIL.vector3_from_xml(body_tag, 'w', body.w)

    print('Done loading restart data for time=', engine.solver_params.current_time)


def __load_connector(engine, body, connector_tag):
    if body is None:
        raise RuntimeError('__load_connector() body was None')

    # TODO Need to carefully design control over rigging frame

    connector = JointConnector(body)
    connector.transform.r = UTIL.vector3_from_xml(connector_tag, 'r', connector.transform.r)
    connector.transform.q = UTIL.quaternion_from_xml(connector_tag, 'q', None)

    if connector.transform.q is None:
        axis3 = UTIL.vector3_from_xml(connector_tag, 'axis', V3.k())
        axis1, axis2, axis3 = V3.make_orthonormal_vectors(axis3)
        R = M3.make_from_rows(axis1, axis2, axis3)
        connector.transform.q = Q.from_matrix(R)

    return connector


def __save_connector(socket, type, ball_joint_tag):
    if type in ['plug', 'socket']:
        connector_tag = ET.SubElement(ball_joint_tag, type)
    else:
        raise RuntimeError('save_connector() illegal type encountered must be plug or socket')

    connector_tag.attrib['r'] = UTIL.array2string(socket.transform.r)
    connector_tag.attrib['q'] = UTIL.array2string(socket.transform.q)


def __load_ball_joint(engine, ball_joint_tag):
    name = UTIL.string_from_xml(ball_joint_tag, 'name', None)

    if name is None:
        raise RuntimeError('__load_ball_joint() missing name attribute')

    if name in engine.joints:
        raise RuntimeError('__load_ball_joint(): ball joint with name =' + name + ' already exists')

    socket_body_name = UTIL.string_from_xml(ball_joint_tag, 'socket_body', None)

    if socket_body_name is None:
        raise RuntimeError('__load_ball_joint() missing socket_body attribute')

    if socket_body_name not in engine.rigid_bodies:
        raise RuntimeError('__load_ball_joint(): socket_body with name =' + socket_body_name + ' could not be found')

    socket_body = engine.rigid_bodies[socket_body_name]

    plug_body_name = UTIL.string_from_xml(ball_joint_tag, 'plug_body', None)

    if plug_body_name is None:
        raise RuntimeError('__load_ball_joint() missing plug_body attribute')

    if plug_body_name not in engine.rigid_bodies:
        raise RuntimeError('__load_ball_joint(): plug body with name =' + plug_body_name + ' could not be found')

    plug_body = engine.rigid_bodies[plug_body_name]

    joint = BallJoint(name, socket_body, plug_body)

    socket_tag = ball_joint_tag.find('socket')
    plug_tag = ball_joint_tag.find('plug')

    socket = __load_connector(engine, socket_body, socket_tag)
    plug = __load_connector(engine, plug_body, plug_tag)

    joint.socket = socket
    joint.plug = plug

    joint.error_reduction = UTIL.float_from_xml(ball_joint_tag, 'error_reduction', joint.error_reduction)

    engine.joints[joint.name] = joint


def __save_ball_joint(ball_joint, engine_tag):
    ball_joint_tag = ET.SubElement(engine_tag, 'ball')

    ball_joint_tag.attrib['name'] = ball_joint.name
    ball_joint_tag.attrib['socket_body'] = ball_joint.socket.body.name
    ball_joint_tag.attrib['plug_body'] = ball_joint.plug.body.name
    ball_joint_tag.attrib['error_reduction'] = str(ball_joint.error_reduction)
    __save_connector(ball_joint.socket, 'socket', ball_joint_tag)
    __save_connector(ball_joint.plug, 'plug', ball_joint_tag)


def load_from_elementtree(root):

    engine_tag = root.find('engine')
    if engine_tag is None:
        raise RuntimeError('load_from_elementtree(): no engine element')

    engine = Engine()

    __load_profiler(engine, engine_tag)
    __load_solver_parameters(engine, engine_tag)
    __load_material_library(engine, engine_tag)

    for gravity_tag in engine_tag.iter('gravity'):
        __load_gravity_force(engine, gravity_tag)

    for damping_tag in engine_tag.iter('damping'):
        __load_damping_force(engine, damping_tag)

    for shape_tag in engine_tag.iter('shape'):
        __load_shape(engine, shape_tag)

    for body_tag in engine_tag.iter('body'):
        __load_body(engine, body_tag)

    for ball_joint_tag in engine_tag.iter('ball'):
        __load_ball_joint(engine, ball_joint_tag)

    for python_tag in engine_tag.iter('python'):
        code = python_tag.text
        code = textwrap.dedent(code)
        try:
            exec(code)
        except Exception as e:
            print(e)

    __load_motion_recorder(engine, engine_tag)

    if engine.solver_params.mode == "simulate" and engine.motion_recorder.load and engine.motion_recorder.on:
        raise RuntimeWarning('Loading motion data and simulating and recording motion data')

    __load_restart_data(engine)

    return engine


def save_to_elementtree(engine, root):
    if engine is None:
        return

    engine_tag = ET.SubElement(root, 'engine')

    __save_profiler(engine, engine_tag)
    __save_solver_parameters(engine, engine_tag)
    __save_material_library(engine, engine_tag)

    for force in engine.forces.values():
        if force.force_type == 'Gravity':
            __save_gravity_force(force, engine_tag)
        if force.force_type == 'Damping':
            __save_damping_force(force, engine_tag)

    for shape in engine.shapes.values():
        __save_shape(shape, engine_tag)

    for body in engine.rigid_bodies.values():
        __save_body(body, engine_tag)

    for joint in engine.joints.values():
        if joint.type == 'ball':
            __save_ball_joint(joint, engine_tag)

    __save_motion_recorder(engine, engine_tag)


def load(filename):
    ext = os.path.splitext(filename)[-1].lower()

    if ext != '.xml':
        raise RuntimeError('load: file was not a xml file'+filename)

    if not os.path.isfile(filename):
        raise RuntimeError('load: filename was not a file'+filename)

    if not os.path.exists(filename):
        raise RuntimeError('load: file did not exist'+filename)

    xml = ET.parse(filename)
    root = xml.getroot()

    return load_from_elementtree(root)


def save(engine, filename):
    ext = os.path.splitext(filename)[-1].lower()

    if ext != '.xml':
        raise RuntimeError('save: file was not a xml file'+filename)

    root = ET.Element('scene')

    save_to_elementtree(engine, root)

    UTIL.xml_pretty_indent(root)

    tree = ET.ElementTree(root)
    tree.write(filename)
