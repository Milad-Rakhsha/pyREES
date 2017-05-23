import xml.etree.ElementTree as ET
import REESUtility.util as UTIL
import REESMath.vector3 as V3
import REESMath.angle as ANGLE
from REESVisualization.types import *


def __load_movie_recorder(visualization_tag):
    movie_recorder_tag = visualization_tag.find('movie_recorder')

    mov = MovieRecorder()

    if movie_recorder_tag is None:
        return mov

    mov.path = UTIL.string_from_xml(movie_recorder_tag, 'path', mov.path)
    mov.base_name = UTIL.string_from_xml(movie_recorder_tag, 'basename', mov.base_name)
    mov.image_type = UTIL.string_from_xml(movie_recorder_tag, 'type', mov.type)
    mov.set_state(UTIL.bool_from_xml(movie_recorder_tag, 'on', mov.on))

    return mov


def __save_movie_recorder(movie_recorder, visualization_tag):
    mov_tag = ET.SubElement(visualization_tag, 'movie_recorder')
    mov_tag.attrib['path'] = movie_recorder.path
    mov_tag.attrib['basename'] = movie_recorder.base_name
    mov_tag.attrib['type'] = movie_recorder.image_type
    mov_tag.attrib['on'] = str(movie_recorder.on)


def __load_camera(visualization_tag):
    camera = Camera()

    camera_tag = visualization_tag.find('camera')

    if camera_tag is None:
        return camera

    camera.z_near = UTIL.float_from_xml(camera_tag, 'near', camera.z_near)
    camera.z_far = UTIL.float_from_xml(camera_tag, 'far', camera.z_far)
    fovy_in_degrees = UTIL.float_from_xml(camera_tag, 'fovy', ANGLE.radians_to_degrees(camera.fovy))
    camera.fovy = ANGLE.degrees_to_radians(fovy_in_degrees)
    eye = UTIL.vector3_from_xml(camera_tag, 'eye', camera.eye)
    up = UTIL.vector3_from_xml(camera_tag, 'up', camera.up)
    center = UTIL.vector3_from_xml(camera_tag, 'center', camera.center)
    camera.update(eye, center, up)

    return camera


def __save_camera(camera, visualization_tag):
    cam_tag = ET.SubElement(visualization_tag, 'camera')
    cam_tag.attrib['eye'] = UTIL.array2string(camera.eye)
    cam_tag.attrib['up'] = UTIL.array2string(camera.up)
    cam_tag.attrib['center'] = UTIL.array2string(camera.center)
    fovy_degrees = ANGLE.radians_to_degrees(camera.fovy)
    cam_tag.attrib['fovy'] = str(fovy_degrees)
    cam_tag.attrib['near'] = str(camera.z_near)
    cam_tag.attrib['far'] = str(camera.z_far)


def __load_light(light_setup, light_tag):
    if light_tag is None:
        raise RuntimeError('__load_light can not load from none?')

    index = UTIL.int_from_xml(light_tag, 'idx', None)

    if index is None:
        raise RuntimeError('__load_light missing index attribute')

    light = SpotLight(index)

    light.position = UTIL.vector3_from_xml(light_tag, 'position', light.position)
    light.target = UTIL.vector3_from_xml(light_tag, 'target', light.target)
    light.cutoff_angle = UTIL.float_from_xml(light_tag, 'cutoff_angle', light.cutoff_angle)
    light.attenuation = UTIL.float_from_xml(light_tag, 'attenuation', light.attenuation)
    light.Is = UTIL.vector3_from_xml(light_tag, 'Is', light.Is)
    light.Id = UTIL.vector3_from_xml(light_tag, 'Id', light.Id)
    light.Ia = UTIL.vector3_from_xml(light_tag, 'Ia', light.Ia)

    light_setup.lights.append(light)


def __save_light(light, visualization_tag):
    light_tag = ET.SubElement(visualization_tag, 'light')
    light_tag.attrib['idx'] = str(light.index)
    light_tag.attrib['position'] = UTIL.array2string(light.position)
    light_tag.attrib['target'] = UTIL.array2string(light.target)
    light_tag.attrib['cutoff_angle'] = str(light.cutoff_angle)
    light_tag.attrib['attenuation'] = str(light.attenuation)
    light_tag.attrib['Is'] = UTIL.array2string(light.Is)
    light_tag.attrib['Id'] = UTIL.array2string(light.Id)
    light_tag.attrib['Ia'] = UTIL.array2string(light.Ia)


def __load_visual_material(materials, visual_material_tag):
    if visual_material_tag is None:
        raise RuntimeError('__load_material can not load from none?')

    name = UTIL.string_from_xml(visual_material_tag, 'name', None)

    if name is None:
        raise RuntimeError('__load_material missing name on material')

    material = VisualMaterial(name)

    material.wire_color = UTIL.vector3_from_xml(visual_material_tag, 'wire_color', material.wire_color)
    material.wire_thickness = UTIL.float_from_xml(visual_material_tag, 'wire_thickness', material.wire_thickness)
    material.specular_exponent = UTIL.float_from_xml(visual_material_tag, 'specular_exponent', material.specular_exponent)
    material.Ks = UTIL.vector3_from_xml(visual_material_tag, 'Ks', material.Ks)
    material.Kd = UTIL.vector3_from_xml(visual_material_tag, 'Kd', material.Kd)
    material.Ka = UTIL.vector3_from_xml(visual_material_tag, 'Ka', material.Ka)

    materials[material.name] = material


def __save_visual_material(visual_material, visualization_tag):
    material_tag = ET.SubElement(visualization_tag, 'visual_material')
    material_tag.attrib['name'] = visual_material.name
    material_tag.attrib['Ks'] = UTIL.array2string(visual_material.Ks)
    material_tag.attrib['Kd'] = UTIL.array2string(visual_material.Kd)
    material_tag.attrib['Ka'] = UTIL.array2string(visual_material.Ka)
    material_tag.attrib['specular_exponent'] = str(visual_material.specular_exponent)
    material_tag.attrib['wire_color'] = UTIL.array2string(visual_material.wire_color)
    material_tag.attrib['wire_thickness'] = str(visual_material.wire_thickness)


def __load_grid(visualization_tag):

    grid_tag = visualization_tag.find('grid')

    grid = Grid()

    if grid_tag is None:
        return grid

    grid.start_x = UTIL.float_from_xml(grid_tag, 'start_x', grid.start_x)
    grid.end_x = UTIL.float_from_xml(grid_tag, 'end_x',  grid.end_x)
    grid.cells_x = UTIL.int_from_xml(grid_tag, 'cells_x', grid.cells_x)
    grid.start_z = UTIL.float_from_xml(grid_tag, 'start_z', grid.start_z)
    grid.end_z = UTIL.float_from_xml(grid_tag, 'end_z',  grid.end_z)
    grid.cells_z = UTIL.int_from_xml(grid_tag, 'cells_z', grid.cells_z)
    grid.on = UTIL.bool_from_xml(grid_tag, 'on', grid.on)

    if grid.cells_z <= 0:
        raise RuntimeError('load_grid_from_xml(): illegal cells z value')
    if grid.cells_x <= 0:
        raise RuntimeError('load_grid_from_xml(): illegal cells x value')

    return grid


def __save_grid(grid, visualization_tag):
    grid_tag = ET.SubElement(visualization_tag, 'grid')
    grid_tag.attrib['start_x'] = str(grid.start_x)
    grid_tag.attrib['end_x'] = str(grid.end_x)
    grid_tag.attrib['cells_x'] = str(grid.cells_x)
    grid_tag.attrib['start_z'] = str(grid.start_z)
    grid_tag.attrib['end_z'] = str(grid.end_z)
    grid_tag.attrib['cells_z'] = str(grid.cells_z)
    grid_tag.attrib['on'] = str(grid.on)


def load_from_elementtree(root):

    visualization_tag = root.find('visualization')
    if visualization_tag is None:
        raise RuntimeError('load(): no visualization element')

    light_setup = LightSetup()
    for light_tag in visualization_tag.iter('light'):
        __load_light(light_setup, light_tag)

    materials = dict()
    for visual_material_tag in visualization_tag.iter('visual_material'):
        __load_visual_material(materials, visual_material_tag)

    movie_recorder = __load_movie_recorder(visualization_tag)
    camera = __load_camera(visualization_tag)
    grid = __load_grid(visualization_tag)

    clear_color = UTIL.vector3_from_xml(visualization_tag, 'clear_color', V3.make(0.3, 0.3, 0.3))
    contact_scale = UTIL.float_from_xml(visualization_tag, 'contact_scale', 0.5)

    return camera, light_setup, materials, clear_color, grid, movie_recorder, contact_scale


def save_to_elementtree(camera, light_setup, materials, clear_color, grid, movie_recorder, contact_scale, root):
    visualization_tag = ET.SubElement(root, 'visualization')
    visualization_tag.attrib['clear_color'] = UTIL.array2string(clear_color)
    visualization_tag.attrib['contact_scale'] = str(contact_scale)
    for visual_material in materials.values():
        __save_visual_material(visual_material, visualization_tag)
    for light in light_setup.lights:
        __save_light(light, visualization_tag)
    __save_camera(camera, visualization_tag)
    __save_movie_recorder(movie_recorder, visualization_tag)
    __save_grid(grid, visualization_tag)


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


def save(camera, light_setup, materials, clear_color, grid_info, movie_recorder, contact_scale, filename):
    ext = os.path.splitext(filename)[-1].lower()

    if ext != '.xml':
        raise RuntimeError('save: file was not a xml file'+filename)

    root = ET.Element('scene')
    save_to_elementtree(camera, light_setup, materials, clear_color, grid_info, movie_recorder, contact_scale, root)
    UTIL.xml_pretty_indent(root)
    tree = ET.ElementTree(root)
    tree.write(filename)

