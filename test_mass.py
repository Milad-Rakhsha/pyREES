import REESMass.mass as MASS
import REESMesh.mesh as MESH
import REESMath.coordsys as C
import REESMath.vector3 as V3
import REESMath.quaternion as Q


if __name__ == '__main__':
    mesh = MESH.read_obj('resources/objs/box.obj')

    prop = MASS.compute_mass_properties(mesh, 1.0)
    print(prop)

    print('--- Original Coordinates ----------------------------')
    for vh in mesh.vertices():
        print(MESH.get_vertex_coords(mesh, vh))
    print('--- Transformed Coordinates -------------------------')
    MESH.translate(mesh, V3.i())
    for vh in mesh.vertices():
        print(MESH.get_vertex_coords(mesh, vh))
    print('-----------------------------------------------------')


    prop = MASS.compute_mass_properties(mesh, 1.0)
    print(prop)

    (r, q, m, I_body) = MASS.xform_2_body_space(prop)
    print('model to body translation', r)
    print('model to body rotation', q)

    MESH.translate(mesh, r)
    MESH.rotate(mesh, q)

    print('--- Body Coordinates -------------------------------')
    for vh in mesh.vertices():
        print(MESH.get_vertex_coords(mesh, vh))
    print('----------------------------------------------------')
