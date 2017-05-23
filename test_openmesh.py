import sys

sys.path.append('/usr/local/lib/python/')


import openmesh as OM
import REESMesh.mesh as MESH


if __name__ == '__main__':

    mesh = OM.TriMesh()

    vh0 = mesh.add_vertex(OM.TriMesh.Point(0, 1, 0))
    #vh1 = mesh.add_vertex(OM.TriMesh.Point(1, 0, 0))
    #vh2 = mesh.add_vertex(OM.TriMesh.Point(2, 1, 0))
    #vh3 = mesh.add_vertex(OM.TriMesh.Point(0, -1, 0))
    #vh4 = mesh.add_vertex(OM.TriMesh.Point(2, -1, 0))
    #fh0 = mesh.add_face(vh0, vh1, vh2)
    #fh1 = mesh.add_face(vh1, vh3, vh4)
    #fh2 = mesh.add_face(vh0, vh3, vh1)
    #mesh.garbage_collection()
    #mesh.delete_face(fh1)
    mesh.request_face_status()
    mesh.request_vertex_status()
    mesh.request_halfedge_status()
    mesh.delete_vertex(vh0)
    mesh.garbage_collection()

    mesh = MESH.read_obj('resources/objs/box.obj')

    print('vertices =', mesh.n_vertices())
    print('edges =', mesh.n_edges())
    print('half edges =', mesh.n_halfedges())
    print('faces =', mesh.n_faces())
