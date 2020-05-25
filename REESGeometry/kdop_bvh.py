import REESGeometry.kdop as kdop
import REESMath.vector3 as V3


class ContactPoint:
    def __init__(self):
        self.p = V3.zero()
        self.n = V3.zero()
        self.d = 0.0


class Node:

    UNDEFINED = -1

    def __init__(self):
        self.volume = kdop.KDOP()
        self.parent = Node.UNDEFINED
        self.start = Node.UNDEFINED
        self.end = Node.UNDEFINED


class SubTree:
    def __init__(self):
        self.nodes = []


class Tree:
    def __init__(self, mesh):
        self.root = kdop.KDOP()
        self.branches = []
        self.mesh = mesh


def __is_root(node):
    return node.parent == Node.UNDEFINED and node.start!=Node.UNDEFINED and node.end!=Node.UNDEFINED


def __is_leaf(node):
    return node.start!=Node.UNDEFINED and node.start==node.end


def __is_undefined(node):
    return node.parent == Node.UNDEFINED and node.start==Node.UNDEFINED and node.end==Node.UNDEFINED


def __subtree_traversal(mesh_A, mesh_B, subtree_A, subtree_B, idx_A, idx_B, results):
    node_A = subtree_A.nodes[idx_A]
    node_B = subtree_B.nodes[idx_B]

    if not kdop.overlap(node_A.volume, node_B.volume):
        return

    A_is_leaf = __is_leaf(node_A)
    B_is_leaf = __is_leaf(node_B)

    if A_is_leaf and B_is_leaf:
        fah = mesh_A.face_handle(node_A.start)
        fbh = mesh_A.face_handle(node_B.start)
        # test actual geometry for overlap and generate contact points

        results.append(ContactPoint())

    elif not A_is_leaf and not B_is_leaf:

        for child_A in range(node_A.start, node_A.end + 1):
            for child_B in range(node_B.start, node_B.end + 1):
                traversal(mesh_A, mesh_B, subtree_A, subtree_B, child_A, child_B)

    elif not A_is_leaf and B_is_leaf:

        for child_A in range(node_A.start, node_A.end + 1):
            traversal(mesh_A, mesh_B, subtree_A, subtree_B, child_A, idx_B)

    elif A_is_leaf and not B_is_leaf:

        for child_B in range(node_B.start, node_B.end + 1):
            traversal(mesh_A, mesh_B, subtree_A, subtree_B, idx_A, child_B)


def traversal(tree_A, tree_B, results):

    if not kdop.overlap(tree_A.root, tree_B.root):
        return

    for subtree_A in tree_A.branches:
        for subtree_B in tree_B.branches:
            __subtree_traversal(tree_A.mesh, tree_B.mesh, subtree_A, subtree_B, 0, 0, results)


def __refit_subtree(X, mesh, subtree):
    for node in reversed(subtree.nodes):
        if __is_undefined(node):
            continue
        if __is_leaf(node):
            fh = mesh.face_handle(node.start)
            node.volume = kdop.make() # add geometry arguments
        else:
            left_child = subtree.nodes[node.start]
            right_child = subtree.nodes[node.end]
            node.volume = kdop.union(left_child.volume, right_child.volume)


def refit_tree(X, tree):

    for subtree in tree.branches:
        ___refit_subtree(X, tree.mesh, subtree)

    tree.root kdop.KDOP()

    for subtree in tree.branches:
        tree.root = kdop.union(tree.root, subtree.nodes[0].volume)


def make_tree(mesh):
    return 0


