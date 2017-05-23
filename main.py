import sys
import getopt
import REESSimulation.api as API
import REESSimulation.xml as XML
from REESVisualization.main_window import *


def no_gui(argv):
    input_file = None
    output_file = None
    try:
        opts, args = getopt.getopt(argv, 'hi:o:', ['help', 'input=', 'output='])
    except getopt.GetoptError:
        print('Usage: main.py -i <input file> -o <output file>')
        sys.exit(2)

    for opt, arg in opts:
        if opt in ('-h', '--help'):
            print('Usage: main.py -i <input file> -o <output file>')
            sys.exit()
        elif opt in ('-i', '--input'):
            input_file = arg
        elif opt in ('-o', '--output'):
            output_file = arg

    if input_file is None:
        print('Error input file was None')
        print('Usage: main.py -i <input file> -o <output file>')
        sys.exit(2)

    if output_file is None:
        print('Error output file was None')
        print('Usage: main.py -i <input file> -o <output file>')
        sys.exit(2)

    print('input file =', input_file)
    print('output file =', output_file)

    engine = XML.load(input_file)
    while API.simulate(engine):
        print('.')
    XML.save(engine, output_file)


def with_gui():
    initialize_opengl()

    app = QApplication(sys.argv)
    main_window = MainWindow()
    sys.exit(app.exec_())


if __name__ == '__main__':
    if len(sys.argv[1:]) > 0:
        no_gui(sys.argv[1:])
    else:
        with_gui()


"""

How to pass new positions to KDOPBVH?

Support for creating a BVH tandem traversal using CCD (self-made)

Support for making procedural configurations? (self-made)

Support for making a Prox-solver (self-made)

"""


'''

G = nx.Graph()

G.add_node(0, body=RigidBody('a'))
G.add_node(1, body=RigidBody('b'))
G.add_node(2, body=RigidBody('c'))
G.add_node(3, body=RigidBody('d'))
G.add_node(4, body=RigidBody('e'))
G.add_edge(0, 1, info=Connection())
G.add_edge(0, 2, info=Connection())
G.add_edge(0, 3, info=Connection())
G.add_edge(2, 4, info=Connection())
G.add_edge(3, 4, info=Connection())

print(G.nodes())
print(G.edges())

print('----------------------')
for n in G.neighbors_iter(0):
    print(n)

print('----------------------')
for n in G.neighbors_iter(4):
    body = G.node[n]['body']
    print(body.name)
    body.name *= 2

print('----------------------')
for n in G.neighbors_iter(4):
    body = G.node[n]['body']
    print(body.name)

print('----------------------')
for n in G.neighbors_iter(0):
    info = G.get_edge_data(0, n)['info']
    print(info.value)


print('----------------------')
if G.has_edge(0,1):
    print('has edge 0-1')


print('----------------------')
if not G.has_edge(0,4):
    print('did not have edge 0-4')

nx.draw(G)
plt.show()
'''

