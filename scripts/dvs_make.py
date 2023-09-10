import math
from matplotlib import pyplot as plt

N_layer = 3
N_node = 15

NODES = []
CELLS = []


class Node:

    def __init__(self, nid, r, a):
        self.nid = nid
        self.r = r
        self.a = a

    def coordinate_x(self):
        return self.r * math.cos(self.a)

    def coordinate_y(self):
        return self.r * math.sin(self.a)


QUAD = 2
TRIA = 3


class Cell:

    def __init__(self, cid, elem_type, *args):
        self.cid = cid
        self.elem_type = elem_type
        self.node_list = args


def make_layer(radius, dr, node_num, offset, is_first_layer=False):
    in_layer = []
    out_layer = []
    if is_first_layer:
        nid = len(NODES)
        node = Node(nid, 0.0, 0.0)
        in_layer.append(nid)
        NODES.append(node)
    else:
        for i in range(node_num):
            nid = len(NODES)
            node = Node(nid, radius, 2 * (i / node_num) * math.pi + offset)
            in_layer.append(nid)
            NODES.append(node)
    for i in range(node_num):
        nid = len(NODES)
        node = Node(nid, radius + dr, 2 * (i / node_num) * math.pi + offset)
        out_layer.append(nid)
        NODES.append(node)
    for i in range(node_num):
        cid = len(CELLS)
        if is_first_layer:
            cell = Cell(cid, TRIA, in_layer[0], out_layer[i - 1], out_layer[i])
        else:
            cell = Cell(cid, QUAD, in_layer[i - 1], out_layer[i - 1], out_layer[i], in_layer[i])
        CELLS.append(cell)


def plot():
    for cell in CELLS:
        node_x = [NODES[i].coordinate_x() for i in cell.node_list]
        node_y = [NODES[i].coordinate_y() for i in cell.node_list]
        node_x.append(node_x[0])
        node_y.append(node_y[0])
        plt.plot(node_x, node_y, color="black", lw=0.4)


def write_to_neu():
    with open("dvs_make.neu", "w") as f:
        # control info
        f.write("CONTROL INFO\n"
                "     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL\n"
                f"     {len(NODES)}    {len(CELLS)}    1         2         2         2\n"
                "ENDOFSECTION\n")
        # node
        f.write("NODAL COORDINATES\n")
        for node in NODES:
            f.write(f"  {node.nid + 1} {node.coordinate_x()} {node.coordinate_y()}\n")
        f.write("ENDOFSECTION\n")
        # ELEMENTS/CELLS
        f.write("ELEMENTS/CELLS\n")
        for cell in CELLS:
            f.write(f"  {cell.cid + 1}  {cell.elem_type}  {len(cell.node_list)}")
            for n in cell.node_list:
                f.write(f"  {n + 1}")
            f.write("\n")
        f.write("ENDOFSECTION\n")
    print("write file.")


if __name__ == '__main__':
    r = 0.005
    dr = 0.05
    node_0 = 31
    for layer in range(50):
        node_num = int(node_0 + 6.0 * layer)
        make_layer(r, dr, node_num, 6.5 * math.pi / node_num * layer, layer == 0)
        r += dr
        dr += math.log(1.0 + layer * 0.0025 / 10)
    print("Cell num = ", len(CELLS))
    plot()
    plt.gca().set_aspect('equal')
    plt.show()
    write_to_neu()
