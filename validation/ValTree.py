from graphviz import Digraph

def create_val_dec_tree(out_file):

    #nodes
    #edited to add "false" for nodes that should never be reached
    binary_tree_attrib = ['both_primers_found?', 'same_contig?', 'ehybrid ?', 'ehybrid.?', 'ehybrid?', 'cutoff by contig break ?', 'CASE 6',
                          'epcr?', 'CASE 8', 'cutoff by contig break?']
    binary_tree_attrib.insert(10, 'CASE 10')
    binary_tree_attrib.insert(11, 'False1')
    binary_tree_attrib.insert(12, 'CASE 12')
    for i in range(13, 15):
        binary_tree_attrib.insert(i, 'Fail')
    binary_tree_attrib.insert(15, 'False2')
    binary_tree_attrib.insert(16, 'primers facing each other?')
    for i in range(17, 19):
        binary_tree_attrib.insert(i, 'Fail')
    binary_tree_attrib.insert(19, 'False3')
    binary_tree_attrib.insert(20, 'CASE 20')
    for i in range(21, 23):
        binary_tree_attrib.insert(i, 'Fail')
    for i in range(23, 33):
        binary_tree_attrib.insert(i, 'Fail')
    binary_tree_attrib.insert(33, 'SNP?')
    binary_tree_attrib.insert(34, 'CASE 34')
    for i in range(35, 67):
        binary_tree_attrib.insert(i, 'Fail')
    binary_tree_attrib.insert(67, 'CASE 67')
    binary_tree_attrib.insert(68, 'correct distance btwn primers?')
    for i in range(69, 137):
        binary_tree_attrib.insert(i, 'Fail')
    binary_tree_attrib.insert(137, 'False4')
    for i in range(138, 139):
        binary_tree_attrib.insert(i, 'CASE 138')

    #edges
    edges_lables = get_edges(binary_tree_attrib)
    edges = edges_lables[0]
    labels = edges_lables[1]

    dot = Digraph(comment='F- Validation Decision Tree')
    for node in binary_tree_attrib:
        dot.node(node)
    for i, edge in enumerate(edges):
        dot.edge(edge[0], edge[1], labels[i])
    print(dot.source)
    dot.render(out_file, view=True)

def get_edges(attr_bin_tree):
    edges = []
    labels = []
    for i in range(0, len(attr_bin_tree) - 1):
        if attr_bin_tree[i] != 'Fail' and 'False' not in attr_bin_tree[i] and 'CASE' not in attr_bin_tree[i]:
            if attr_bin_tree[(2 * i) + 1] != 'Fail':
                attr_i = attr_bin_tree[i]
                left = attr_bin_tree[(2 * i) + 1]
                edges.append((attr_i, left))
                labels.append("True")

            if attr_bin_tree[(2 * i) + 2] != 'Fail':
                attr_i = attr_bin_tree[i]
                right = attr_bin_tree[(2 * i) + 2]
                edges.append((attr_i, right))
                labels.append("False")
    return [edges, labels]

if __name__ == '__main__':
    out_file = '/home/sfisher/eCGF/f_neg_decision_tree.gv'
    create_val_dec_tree(out_file)