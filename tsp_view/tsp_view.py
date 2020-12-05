#!/usr/bin/python

import os
import sys
import subprocess


def main():
    # number of nodes (first line)
    number_of_nodes = int(sys.stdin.readline())

    # coordinates of nodes
    nodes = [{'x': 0.0, 'y': 0.0} for _ in range(number_of_nodes)]
    for i in range(number_of_nodes):
        line = sys.stdin.readline()
        x_str, y_str = line.split(' ')
        nodes[i]['x'] = float(x_str)
        nodes[i]['y'] = float(y_str)

    # route
    route = []
    while True:
        line = sys.stdin.readline()
        if not line:
            break
        else:
            route.append(int(line[:-1]))
    route.append(route[0])

    # generate plot data file
    file_vertex = 'data/points.dat'
    file_edge = 'data/lines.dat'

    with open(file_vertex, 'w') as f:
        for node in nodes:
            f.write('%f %f\n' % (node['x'], node['y']))

    with open(file_edge, 'w') as f:
        for i in route:
            f.write('%f %f\n' % (nodes[i]['x'], nodes[i]['y']))

    # execute gnuplot
    subprocess.call(['gnuplot', './plot.gp'])


if __name__ == '__main__':
    main()
