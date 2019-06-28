import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np


def parse_rsml(organ :ET, polylines :list, properties :dict, functions :dict) -> (list, dict, dict):
    """Recursivly parses the rsml file """

    for poly in organ.iterfind('geometry'):  # only one
        polyline = []
        for p in poly[0]:  # 0 is the polyline
            n = p.attrib
            newnode = [ float(n['x']), float(n['y']), float(n['z']) ]
            polyline.append(newnode)
        polylines.append(polyline)

    for prop in organ.iterfind('properties'):
        for p in prop:  # i.e legnth, type, etc..
            try:
                properties[str(p.tag)].append(float(p.attrib['value']))
            except:
                properties[str(p.tag)] = [float(p.attrib['value'])]

    for funcs in organ.iterfind('functions'):
        for fun in funcs:
            samples = [ ]
            for sample in fun.iterfind('sample'):
                samples.append(float(sample.attrib['value']))
            try:
                functions[str(fun.attrib['name'])].append(samples)
            except:
                functions[str(fun.attrib['name'])] = [samples]

    for elem in organ.iterfind('root'):  # and all laterals
        polylines, properties, functions = parse_rsml(elem, polylines, properties, functions)

    return polylines, properties, functions


def read_rsml(name :str) -> (list, dict, dict):
    """Parses the RSML file into: 
    a list of polylines, with one polyline per root
    a dictionary of properties one per root
    a dictionary of functions     
    """
    root = ET.parse(name).getroot()
    plant = root[1][0]
    polylines = []
    properties = { }
    functions = { }
    for elem in plant.iterfind('root'):
        (polylines, properties, functions) = parse_rsml(elem, polylines, properties, functions)

    return polylines, properties, functions


def plot_rsml(polylines :list):
    """Plots the polylines in y-z axis"""
    for i, pl in enumerate(polylines):
        nodes = np.array(pl)
        plt.plot(nodes[:, 1], nodes[:, 2])  # y,z plot
    plt.show()


if __name__ == '__main__':
        polylines, properties, functions = read_rsml("root_grid/RootSystem.rsml")
        print("Properties:")
        for key, v in properties.items() :
            print("\t", key)
        print("Functions:")
        for key, v in functions.items() :
            print("\t", key)
        plot_rsml(polylines)
