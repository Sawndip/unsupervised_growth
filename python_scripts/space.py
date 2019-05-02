# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 13:25:50 2017

@author: jingroup

Script defines distances and relationships between spatial arrangement and number
of coordinates used to define neuronal location in space
"""
import numpy as np
from functools import partial


def distance_on_sphere(v1, v2, R):
    """
    Computes distance between two points on sphere of radius R
    
    Input: v1,v2 - 3d vectors containing x,y,z coordinates; 
           R - sphere radius
           
    Output: distance between two points on sphere
    """
    return np.arccos(np.dot(v1, v2)) / R


def distance(v1, v2):
    """
    Computes distance between two points in space
    
    Input: v1,v2 - vectors containing coordinates; 
           
    Output: distance between two points
    """
    return np.linalg.norm(v1-v2)

num_coordinates = {"square" : 2, "sphere" : 3, "cube" : 3}
distance_function = {"sphere" : partial(distance_on_sphere, R=1), "square" : distance, "cube" : distance}

if __name__ == "__main__":
    print num_coordinates["square"]
    print num_coordinates["cube"]
    print num_coordinates["sphere"]
    
    x = np.array([-1,0,0])
    y = np.array([1,0,0])
    
    print distance_function["cube"](x,y)
    

