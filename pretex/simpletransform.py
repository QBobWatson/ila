'''
Copyright (C) 2006 Jean-Francois Barraud, barraud@math.univ-lille1.fr
Copyright (C) 2010 Alvin Penner, penner@vaxxine.com

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
barraud@math.univ-lille1.fr

This code defines several functions to make handling of transform
attribute easier.
'''

# JDR: copied from inkscape extensions

import math
import re

def parse_transform(transf, mat=None):
    if mat is None:
        mat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]
    if transf == "" or transf is None:
        return mat
    stransf = transf.strip()
    result = re.match(r"(translate|scale|rotate|skewX|skewY|matrix)"
                      r"\s*\(([^)]*)\)\s*,?", stransf)
    #-- translate --
    if result.group(1) == "translate":
        args = result.group(2).replace(',', ' ').split()
        dx = float(args[0])
        if len(args) == 1:
            dy = 0.0
        else:
            dy = float(args[1])
        matrix = [[1, 0, dx], [0, 1, dy]]
    #-- scale --
    if result.group(1) == "scale":
        args = result.group(2).replace(',', ' ').split()
        sx = float(args[0])
        if len(args) == 1:
            sy = sx
        else:
            sy = float(args[1])
        matrix = [[sx, 0, 0], [0, sy, 0]]
    #-- rotate --
    if result.group(1) == "rotate":
        args = result.group(2).replace(',', ' ').split()
        a = float(args[0])*math.pi/180
        if len(args) == 1:
            cx, cy = (0.0, 0.0)
        else:
            cx, cy = map(float, args[1:])
        matrix = [[math.cos(a), -math.sin(a), cx],
                  [math.sin(a),  math.cos(a), cy]]
        matrix = compose_transform(matrix, [[1, 0, -cx], [0, 1, -cy]])
    #-- skewX --
    if result.group(1) == "skewX":
        a = float(result.group(2))*math.pi/180
        matrix = [[1, math.tan(a), 0], [0, 1, 0]]
    #-- skewY --
    if result.group(1) == "skewY":
        a = float(result.group(2))*math.pi/180
        matrix = [[1, 0, 0], [math.tan(a), 1, 0]]
    #-- matrix --
    if result.group(1) == "matrix":
        a11, a21, a12, a22, v1, v2 = result.group(2).replace(',', ' ').split()
        matrix = [[float(a11), float(a12), float(v1)],
                  [float(a21), float(a22), float(v2)]]

    matrix = compose_transform(mat, matrix)
    if result.end() < len(stransf):
        return parse_transform(stransf[result.end():], matrix)
    return matrix

def format_transform(mat):
    return "matrix({},{},{},{},{},{})".format(
        mat[0][0], mat[1][0], mat[0][1], mat[1][1], mat[0][2], mat[1][2])

def invert_transform(mat):
    det = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0]
    if det !=0:  # det is 0 only in case of 0 scaling
        # invert the rotation/scaling part
        a11 =  mat[1][1]/det
        a12 = -mat[0][1]/det
        a21 = -mat[1][0]/det
        a22 =  mat[0][0]/det
        # invert the translational part
        a13 = -(a11*mat[0][2] + a12*mat[1][2])
        a23 = -(a21*mat[0][2] + a22*mat[1][2])
        return [[a11,a12,a13],[a21,a22,a23]]
    else:
        return[[0,0,-mat[0][2]],[0,0,-mat[1][2]]]

def compose_transform(M1, M2):
    a11 = M1[0][0]*M2[0][0] + M1[0][1]*M2[1][0]
    a12 = M1[0][0]*M2[0][1] + M1[0][1]*M2[1][1]
    a21 = M1[1][0]*M2[0][0] + M1[1][1]*M2[1][0]
    a22 = M1[1][0]*M2[0][1] + M1[1][1]*M2[1][1]

    v1 = M1[0][0]*M2[0][2] + M1[0][1]*M2[1][2] + M1[0][2]
    v2 = M1[1][0]*M2[0][2] + M1[1][1]*M2[1][2] + M1[1][2]
    return [[a11,a12,v1],[a21,a22,v2]]

