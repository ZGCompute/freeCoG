import scipy.io as io
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import os
import numpy as np
import sys

def plot_cortex(subj):
    
    # set up shader
    shader = gl.shaders.ShaderProgram('my_shader', [   
        gl.shaders.VertexShader("""
            varying vec3 normal;
            void main() {
                // compute here for use in fragment shader
                normal = normalize(gl_NormalMatrix * gl_Normal * 0.5);
                gl_FrontColor = gl_Color;
                gl_BackColor = gl_Color;
                gl_Position = ftransform();
            }
        """),
        gl.shaders.FragmentShader("""
            varying vec3 normal;
            void main() {
                float p = dot(normal, normalize(vec3(1.0, -4.0, -4.0)));
                p = p < 0. ? p * -0.8 : p * 0.8;
                vec4 color = gl_Color;
                color.x = color.x * (0.3 + p);
                color.y = color.y * (0.3 + p);
                color.z = color.z * (0.3 + p);
                gl_FragColor = color;
            }
        """)
    ])
    
    root = '/Users/zack/Documents/Imaging/subjects/';
    subj_root = root + subj + '/Meshes';
    os.chdir(subj_root);
    
    pg.setConfigOption('background', 'w');

    lh_tri = io.loadmat('lh_tri');
    lh_faces = lh_tri.get('tri');
    lh_faces = lh_faces -1;
    lh_vert = io.loadmat('lh_vert')
    lh_verts = lh_vert.get('vert')
    
    rh_tri = io.loadmat('rh_tri');
    rh_faces = rh_tri.get('tri');
    rh_faces = rh_faces -1;
    rh_vert = io.loadmat('rh_vert')
    rh_verts = rh_vert.get('vert')
    
    app = QtGui.QApplication([])
    w = gl.GLViewWidget()
    w.show()
    w.setWindowTitle(('Patient Pial Surface and Depths for: ' + subj))
    w.setCameraPosition(distance=250);
    
    
    lhem = gl.GLMeshItem(vertexes = lh_verts, faces=lh_faces, drawFaces=True, drawEdges=False, faceColor=(1,0,0,1),edgeColor=(0,0,0,0), smooth=True, shader='my_shader');
    lhem.translate(5, 5, 0)
    lhem.setGLOptions('additive');
    w.addItem(lhem);
    
    rhem = gl.GLMeshItem(vertexes = rh_verts, faces=rh_faces, drawFaces=True, drawEdges=False, faceColor=(1,0,0,1),edgeColor=(0,0,0,0), smooth=True, shader='my_shader');
    rhem.translate(5, 5, 0);
    rhem.setGLOptions('additive');
    w.addItem(rhem);
    
    
    index = 0;
    def update():
        
        global lhem, rhem, index
        lhem.rotate(1,0,0,1);
        rhem.rotate(1,0,0,1);
        index -= 1;
        
        w.setCameraPosition(distance=750 + index);
        if (750 + index) <= 220:
            w.setCameraPosition(distance=220);
        index = index + 1;
            
        t = QtCore.QTimer()
        t.timeout.connect(update); 
        t.start(50);
        
        #if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()