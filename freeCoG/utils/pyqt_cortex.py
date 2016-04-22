import scipy.io as io
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
import pyqtgraph.opengl as gl
import os
import numpy as np


def pyqt_cort(subj, mesh_dir, subcort_dir, elecs_dir, Grid, depths, Faces, Edge):

    # for white background
    pg.setConfigOption('background', 'w');

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


    #define subjects mesh dir
    os.chdir(mesh_dir);

    #nave to mesh and loat it
    lh_tri = io.loadmat('lh_tri');
    lh_faces = lh_tri.get('tri');
    lh_faces = lh_faces -1;

    lh_vert = io.loadmat('lh_vert')
    lh_verts = lh_vert.get('vert')

    #nave to mesh and loat it
    rh_tri = io.loadmat('rh_tri');
    rh_faces = rh_tri.get('tri');
    rh_faces = rh_faces -1;

    rh_vert = io.loadmat('rh_vert')
    rh_verts = rh_vert.get('vert')

    # go to subcort dir
    os.chdir(subcort_dir);

    # set up hippocampus mesh
    lHipp_tri = io.loadmat('subcort_lHipp_tri');
    lHipp_faces = lHipp_tri.get('tri');
    lHipp_faces = lHipp_faces -1;

    lHipp_vert = io.loadmat('subcort_lHipp_vert')
    lHipp_verts = lHipp_vert.get('vert')

    # set up amgd mesh
    lAmgd_tri = io.loadmat('subcort_lAmgd_tri');
    lAmgd_faces = lAmgd_tri.get('tri');
    lAmgd_faces = lAmgd_faces -1;

    lAmgd_vert = io.loadmat('subcort_lAmgd_vert')
    lAmgd_verts = lAmgd_vert.get('vert')

    # set up GP mesh
    lGP_tri = io.loadmat('subcort_lGP_tri');
    lGP_faces = lGP_tri.get('tri');
    lGP_faces = lGP_faces -1;

    lGP_vert = io.loadmat('subcort_lGP_vert')
    lGP_verts = lAmgd_vert.get('vert')

    # set up Caud mesh
    lCaud_tri = io.loadmat('subcort_lCaud_tri');
    lCaud_faces = lCaud_tri.get('tri');
    lCaud_faces = lCaud_faces -1;

    lCaud_vert = io.loadmat('subcort_lCaud_vert')
    lCaud_verts = lCaud_vert.get('vert')

    # set up Put mesh
    lPut_tri = io.loadmat('subcort_lPut_tri');
    lPut_faces = lPut_tri.get('tri');
    lPut_faces = lPut_faces -1;

    lPut_vert = io.loadmat('subcort_lPut_vert')
    lPut_verts = lPut_vert.get('vert')

    # set up Thal mesh
    lThal_tri = io.loadmat('subcort_lThal_tri');
    lThal_faces = lThal_tri.get('tri');
    lThal_faces = lThal_faces -1;

    lThal_vert = io.loadmat('subcort_lThal_vert')
    lThal_verts = lThal_vert.get('vert')

    # set up hippocampus mesh
    rHipp_tri = io.loadmat('subcort_rHipp_tri');
    rHipp_faces = rHipp_tri.get('tri');
    rHipp_faces = rHipp_faces -1;

    rHipp_vert = io.loadmat('subcort_rHipp_vert')
    rHipp_verts = rHipp_vert.get('vert')

    # set up amgd mesh
    rAmgd_tri = io.loadmat('subcort_rAmgd_tri');
    rAmgd_faces = rAmgd_tri.get('tri');
    rAmgd_faces = rAmgd_faces -1;

    rAmgd_vert = io.loadmat('subcort_rAmgd_vert')
    rAmgd_verts = rAmgd_vert.get('vert')

    # set up GP mesh
    rGP_tri = io.loadmat('subcort_rGP_tri');
    rGP_faces = rGP_tri.get('tri');
    rGP_faces = rGP_faces -1;

    rGP_vert = io.loadmat('subcort_rGP_vert')
    rGP_verts = rGP_vert.get('vert')

    # set up Caud mesh
    rCaud_tri = io.loadmat('subcort_rCaud_tri');
    rCaud_faces = rCaud_tri.get('tri');
    rCaud_faces = rCaud_faces -1;

    rCaud_vert = io.loadmat('subcort_rCaud_vert')
    rCaud_verts = rCaud_vert.get('vert')

    # set up Put mesh
    rPut_tri = io.loadmat('subcort_rPut_tri');
    rPut_faces = rPut_tri.get('tri');
    rPut_faces = rPut_faces -1;

    rPut_vert = io.loadmat('subcort_rPut_vert')
    rPut_verts = rPut_vert.get('vert')

    # set up Thal mesh
    rThal_tri = io.loadmat('subcort_rThal_tri');
    rThal_faces = rThal_tri.get('tri');
    rThal_faces = rThal_faces -1;

    rThal_vert = io.loadmat('subcort_rThal_vert')
    rThal_verts = rThal_vert.get('vert')

    # go to fiber mesh dir
    #subj_fiberMesh_dir = '/Users/Zach/Desktop/Data/MR_data/EC69/DTI_HARDI_55_dirs/dti55trilin/fiberMeshes/lARC';
    #os.chdir(subj_fiberMesh_dir);

    # set up fiber mesh
    #lARC_tri = io.loadmat('lARC_fiber_1');
    #lARC_faces = lARC_tri.get('fiberMesh');
    #lARC_faces = lARC_faces -1;


    # launch opengl app
    app = QtGui.QApplication([])
    w = gl.GLViewWidget()
    w.show()
    w.setWindowTitle(('Patient Pial Surface and Depths for: ' + subj))
    w.setCameraPosition(distance=750);

    # set grid on
    g = gl.GLGridItem()
    g.scale(100,100,50)
    #w.addItem(g)


    # add electrodes as spheres
    os.chdir(elecs_dir);
    grid = io.loadmat('depth_probes.mat');
    grid = grid.get('elecmatrix');
    m1 = gl.GLScatterPlotItem(pos=grid, size = 15., color=(1,0,0,1), pxMode=True);
    m1.translate(5, 5, 0);
    m1.setDepthValue(1);
    m1.setGLOptions('additive');

    #grid
    grid = io.loadmat('hd_grid.mat');
    grid = grid.get('elecmatrix');
    m4 = gl.GLScatterPlotItem(pos=grid, size = 15., color=(1,0,0,1), pxMode=True);
    m4.translate(5, 5, 0);
    m4.setDepthValue(1);
    m4.setGLOptions('translucent');

    phase = 0.;

    #load in ecog data
    #subj_data_dir = '/Users/Zach/Desktop/Data/MR_data/EC69/ECoG';
    #os.chdir(subj_data_dir);
    #data = io.loadmat('hg_listen_avg.mat');
    #avg_trial = data.get('avg_trial');
    #d2 = (grid**2).sum(axis=1)**0.5;
    
    # add depths and grid to plot
    w.addItem(m1);
    if Grid == True:
       w.addItem(m4);

    #z = avg_trial;
    #p4 = gl.GLSurfacePlotItem(x=grid[:,0], y = grid[:,1], shader='heightColor', computeNormals=False, smooth=False); #
    #p4.shader()['colorMap'] = np.array([0.2, 2, 0.5, 0.2, 1, 1, 0.2, 0, 2]); #avg_trial[:,0];
    #p4.translate(10, 10, 0)
    #w.addItem(p4)

    # An alternative to the custom shader is to invert the mesh normals:
    #meshdata = gl.MeshData(vertexes=faces);
    #norms = meshdata.faceNormals();
    #cd norms *= - 1;

    # set up lh mesh
    if Faces == True:
        edgeColor = (0,0,0,1);
    else:
        edgeColor = (1,1,1,0.2);
    colors1 = np.random.random(size=(lh_verts.shape[0], 3, 4));
    m2 = gl.GLMeshItem(vertexes = lh_verts, faces=lh_faces, drawFaces=Faces, drawEdges=Edge, faceColor=(1,0,0,0),edgeColor=edgeColor, smooth=True, shader='my_shader');
    m2.translate(5, 5, 0)
    m2.setGLOptions('translucent');#m2.setGLOptions('additive')
    m2.setDepthValue(0.5);
    w.addItem(m2);

    # add rh mesh
    colors2 = np.random.random(size=(rh_verts.shape[0], 3, 4));
    m3 = gl.GLMeshItem(vertexes = rh_verts, faces=rh_faces, drawFaces=Faces, drawEdges=Edge, faceColor=(1,0,0,0),edgeColor=edgeColor, smooth=True, shader='my_shader');
    m3.translate(5, 5, 0)
    m3.setGLOptions('translucent');#m2.setGLOptions('additive')
    m3.setDepthValue(0.5);
    w.addItem(m3);

    # add the hippocampus
    lHipp = gl.GLMeshItem(vertexes = lHipp_verts, faces=lHipp_faces, drawFaces=False, drawEdges=True, faceColor=(1,0,1,0.4),edgeColor=(1,0,0.5,0.2), smooth=True, shader='my_shader');
    lHipp.translate(5, 5, 0)
    lHipp.setGLOptions('translucent');#m2.setGLOptions('additive')
    w.addItem(lHipp);

    # add the Amygdala
    lAmgd = gl.GLMeshItem(vertexes = lAmgd_verts, faces=lAmgd_faces, drawFaces=False, drawEdges=True, faceColor=(0,0.4,0.4,0.4),edgeColor=(1,0.4,0,0.2), smooth=True, shader='my_shader');
    lAmgd.translate(5, 5, 0)
    lAmgd.setGLOptions('translucent');#m2.setGLOptions('additive')
    w.addItem(lAmgd);

    # add the Globus Pallidus
    lGP = gl.GLMeshItem(vertexes = lGP_verts, faces=lGP_faces, drawFaces=False, drawEdges=True, faceColor=(0,0.4,0.4,0.4),edgeColor=(1,0.8,0.2,0.2), smooth=True, shader='my_shader');
    lGP.translate(5, 5, 0)
    lGP.setGLOptions('translucent');#m2.setGLOptions('additive')
    w.addItem(lGP);

    # add the Caudate
    lCaud = gl.GLMeshItem(vertexes = lCaud_verts, faces=lCaud_faces, drawFaces=False, drawEdges=True, faceColor=(0,0.4,0.4,0.4),edgeColor=(1,0.2,0.2,0.2), smooth=True, shader='my_shader');
    lCaud.translate(5, 5, 0)
    lCaud.setGLOptions('translucent');#m2.setGLOptions('additive')
    w.addItem(lCaud);

    # add the Putamen
    lPut = gl.GLMeshItem(vertexes = lPut_verts, faces=lPut_faces, drawFaces=False, drawEdges=True, faceColor=(0,0.4,0.4,0.4),edgeColor=(0.2,0.3,0.5,0.2), smooth=True, shader='my_shader');
    lPut.translate(5, 5, 0)
    lPut.setGLOptions('translucent');#m2.setGLOptions('additive')
    w.addItem(lPut);

    # add the Thalimus
    lThal = gl.GLMeshItem(vertexes = lThal_verts, faces=lThal_faces, drawFaces=False, drawEdges=True, faceColor=(0,0.4,0.4,0.4),edgeColor=(0.1,0.7,0.4,0.2), smooth=True, shader='my_shader');
    lThal.translate(5, 5, 0)
    lThal.setGLOptions('translucent');#m2.setGLOptions('additive')
    w.addItem(lThal);

    # add the Thalimus
    #lARC = gl.GLMeshItem(faces=lARC_faces, drawFaces=False, drawEdges=True, faceColor=(1,0.4,0.4,0.4),edgeColor=(1,0.7,0.4,0.2), smooth=True, shader='my_shader');
    #lARC.translate(5, 5, 0)
    #lARC.setGLOptions('translucent');#m2.setGLOptions('additive')
    #w.addItem(lARC);
    
    # add the hippocampus
    rHipp = gl.GLMeshItem(vertexes = rHipp_verts, faces=rHipp_faces, drawFaces=False, drawEdges=True, faceColor=(1,0,1,0.4),edgeColor=(1,0,0.5,0.2), smooth=True, shader='my_shader');
    rHipp.translate(5, 5, 0)
    rHipp.setGLOptions('translucent');#m2.setGLOptions('additive')
    w.addItem(rHipp);

    # add the Amygdala
    rAmgd = gl.GLMeshItem(vertexes = rAmgd_verts, faces=rAmgd_faces, drawFaces=False, drawEdges=True, faceColor=(0,0.4,0.4,0.4),edgeColor=(1,0.4,0,0.2), smooth=True, shader='my_shader');
    rAmgd.translate(5, 5, 0)
    rAmgd.setGLOptions('translucent');#m2.setGLOptions('additive')
    w.addItem(rAmgd);

    # add the Globus Pallidus
    rGP = gl.GLMeshItem(vertexes = rGP_verts, faces=rGP_faces, drawFaces=False, drawEdges=True, faceColor=(0,0.4,0.4,0.4),edgeColor=(1,0.8,0.2,0.2), smooth=True, shader='my_shader');
    rGP.translate(5, 5, 0)
    rGP.setGLOptions('translucent');#m2.setGLOptions('additive')
    w.addItem(rGP);

    # add the Caudate
    rCaud = gl.GLMeshItem(vertexes = rCaud_verts, faces=rCaud_faces, drawFaces=False, drawEdges=True, faceColor=(0,0.4,0.4,0.4),edgeColor=(1,0.2,0.2,0.2), smooth=True, shader='my_shader');
    rCaud.translate(5, 5, 0)
    rCaud.setGLOptions('translucent');#m2.setGLOptions('additive')
    w.addItem(rCaud);

    # add the Putamen
    rPut = gl.GLMeshItem(vertexes = rPut_verts, faces=rPut_faces, drawFaces=False, drawEdges=True, faceColor=(0,0.4,0.4,0.4),edgeColor=(0.2,0.3,0.5,0.2), smooth=True, shader='my_shader');
    rPut.translate(5, 5, 0)
    rPut.setGLOptions('translucent');#m2.setGLOptions('additive')
    w.addItem(rPut);

    # add the Thalimus
    rThal = gl.GLMeshItem(vertexes = rThal_verts, faces=rThal_faces, drawFaces=False, drawEdges=True, faceColor=(0,0.4,0.4,0.4),edgeColor=(0.1,0.7,0.4,0.2), smooth=True, shader='my_shader');
    rThal.translate(5, 5, 0)
    rThal.setGLOptions('translucent');#m2.setGLOptions('additive')
    w.addItem(rThal);

    # add the Thalimus
    #lARC = gl.GLMeshItem(faces=lARC_faces, drawFaces=False, drawEdges=True, faceColor=(1,0.4,0.4,0.4),edgeColor=(1,0.7,0.4,0.2), smooth=True, shader='my_shader');
    #lARC.translate(5, 5, 0)
    #lARC.setGLOptions('translucent');#m2.setGLOptions('additive')
    #w.addItem(lARC);


    # update over 50 iterations the color value plotted at each electrode
    index = 0;
    def update():
        ## update volume  colors
        global m1, m4, d2, index
        m1.rotate(1,0,0,1);
        m2.rotate(1,0,0,1);
        m3.rotate(1,0,0,1);
        if Grid == True:
            m4.rotate(1,0,0,1);
        index -= 1;
        #p4.setData(z=z[:,index%z.shape[1]]);
        #g.rotate(1,0,0,1);
        lHipp.rotate(1,0,0,1);
        lAmgd.rotate(1,0,0,1);
        lGP.rotate(1,0,0,1);
        lCaud.rotate(1,0,0,1);
        lPut.rotate(1,0,0,1);
        lThal.rotate(1,0,0,1);
        rHipp.rotate(1,0,0,1);
        rAmgd.rotate(1,0,0,1);
        rGP.rotate(1,0,0,1);
        rCaud.rotate(1,0,0,1);
        rPut.rotate(1,0,0,1);
        rThal.rotate(1,0,0,1);
        #lARC.rotate(1,0,0,1);

        w.setCameraPosition(distance=750 + index);
        if (750 + index) <= 220:
            w.setCameraPosition(distance=220);
        
        index = index + 1;

            
        
        
    t = QtCore.QTimer()
    t.timeout.connect(update)
    t.start(50)

    # launch mesh plot
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
