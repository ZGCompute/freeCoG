# load and plot blender electrodes
import numpy as np

coords = np.loadtxt('hd_grid.txt');
for array in coords:
    X = array[0];
    Y = array[1];
    Z = array[2];
    bpy.ops.mesh.primitive_uv_sphere_add(size=1, view_align=False, enter_editmode=False, location=( X, Y, Z));
    