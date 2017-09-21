#seg_3d.py

# function to manually clean outputs of 3D 
# electrode detection

import time
import scipy.io
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def seg_3D(CT_dir,fname, out_fname):

    # load in 3D output
    coords = scipy.io.loadmat(CT_dir + '/' + fname).get('elecmatrix');
    x = coords[:,0];
    y = coords[:,1];
    z = coords[:,2];

    # plot dirty elec results
    s = 150;
    fig = plt.figure(facecolor="white");
    ax = fig.add_subplot(111,projection='3d');
    

    # create point and click function for selecting false alarm points
    seg_elecs = [];
    def onpick(event):

        # get event indices and x y z coord for click
        ind = event.ind[0];
        x_e, y_e, z_e = event.artist._offsets3d;
        print x_e[ind], y_e[ind], z_e[ind];
        seg_coord = [x_e[ind], y_e[ind], z_e[ind]];
        seg_elecs.append(seg_coord);

        # refresh plot with red dots picked and blue dots clean
        #fig.clf();
        seg_arr = np.array(seg_elecs);
        x_f = seg_arr[:,0];
        y_f = seg_arr[:,1];
        z_f = seg_arr[:,2];
        plt.cla();
        Axes3D.scatter3D(ax,x_f,y_f,z_f,s=150,c='r', picker=5);

        # get array of only clean elecs to re-plot as blue dots
        clean_list = list(coords);
        clean_list = [list(i) for i in clean_list];
        for coordin in clean_list:
                if list(coordin) in seg_elecs:
                    clean_list.pop(clean_list.index(coordin));
        clean_coordis = np.array(clean_list);
        x_c = clean_coordis[:,0];
        y_c = clean_coordis[:,1];
        z_c = clean_coordis[:,2];
        Axes3D.scatter3D(ax,x_c, y_c, z_c, s=150,c='b', picker=5);
        time.sleep(0.5);
        plt.draw();
    
    Axes3D.scatter3D(ax,x,y,z,s=s, picker=5);
    #title_font = {'fontname':'serif', 'size':'24', 'color':'black', 'weight':'bold',                            
    #'verticalalignment':'bottom'};                                                                        
    plt.title('3D Output Coordinates: ');
    fig.canvas.mpl_connect('pick_event', onpick);
    plt.ion();
    plt.show(fig);

    # wait for user to press enter to continue
    while True:
        
        fig.canvas.mpl_connect('pick_event', onpick);
        plt.draw();
        i = raw_input("Enter text or Enter to quit: ");
        if i == '':
            
            # remove all false coords from elec centers list
            fig.clf();
            coords_list = list(coords);
            coords_list = [list(i) for i in coords_list];
            for coordi in coords_list:
                if list(coordi) in seg_elecs:
                    coords_list.pop(coords_list.index(coordi));
            clean_coords = np.array(seg_elecs)

            # plot result
            X = clean_coords[:,0];
            Y = clean_coords[:,1];
            Z = clean_coords[:,2];


            s = 150;
            fig = plt.figure(facecolor="white");
            ax = fig.add_subplot(111,projection='3d');
            Axes3D.scatter3D(ax,X,Y,Z,s=s);
            plt.show(fig);
            
            # save the cleaned elecs file
            new_fname = CT_dir + '/' + out_fname;            
            scipy.io.savemat(new_fname, {'elecmatrix':clean_coords});
            
            return clean_coords;
            break;

        
   
    
        
    
