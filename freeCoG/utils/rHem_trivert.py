from subcortFs2mlab import *
import os

def get_trivert(subj):

   subjects_dir = '/crtx1/zgreenberg/subjects';
   os.system('cd %s/%s/ascii/' %(subjects_dir, subj));

   # convert all ascii subcortical meshes to matlab vert, tri coords
   subcort = 'aseg_rAcumb.asc';
   nuc = 'rAcumb';
   subcortFs2mlab(subcort,nuc);

   subcort = 'aseg_rAmgd.asc';
   nuc = 'rAmgd';
   subcortFs2mlab(subcort,nuc);

   subcort = 'aseg_rCaud.asc';
   nuc = 'rCaud';
   subcortFs2mlab(subcort,nuc);

   nuc = 'rGP';
   subcort = 'aseg_rGP.asc';
   subcortFs2mlab(subcort,nuc);

   nuc = 'rHipp';
   subcort = 'aseg_rHipp.asc';
   subcortFs2mlab(subcort,nuc);

   subcort = 'aseg_rPut.asc';
   subcort = 'aseg_rPut.asc';
   subcortFs2mlab(subcort,nuc);

   subcort = 'aseg_rThal.asc';
   nuc = 'rThal';
   subcortFs2mlab(subcort,nuc);

   nuc = 'rLatVent';
   subcort = 'aseg_rLatVent.asc';
   subcortFs2mlab(subcort,nuc);

   subcort = 'aseg_rInfLatVent.asc';
   nuc = 'rInfLatVent';
   subcortFs2mlab(subcort,nuc);
