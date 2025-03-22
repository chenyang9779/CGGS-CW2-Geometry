#! /bin/bash

cd /home/chenyang/repos/CGGS-CW2-Geometry

mv ~/Downloads/CGGS_cw2_Geometry.pdf s2693586/s2693586.pdf


cp code/include/compute_laplacian.h code/include/compute_mean_curvature_normal.h code/include/compute_angle_defect.h s2693586/section1/

cp code/include/compute_mean_curvature_normal.h code/include/compute_angle_defect.h code/include/compute_laplacian.h code/include/mean_curvature_flow.h s2693586/section2/

cp code/include/tutte_parameterization.h code/include/compute_laplacian.h s2693586/section3/

cp -r code/Section4 s2693586