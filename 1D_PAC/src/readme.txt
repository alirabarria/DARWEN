{\rtf1\ansi\ansicpg1252\cocoartf1038\cocoasubrtf360
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww9000\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\ql\qnatural\pardirnatural

\f0\fs24 \cf0 2013_06_05:\
   -CPU time is reduced much using ifort instead of gfortran\
    only available on avakas or venus clusters\
   -various tests done for GJ436b indicate that using the network "olivia"\
    the integration goes until much late times while using the network "room"\
    the integration arrives just to earlier times (and takes much more CPU time).\
    Here, however, the network "olivia" should be bad for temperatures around\
    300-500 K due to several too high rate constants, so that we keep using the\
    network "room" and wait until a new network is created, hopefully before ending\
    the post-doc at Bordeaux.\
\
 #  Need to update read_reac to include: \
     - photoreactions\
     - collisional reactions with fall-off\
     - reverse with rate computed by chemical equilibrium\
# Need to include lj_length for important species such as H, and other atoms\
\
2012_04_13:\
backup before implementing the Newton-Raphson method neglecting molecular diffusion (which allows to use conservation equations for each individual layer and results in a tridiagonal by blocks jacobian matrix)\
\
2012_04_16:\
backup after implementing Newton-Raphson method negleting molecular diffusion.\
Works, but diverges for moderately complicated systems\
}

2018: Robin Baeyens
modification made to include 2D-eddy input grids. 
The 'old' files are still present under the name '_before_I_made_Kzz_2D'.
Only tested for the solid-body rotation case. 
