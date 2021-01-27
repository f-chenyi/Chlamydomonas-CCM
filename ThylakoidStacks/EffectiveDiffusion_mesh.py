# =============================================== #
# ============ Import head packages  ============ #
from fenics import *
import pygmsh as pg
import os
import meshio
import numpy as np
import csv
# =============================================== #
# =============================================== #




def EqualSpacingIndex(m,n):
    return [i*n//m + n//(2*m) for i in range(m)]




def createMeshFull(N_STACK, N_nrrw, d_t, d_h, d_s, d_l, Delta_n, Delta_w, lmesh1, lmesh2, geoname, mshname, xmlname):
    
    '''
    This function outputs the mesh file for certain geometry of the thylakoid stacks,
    which is parameterized by the input variables:
    
        > N_STACKS = number of thylakoid stacks in the simulation domain
        
        > N_nrrw   = number of layers with narrow gaps (see our paper for details)
        
        > d_t      = thickness of the thylakoid membranes
        
        > d_h      = height of the thylakoid lumen
        
        > d_s      = smaller spacing between the thylakoid stacks
        
        > d_l      = larger spacing between the thylakoid stacks
        
        > Delta_n  = width of narrow gaps
        
        > Delta_w  = width of wide gaps
        
    See Fig. S4 of our paper for more details.
    
    
    Other inputs are:
    
        > lmesh1, lmesh2: mesh size of the thylakoid membranes and the space of lumen/stroma
        
        > geoname, mshname, xmlname: directories for output .geo .msh and .xml files
    '''
    
    # y-coordinates of the thylakoid stacks
    ynull  = np.linspace( 2*d_t + d_h + d_l/2 + d_s/2, 1 - (2*d_t + d_h + d_l/2 + d_s/2), N_STACK)
    yminus = ynull - d_t - d_s/2 - d_h/2
    yplus  = ynull + d_t + d_s/2 + d_h/2
    yc     = list(yminus) + list(yplus)
    yc.sort()
    yc_STACK    = np.array(yc)
    
    # x-coordinates of the center of the gaps. Note that we model a geometry where all the gaps are aligned.
    xc_OPEN     = 0.5*np.ones(len(yc_STACK))
    
    
    # gap size
    id_nrrw = EqualSpacingIndex(N_nrrw,N_STACK)
    Delta_STACK = Delta_w  * np.ones(len(yc_STACK))
    Delta_STACK[np.array(id_nrrw)*2]   = Delta_n
    Delta_STACK[np.array(id_nrrw)*2+1] = Delta_n
    
    
    geom            = pg.built_in.Geometry()
    
    # counters of geometric objects
    domain_count = 1
    point_count  = 1
    line_count   = 1
    surf_count   = 1
    
    # draw the bottom line
    p_init   = geom.add_point([1., 0., 0.],lcar = lmesh1)
    p_now    = p_init
    p_next   = geom.add_point([0., 0., 0.],lcar = lmesh1)
    line_now = geom.add_line(p_now,p_next)
    line_sum = [line_now]
    p_now    = p_next
    
    # draw the left part of the domain
    for i in np.arange(len(yc_STACK)):
        p1 = geom.add_point([ 0. ,                                   yc_STACK[i]-d_h/2,   0.0], lcar=lmesh2)
        p2 = geom.add_point([xc_OPEN[i] - Delta_STACK[i]/2 - d_t,  yc_STACK[i]-d_h/2,   0.0], lcar=lmesh2)
        p3 = geom.add_point([xc_OPEN[i] - Delta_STACK[i]/2 - d_t,  yc_STACK[i]+d_h/2,   0.0], lcar=lmesh2)
        p4 = geom.add_point([ 0.,                                    yc_STACK[i]+d_h/2,   0.0], lcar=lmesh2)
        p5 = geom.add_point([ 0.,                                    yc_STACK[i]+d_h/2+d_t, 0.0], lcar=lmesh2)
        p6 = geom.add_point([ xc_OPEN[i] - Delta_STACK[i]/2,         yc_STACK[i]+d_h/2+d_t, 0.0], lcar=lmesh2)
        p7 = geom.add_point([ xc_OPEN[i] - Delta_STACK[i]/2,         yc_STACK[i]-d_h/2-d_t, 0.0], lcar=lmesh2)
        p8 = geom.add_point([ 0.,                                    yc_STACK[i]-d_h/2-d_t, 0.0], lcar=lmesh2)

        line0 = geom.add_line(p_now,p8)
        line1 = geom.add_line(p8,p7)
        line2 = geom.add_line(p7,p6)
        line3 = geom.add_line(p6,p5)
        line4 = geom.add_line(p5,p4)
        line5 = geom.add_line(p4,p3)
        line6 = geom.add_line(p3,p2)
        line7 = geom.add_line(p2,p1)
        line8 = geom.add_line(p1,p8)
        line9 = geom.add_line(p1,p4)

        lineloop1 = geom.add_line_loop([line1,line2,line3,line4,\
                                        line5,line6,line7,line8])
        lineloop2 = geom.add_line_loop([line5,line6,line7,line9])

        line_sum  = line_sum + [line0]
        line_sum  = line_sum + [line1]
        line_sum  = line_sum + [line2]
        line_sum  = line_sum + [line3]

        surf1     = geom.add_plane_surface(lineloop2)
        psurf1    = geom.add_physical(surf1,label=domain_count)
        domain_count += 1
        surf2     = geom.add_plane_surface(lineloop1)
        psurf2    = geom.add_physical(surf2,label=domain_count)
        domain_count += 1

        p_now = p5
    
    # draw the top lines
    p_next = geom.add_point([0., 1., 0.],lcar = lmesh1)
    line_sum = line_sum + [geom.add_line(p_now,p_next)]
    p_now  = p_next

    p_next = geom.add_point([1., 1., 0.],lcar = lmesh1)
    line_sum = line_sum + [geom.add_line(p_now,p_next)]
    p_now  = p_next
    
    # draw the right part
    i = 0
    for i in np.arange(len(yc_STACK))[::-1]:

        p1 = geom.add_point([ 1. ,                                   yc_STACK[i]-d_h/2,   0.0], lcar=lmesh2)
        p2 = geom.add_point([xc_OPEN[i] + Delta_STACK[i]/2 + d_t,  yc_STACK[i]-d_h/2,   0.0], lcar=lmesh2)
        p3 = geom.add_point([xc_OPEN[i] + Delta_STACK[i]/2 + d_t,  yc_STACK[i]+d_h/2,   0.0], lcar=lmesh2)
        p4 = geom.add_point([ 1.,                                    yc_STACK[i]+d_h/2,   0.0], lcar=lmesh2)
        p5 = geom.add_point([ 1.,                                    yc_STACK[i]+d_h/2+d_t, 0.0], lcar=lmesh2)
        p6 = geom.add_point([ xc_OPEN[i] + Delta_STACK[i]/2,         yc_STACK[i]+d_h/2+d_t, 0.0], lcar=lmesh2)
        p7 = geom.add_point([ xc_OPEN[i] + Delta_STACK[i]/2,         yc_STACK[i]-d_h/2-d_t, 0.0], lcar=lmesh2)
        p8 = geom.add_point([ 1.,                                    yc_STACK[i]-d_h/2-d_t, 0.0], lcar=lmesh2)

        line0 = geom.add_line(p_now,p5)
        line1 = geom.add_line(p5,p6)
        line2 = geom.add_line(p6,p7)
        line3 = geom.add_line(p7,p8)
        line4 = geom.add_line(p8,p1)
        line5 = geom.add_line(p1,p2)
        line6 = geom.add_line(p2,p3)
        line7 = geom.add_line(p3,p4)
        line8 = geom.add_line(p4,p5)
        line9 = geom.add_line(p4,p1)

        lineloop1 = geom.add_line_loop([line1,line2,line3,line4,\
                                        line5,line6,line7,line8])
        lineloop2 = geom.add_line_loop([line5,line6,line7,line9])

        line_sum  = line_sum + [line0]
        line_sum  = line_sum + [line1]
        line_sum  = line_sum + [line2]
        line_sum  = line_sum + [line3]

        surf1     = geom.add_plane_surface(lineloop2)
        psurf1    = geom.add_physical(surf1,label=domain_count)
        domain_count += 1
        surf2     = geom.add_plane_surface(lineloop1)
        psurf2    = geom.add_physical(surf2,label=domain_count)
        domain_count += 1

        p_now = p8
    
    # add the stroma domain
    line_sum    = line_sum + [geom.add_line(p_now,p_init)]    # line contour of all the outer boundaries of thyalkoid membranes
    line_stroma = geom.add_line_loop(line_sum)
    surf_stroma = geom.add_plane_surface(line_stroma)
    psurf_stroma= geom.add_physical(surf_stroma,label=domain_count)

    mesh = pg.helpers.generate_mesh(geom,geo_filename=geoname)
    
    # Make sure that Gmsh is installed under the following directory
    os.system('/Applications/Gmsh.app/Contents/MacOS/gmsh %s'%geoname\
               + ' -2 -o %s'%mshname)
    os.system('dolfin-convert %s'%mshname + ' %s'%xmlname)
        
# =============================================== #
# =============================================== #



