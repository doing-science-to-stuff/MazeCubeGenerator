/*
 * maze-export.c
 * MazeCubeGen: maze cube generator
 *
 * Copyright (c) 2020-2024 Bryan Franklin. All rights reserved.
 */
#include <errno.h>
#include <math.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include "maze.h"
#include "maze-export.h"
#include "trig-list.h"

/* circular marker */
static void maze_add_marker1(trig_list_t *list, maze_t *maze, int face, position_t pos, double radius, double scale) {
    trig_list_t marker;
    trig_list_init(&marker);

    int d1 = maze->faces[face].d1;
    int d2 = maze->faces[face].d2;

    /* get row and column for marker */
    int c = pos[d2];
    int r = pos[d1];

    const int numSegsI = 64;
    const int numSegsJ = 8;
    for(int i=0; i<numSegsI; ++i) {
        /* compute points along circle */
        double thetaI1 =  2.0 * M_PI * i / numSegsI;
        double x1 = cos(thetaI1)*(1-radius);
        double y1 = sin(thetaI1)*(1-radius);
        double thetaI2 =  2.0 * M_PI * (i+1) / numSegsI;
        double x2 = cos(thetaI2)*(1-radius);
        double y2 = sin(thetaI2)*(1-radius);

        /* position circle points onto face */
        int col1 = y1+c+0.5;
        int col2 = y2+c+0.5;
        int row1 = x1+r+0.5;
        int row2 = x2+r+0.5;

        /* check if cell under each point is set */
        int cell1 = face_get_cell(&maze->faces[face], row1, col1);
        int cell2 = face_get_cell(&maze->faces[face], row2, col2);

        double xc0 = 0.0, yc0 = 0.0, zc0 = 0.0;
        for(int j=0; j<=numSegsJ; ++j) {
            double thetaJ1 =  M_PI * j / numSegsJ;
            double thetaJ2 =  M_PI * (j+1) / numSegsJ;

            /* deal with back-side */
            if( j == numSegsJ ) { thetaJ2 = 0.0; }

            /* compute raw torus coordinates */
            double x11 = cos(thetaI1) * (1+radius*cos(thetaJ1)) + r;
            double y11 = sin(thetaI1) * (1+radius*cos(thetaJ1)) + c;
            double z11 = radius*sin(thetaJ1) + 0.5;

            double x12 = cos(thetaI1) * (1+radius*cos(thetaJ2)) + r;
            double y12 = sin(thetaI1) * (1+radius*cos(thetaJ2)) + c;
            double z12 = radius*sin(thetaJ2) + 0.5;

            double x21 = cos(thetaI2) * (1+radius*cos(thetaJ1)) + r;
            double y21 = sin(thetaI2) * (1+radius*cos(thetaJ1)) + c;
            double z21 = radius*sin(thetaJ1)+0.5;

            double x22 = cos(thetaI2) * (1+radius*cos(thetaJ2)) + r;
            double y22 = sin(thetaI2) * (1+radius*cos(thetaJ2)) + c;
            double z22 = radius*sin(thetaJ2) + 0.5;

            /* compute normals */
            double x1c = cos(thetaI1) + r;
            double y1c = sin(thetaI1) + c;
            double z1c = 0.5;
            double x2c = cos(thetaI2) + r;
            double y2c = sin(thetaI2) + c;
            double z2c = 0.5;

            double nx11 = x11 - x1c;
            double ny11 = y11 - y1c;
            double nz11 = z11 - z1c;

            double nx12 = x12 - x1c;
            double ny12 = y12 - y1c;
            double nz12 = z12 - z1c;

            double nx21 = x21 - x2c;
            double ny21 = y21 - y2c;
            double nz21 = z21 - z2c;

            double nx22 = x22 - x2c;
            double ny22 = y22 - y2c;
            double nz22 = z22 - z2c;

            /* interpolate crossing point when transitioning to/from
             * an open cell */
            double xc1 = 0.0, yc1 = 0.0, zc1 = 0.0;
            double xc2 = 0.0, yc2 = 0.0, zc2 = 0.0;
            double nxc1 = 0.0, nyc1 = 0.0, nzc1 = 0.0;
            double nxc2 = 0.0, nyc2 = 0.0, nzc2 = 0.0;
            if( row1 != row2 ) {
                /* compute crossing points */
                double t1 = (round(x21*2)/2-x11)/(x21-x11);
                xc1 = x11 + t1*(x21-x11);
                yc1 = y11 + t1*(y21-y11);
                zc1 = z11 + t1*(z21-z11);
                nxc1 = nx11 + t1*(nx21-nx11);
                nyc1 = ny11 + t1*(ny21-ny11);
                nzc1 = nz11 + t1*(nz21-nz11);
                double t2 = (round(x22*2)/2-x12)/(x22-x12);
                xc2 = x12 + t2*(x22-x12);
                yc2 = y12 + t2*(y22-y12);
                zc2 = z12 + t2*(z22-z12);
                nxc2 = nx12 + t2*(nx22-nx12);
                nyc2 = ny12 + t2*(ny22-ny12);
                nzc2 = nz12 + t2*(nz22-nz12);
            }
            if( col1 != col2 ) {
                /* compute crossing points */
                double t1 = (round(y21*2)/2-y11)/(y21-y11);
                xc1 = x11 + t1*(x21-x11);
                yc1 = y11 + t1*(y21-y11);
                zc1 = z11 + t1*(z21-z11);
                nxc1 = nx11 + t1*(nx21-nx11);
                nyc1 = ny11 + t1*(ny21-ny11);
                nzc1 = nz11 + t1*(nz21-nz11);
                double t2 = (round(y22*2)/2-y12)/(y22-y12);
                xc2 = x12 + t2*(x22-x12);
                yc2 = y12 + t2*(y22-y12);
                zc2 = z12 + t2*(z22-z12);
                nxc2 = nx12 + t2*(nx22-nx12);
                nyc2 = ny12 + t2*(ny22-ny12);
                nzc2 = nz12 + t2*(nz22-nz12);
            }

            /* replace appropriate values */
            if( cell1!=0 && cell2==0 ) {
                x21 = xc1; y21 = yc1; z21 = zc1;
                x22 = xc2; y22 = yc2; z22 = zc2;
                nx21 = nxc1; ny21 = nyc1; nz21 = nzc1;
                nx22 = nxc2; ny22 = nyc2; nz22 = nzc2;
            }
            else if( cell1==0 && cell2!=0 ) {
                x11 = xc1; y11 = yc1; z11 = zc1;
                x12 = xc2; y12 = yc2; z12 = zc2;
                nx11 = nxc1; ny11 = nyc1; nz11 = nzc1;
                nx12 = nxc2; ny12 = nyc2; nz12 = nzc2;
            }

            /* record first point along surface for end caps */
            if( cell1 != cell2 && j == 0 ) {
                xc0 = xc1; yc0 = yc1; zc0 = zc1;
            }

            /* output triangles */
            if( cell1 || cell2 ) {
                /* curved surface of marker */
                trig_t t1, t2;
                trig_fill(&t1,
                        x11, y11, z11,
                        x22, y22, z22,
                        x12, y12, z12);
                trig_fill(&t2,
                        x11, y11, z11,
                        x21, y21, z21,
                        x22, y22, z22);
                trig_set_normals(&t1,
                                 nx11, ny11, nz11,
                                 nx22, ny22, nz22,
                                 nx12, ny12, nz12);
                trig_set_normals(&t2,
                                 nx11, ny11, nz11,
                                 nx21, ny21, nz21,
                                 nx22, ny22, nz22);
                trig_list_append(&marker, &t1);
                trig_list_append(&marker, &t2);
            }
            if( cell1==0 && cell2!=0 && j>0 ) {
                /* cap one end */
                trig_list_add(&marker, xc0, yc0, zc0,
                                    xc1, yc1, zc1,
                                    xc2, yc2, zc2);
                /* end-caps are flat, so default normals are fine */
            }
            if( cell1!=0 && cell2==0 && j>0 ) {
                /* cap other end */
                trig_list_add(&marker, xc0, yc0, zc0,
                                    xc2, yc2, zc2,
                                    xc1, yc1, zc1);
                /* end-caps are flat, so default normals are fine */
            }
        }
    }
    
    /* assign marker to a group */
    trig_list_set_groupid(&marker, FACE_MARKER_1);
    
    /* append marker into passed-in list */
    trig_list_concatenate(list, &marker);
    trig_list_free(&marker);
}


/* square marker */
static void maze_add_corner(trig_list_t *list, maze_t *maze, int r, int c, int dr, int dc, int face, double radius, double scale, int rCap, int cCap) {

    /* local triangle list */
    trig_list_t corner;
    trig_list_init(&corner);

    double x10 = 0.0, y10 = 0.0, z10 = 0.0;
    double x30 = 0.0, y30 = 0.0, z30 = 0.0;

    int numSegs = 32;
    for(int i=0; i<=numSegs; ++i) {
        double thetaI1 = M_PI * i / numSegs;
        double thetaI2 = M_PI * (i+1) / numSegs;

        /* deal with back-side */
        if( i == numSegs ) { thetaI2 = 0.0; }

        /* compute raw segment coordinates */
        double x11 = radius*cos(thetaI1)+1.0;
        double y11 = 0.5;
        double z11 = radius*sin(thetaI1)+0.5;

        double x12 = radius*cos(thetaI2)+1.0;
        double y12 = 0.5;
        double z12 = radius*sin(thetaI2)+0.5;

        double x21 = radius*cos(thetaI1)+1.0;
        double y21 = radius*cos(thetaI1)+1.0;
        double z21 = radius*sin(thetaI1)+0.5;

        double x22 = radius*cos(thetaI2)+1.0;
        double y22 = radius*cos(thetaI2)+1.0;
        double z22 = radius*sin(thetaI2)+0.5;

        double x31 = 0.5;
        double y31 = radius*cos(thetaI1)+1.0;
        double z31 = radius*sin(thetaI1)+0.5;

        double x32 = 0.5;
        double y32 = radius*cos(thetaI2)+1.0;
        double z32 = radius*sin(thetaI2)+0.5;

        /* record first point on surface for end caps */
        if( i == 0 ) {
            x10 = x11;
            y10 = y11;
            z10 = z11;
            x30 = x31;
            y30 = y31;
            z30 = z31;
        }

        /* compute normals */
        double x1c = 1.0;
        double y1c = 0.5;
        double z1c = 0.5;

        double x2c = 1.0;
        double y2c = 1.0;
        double z2c = 0.5;

        double x3c = 0.5;
        double y3c = 1.0;
        double z3c = 0.5;

        double nx11 = x11 - x1c;
        double ny11 = y11 - y1c;
        double nz11 = z11 - z1c;
        
        double nx12 = x12 - x1c;
        double ny12 = y12 - y1c;
        double nz12 = z12 - z1c;
            
        double nx21 = x21 - x2c;
        double ny21 = y21 - y2c;
        double nz21 = z21 - z2c;
            
        double nx22 = x22 - x2c;
        double ny22 = y22 - y2c;
        double nz22 = z22 - z2c;
        
        double nx31 = x31 - x3c;
        double ny31 = y31 - y3c;
        double nz31 = z31 - z3c;
            
        double nx32 = x32 - x3c;
        double ny32 = y32 - y3c;
        double nz32 = z32 - z3c;
        
        /* curved surface of marker */
        trig_t t1, t2, t3, t4;
        trig_fill(&t1,
                x11, y11, z11,
                x22, y22, z22,
                x12, y12, z12);
        trig_set_normals(&t1,
                         nx11, ny11, nz11,
                         nx22, ny22, nz22,
                         nx12, ny12, nz12);
        trig_fill(&t2,
                x11, y11, z11,
                x21, y21, z21,
                x22, y22, z22);
        trig_set_normals(&t2,
                         nx11, ny11, nz11,
                         nx21, ny21, nz21,
                         nx22, ny22, nz22);

        trig_fill(&t3,
                x31, y31, z31,
                x32, y32, z32,
                x22, y22, z22);
        trig_set_normals(&t3,
                         nx31, ny31, nz31,
                         nx32, ny32, nz32,
                         nx22, ny22, nz22);
        trig_fill(&t4,
                x31, y31, z31,
                x22, y22, z22,
                x21, y21, z21);
        trig_set_normals(&t4,
                         nx31, ny31, nz31,
                         nx22, ny22, nz22,
                         nx21, ny21, nz21);

        trig_list_append(&corner, &t1);
        trig_list_append(&corner, &t2);
        trig_list_append(&corner, &t3);
        trig_list_append(&corner, &t4);

        /* add end caps */
        if( rCap==0 ) {
            trig_list_add(&corner, x11, y11, z11,
                    x12, y12, z12,
                    x10, y10, z10);
            /* end-caps are flat, so default normals are fine */
        }

        if( cCap==0 ) {
            trig_list_add(&corner, x31, y31, z31,
                    x30, y30, z30,
                    x32, y32, z32);
            /* end-caps are flat, so default normals are fine */
        }
    }

    trig_list_scale(&corner, dr, dc, 1.0);
    trig_list_move(&corner, r, c, 0.0);

    trig_list_concatenate(list, &corner);
    trig_list_free(&corner);
}


static void maze_add_edge(trig_list_t *list, maze_t *maze, int r, int c, int face, double radius, double scale, int rotated, double rOffset, double cOffset) {

    /* local triangle list */
    trig_list_t edge;
    trig_list_init(&edge);

    int numSegs = 32;
    for(int i=0; i<=numSegs; ++i) {
        double thetaI1 = M_PI * i / numSegs;
        double thetaI2 = M_PI * (i+1) / numSegs;

        /* deal with back-side */
        if( i == numSegs ) { thetaI2 = 0.0; }

        /* compute raw segment coordinates */
        double x11 = radius*cos(thetaI1);
        double y11 = 0.5;
        double z11 = radius*sin(thetaI1)+0.5;

        double x12 = radius*cos(thetaI2);
        double y12 = 0.5;
        double z12 = radius*sin(thetaI2)+0.5;

        double x21 = radius*cos(thetaI1);
        double y21 = -0.5;
        double z21 = radius*sin(thetaI1)+0.5;

        double x22 = radius*cos(thetaI2);
        double y22 = -0.5;
        double z22 = radius*sin(thetaI2)+0.5;

        /* compute normals */
        double x1c = 0.0;
        double y1c = 0.5;
        double z1c = 0.5;
        
        double x2c = 0.0;
        double y2c = -0.5;
        double z2c = 0.5;
        
        double nx11 = x11 - x1c;
        double ny11 = y11 - y1c;
        double nz11 = z11 - z1c;

        double nx12 = x12 - x1c;
        double ny12 = y12 - y1c;
        double nz12 = z12 - z1c;

        double nx21 = x21 - x2c;
        double ny21 = y21 - y2c;
        double nz21 = z21 - z2c;

        double nx22 = x22 - x2c;
        double ny22 = y22 - y2c;
        double nz22 = z22 - z2c;
        
        /* outer shell of marker */
        trig_t t1, t2;
        trig_fill(&t1,
                x11, y11, z11,
                x12, y12, z12,
                x22, y22, z22);
        trig_fill(&t2,
                x11, y11, z11,
                x22, y22, z22,
                x21, y21, z21);
        trig_set_normals(&t1,
                         nx11, ny11, nz11,
                         nx12, ny12, nz12,
                         nx22, ny22, nz22);
        trig_set_normals(&t2,
                         nx11, ny11, nz11,
                         nx22, ny22, nz22,
                         nx21, ny21, nz21);
        trig_list_append(&edge, &t1);
        trig_list_append(&edge, &t2);
    }

    if( rotated != 0 )
        trig_list_rotate_axial_around(&edge, 2, M_PI/2.0, 0.0, 0.0, 0.0);
    trig_list_move(&edge, r+rOffset, c+cOffset, 0.0);

    trig_list_concatenate(list, &edge);
    trig_list_free(&edge);
}


static void maze_add_marker2(trig_list_t *list, maze_t *maze, int face, position_t pos, double radius, double scale) {
    trig_list_t marker;
    trig_list_init(&marker);

    int d1 = maze->faces[face].d1;
    int d2 = maze->faces[face].d2;

    /* get row and column for marker */
    int c = pos[d2];
    int r = pos[d1];

    /* check which straight segments are needed */
    int lSide = face_get_cell(&maze->faces[face], r-1, c);
    int rSide = face_get_cell(&maze->faces[face], r+1, c);
    int tSide = face_get_cell(&maze->faces[face], r, c-1);
    int bSide = face_get_cell(&maze->faces[face], r, c+1);

    #if 1
    /* add corners */
    maze_add_corner(&marker, maze, r, c, -1, -1, face, radius, scale, lSide, tSide);
    maze_add_corner(&marker, maze, r, c, -1, 1, face, radius, scale, lSide, bSide);
    maze_add_corner(&marker, maze, r, c, 1, -1, face, radius, scale, rSide, tSide);
    maze_add_corner(&marker, maze, r, c, 1, 1, face, radius, scale, rSide, bSide);
    #endif // 1

    #if 1
    /* add straight segments */
    if( lSide != 0 )
        maze_add_edge(&marker, maze, r, c, face, radius, scale, 0, -1, 0);
    if( rSide != 0 )
        maze_add_edge(&marker, maze, r, c, face, radius, scale, 0, 1, 0);
    if( tSide != 0 )
        maze_add_edge(&marker, maze, r, c, face, radius, scale, 1, 0, -1);
    if( bSide != 0 )
        maze_add_edge(&marker, maze, r, c, face, radius, scale, 1, 0, 1);
    #endif // 1
    
    /* assign marker to a group */
    trig_list_set_groupid(&marker, FACE_MARKER_2);
    
    /* append marker into passed-in list */
    trig_list_concatenate(list, &marker);
    trig_list_free(&marker);
}

static void maze_add_cube(trig_list_t *list, double x, double y, double z, char face_mask, double scaleX, double scaleY, double scaleZ) {
    double dx = scaleX/2.0;
    double dy = scaleY/2.0;
    double dz = scaleZ/2.0;

    #if 0
    /* FOR DEBUGGING PURPOSES ONLY */
    //face_mask = 0x3f;
    face_mask = ((char)(x+11*y)) % 0x40;
    /* set all bits other than face of interest to 1 */
    printf("face_mask = %x, x,y: %g, %g (orig)\t", face_mask, x, y);
    face_mask = ~(0x03 << (2 * 2));
    printf("face_mask = %x (updated)\n", face_mask);
    face_mask = 0x7e;
    face_mask = 0x73;
    face_mask = 0x4f;
    face_mask = 0;
    #endif /* 0 */

    double lx, ly, lz;
    double ux, uy, uz;
    lx = x - dx; ly = y - dy; lz = z - dz;
    ux = x + dx; uy = y + dy; uz = z + dz;

    #if 0
    /* FOR DEBUGGING PURPOSES ONLY */
    lx = x; ly = y; lz = z;
    ux = x + 2*dx; uy = y + 2*dy; uz = z + 2*dz;
    #endif /* 0 */

    #if 0
    printf("%s:\n", __FUNCTION__);
    printf("cube centered at %g,%g,%g\n", x, y, z);
    printf("cube sizes are %g,%g,%g\n", scaleX, scaleY, scaleZ);
    printf("cube extents are %g,%g,%g to %g,%g,%g\n", lx, ly, lz, ux, uy, uz);
    #endif /* 0 */

    if( (face_mask & (1<<0)) == 0) {
        /* left (-x) */
        trig_list_add(list, lx, ly, lz,
                            lx, ly, uz,
                            lx, uy, uz);
        trig_list_add(list, lx, ly, lz,
                            lx, uy, uz,
                            lx, uy, lz);
    }
    if( (face_mask & (1<<1)) == 0) {
        /* right (+x) */
        trig_list_add(list, ux, ly, lz,
                            ux, uy, uz,
                            ux, ly, uz);
        trig_list_add(list, ux, ly, lz,
                            ux, uy, lz,
                            ux, uy, uz);
    }

    if( (face_mask & (1<<2)) == 0) {
        /* front (-y) */
        trig_list_add(list, lx, ly, lz,
                            ux, ly, uz,
                            lx, ly, uz);
        trig_list_add(list, lx, ly, lz,
                            ux, ly, lz,
                            ux, ly, uz);
    }
    if( (face_mask & (1<<3)) == 0) {
        /* back (+y) */
        trig_list_add(list, lx, uy, lz,
                            lx, uy, uz,
                            ux, uy, uz);
        trig_list_add(list, lx, uy, lz,
                            ux, uy, uz,
                            ux, uy, lz);
    }

    if( (face_mask & (1<<4)) == 0) {
        /* bottom (-z) */
        trig_list_add(list, lx, ly, lz,
                            ux, uy, lz,
                            ux, ly, lz);
        trig_list_add(list, lx, uy, lz,
                            ux, uy, lz,
                            lx, ly, lz);
    }
    if( (face_mask & (1<<5)) == 0) {
        /* top (+z face) */
        trig_list_add(list, lx, ly, uz,
                            ux, ly, uz,
                            ux, uy, uz);
        trig_list_add(list, lx, uy, uz,
                            lx, ly, uz,
                            ux, uy, uz);
    }
}


int maze_add_maze_face(maze_t *maze, int face, trig_list_t *list, maze_output_opts_t *opts) {

    double scale = 1.0;
    int separateMarkers = opts->separateMarkers;

    /* for each cell */
    int rows = maze->faces[face].rows;
    int cols = maze->faces[face].cols;
    for(int row=0; row<rows; ++row) {
        for(int col=0; col<cols; ++col) {
            if( face_get_cell(&maze->faces[face], row, col)!=0 ) {
                /* output small cube for high and low faces */
                int x1=0,y1=0,z1=0;
                x1 = row;
                y1 = col;
                z1 = 0;

                #if 0
                /* FOR DEBUGGING PURPOSES ONLY */
                /* spaces out cubes */
                x1 *= 2.0;
                y1 *= 2.0;
                #endif /* 0 */

                /* compute face mask for cube */
                char mask1 = 0;
                if( row>0
                    && face_get_cell(&maze->faces[face], row-1, col) !=0 ) {
                    mask1 |= 1<<0;
                }
                if( row<rows-1
                    && face_get_cell(&maze->faces[face], row+1, col) !=0 ) {
                    mask1 |= 1<<1;
                }
                if( col>0 
                    && face_get_cell(&maze->faces[face], row, col-1) !=0 ) {
                    mask1 |= 1<<2;
                }
                if( col<cols-1
                    && face_get_cell(&maze->faces[face], row, col+1) !=0 ) {
                    mask1 |= 1<<3;
                }

                /* additional masking for border */
                if( row==0 ) {
                    mask1 |= 1<<4 | 1<<0;
                }
                if( row==rows-1 ) {
                    mask1 |= 1<<4 | 1<<1;
                }
                if( col==0 ) {
                    mask1 |= 1<<4 | 1<<2;
                }
                if( col==cols-1 ) {
                    mask1 |= 1<<4 | 1<<3;
                }

                /* export cube, scaling will happen later */
                maze_add_cube(list, x1, y1, z1, mask1, scale, scale, scale);
            }
        }
    }

    if( separateMarkers ) {
        /* add markers to face, to work around slicer glitch (PrusaSlicer 2.8.1) */
        #if 0
        double hiddenRadius = 0.05;
        maze_add_marker2(list, maze, face, maze->startPos, hiddenRadius, scale);
        maze_add_marker1(list, maze, face, maze->endPos, hiddenRadius, scale);
        #endif // ???
        /* re-tag hidden markers as part of the face */
        trig_list_set_groupid(list, FACE_NONE);
    }

    /* re-add end markers again with original group Ids  */
    double markerRadius = 0.1;
    maze_add_marker2(list, maze, face, maze->startPos, markerRadius, scale);
    maze_add_marker1(list, maze, face, maze->endPos, markerRadius, scale);
    if( !separateMarkers ) {
        /* re-tag hidden markers as part of the face */
        trig_list_set_groupid(list, FACE_NONE);
    }

    return 1;
}


int maze_add_maze(maze_t *maze, trig_list_t *list, maze_output_opts_t *opts) {

    /* for each face */
    for(int face=0; face<maze->numFaces; ++face) {

        trig_list_t faceTrigs1, faceTrigs2;
        trig_list_init(&faceTrigs1);
        trig_list_init(&faceTrigs2);
        maze_add_maze_face(maze, face, &faceTrigs1, opts);
        maze_add_maze_face(maze, face, &faceTrigs2, opts);

        /* translate face */
        if( face == 0 ) {
            trig_list_scale(&faceTrigs2, 1.0, 1.0, -1.0);
            trig_list_move(&faceTrigs2, 0.0, 0.0, 0.5);
            trig_list_move(&faceTrigs2, 0.0, 0.0, 0.0);
            trig_list_move(&faceTrigs1, 0.0, 0.0, maze->dimensions[2]-1.0);
        } else if( face == 1 ) {
            trig_list_scale(&faceTrigs2, 1.0, 1.0, -1.0);
            trig_list_rotate_axial(&faceTrigs1, 0, M_PI/2.0);
            trig_list_rotate_axial(&faceTrigs2, 0, M_PI/2.0);
            trig_list_move(&faceTrigs1, 0.0, 0.0, 0.0);
            trig_list_move(&faceTrigs2, 0.0, maze->dimensions[1]-1.0, 0.0);
        } else if( face == 2 ) {
            trig_list_scale(&faceTrigs2, 1.0, 1.0, -1.0);
            trig_list_rotate_axial(&faceTrigs1, 2, M_PI/2.0);
            trig_list_rotate_axial(&faceTrigs2, 2, M_PI/2.0);
            trig_list_rotate_axial(&faceTrigs1, 1, -M_PI/2.0);
            trig_list_rotate_axial(&faceTrigs2, 1, -M_PI/2.0);
            trig_list_move(&faceTrigs2, 0.0, 0.0, 0.0);
            trig_list_move(&faceTrigs1, maze->dimensions[0]-1.0, 0.0, 0.0);
        }
        
        /* update group ids for markers */
        trig_list_replace_groupid(&faceTrigs1, 3*face+1, 3);
        trig_list_replace_groupid(&faceTrigs1, 3*face+2, 4);
        trig_list_replace_groupid(&faceTrigs2, 3*face+1 + 9, 3);
        trig_list_replace_groupid(&faceTrigs2, 3*face+2 + 9, 4);

        /* set group id for the rest of the face */
        trig_list_replace_groupid(&faceTrigs1, 3*face, -1);
        trig_list_replace_groupid(&faceTrigs2, 3*face + 9, -1);

        /* add face to maze list */
        trig_list_concatenate(list, &faceTrigs1);
        trig_list_concatenate(list, &faceTrigs2);

        trig_list_free(&faceTrigs1);
        trig_list_free(&faceTrigs2);
    }

    trig_list_move(list, 0.5, 0.5, 0.5);

    return 0;
}


static int maze_add_flat_border(trig_list_t *list, double xOffset, double yOffset, double xSize, double ySize, double edgeWidth, double scale) {

    #if 0
    printf("%s\n", __FUNCTION__);
    printf("\tedgeWidth: %g\tscale: %g\n", edgeWidth, scale);
    printf("\tx,y: %g,%g; w,h: %g,%g\n", xOffset, yOffset, xSize, ySize);
    #endif /* 0 */

    double xi[4], yi[4], xm[4], ym[4], xo[4], yo[4];

    xo[0] = xOffset;
    xi[0] = xOffset+1;
    yo[0] = yOffset;
    yi[0] = yOffset+1;

    xo[1] = xOffset+xSize;
    xi[1] = xOffset+xSize-1;
    yo[1] = yOffset;
    yi[1] = yOffset+1;

    xo[2] = xOffset+xSize;
    xi[2] = xOffset+xSize-1;
    yo[2] = yOffset+ySize;
    yi[2] = yOffset+ySize-1;

    xo[3] = xOffset;
    xi[3] = xOffset+1;
    yo[3] = yOffset+ySize;
    yi[3] = yOffset+ySize-1;

    double z0 = -0.5*scale + edgeWidth;
    double z1 = 0.5*scale;

    /* adjust scaling */
    for(int i=0; i<4; ++i) {
        xo[i] *= scale;
        yo[i] *= scale;
        xi[i] *= scale;
        yi[i] *= scale;
    }

    if( edgeWidth < 0 ) edgeWidth = 0.0;
    if( edgeWidth >= scale ) edgeWidth = scale;

    xm[0] = xi[0] - edgeWidth;
    ym[0] = yi[0] - edgeWidth;

    xm[1] = xi[1] + edgeWidth;
    ym[1] = yi[1] - edgeWidth;

    xm[2] = xi[2] + edgeWidth;
    ym[2] = yi[2] + edgeWidth;

    xm[3] = xi[3] - edgeWidth;
    ym[3] = yi[3] + edgeWidth;

    for(int i=0; i<4; ++i) {
        int j = (i+1)%4;

        /* outer sloped face */
        #if 1
        trig_list_add(list, xo[i], yo[i], z1,
                            xm[i], ym[i], z0,
                            xm[j], ym[j], z0);
        #endif
        #if 0
        trig_list_add(list, xo[i], yo[i], z1,
                            xm[j], ym[j], z0,
                            xo[j], yo[j], z1);
        #else
        double edgeSize=(11);    // TODO: this needs to be dynamic
        for(int k=0; k<edgeSize; k+=1) {
            double xo_diff = xo[j]-xo[i];
            double yo_diff = yo[j]-yo[i];
            double xo_mid_i = k*xo_diff/edgeSize+xo[i];
            double yo_mid_i = k*yo_diff/edgeSize+yo[i];
            double xo_mid_j = (k+1)*xo_diff/edgeSize+xo[i];
            double yo_mid_j = (k+1)*yo_diff/edgeSize+yo[i];
            //printf("x,y -> %g, %g; %g, %g\n", xo_mid_i, yo_mid_i, xo_mid_j, yo_mid_j);
            trig_list_add(list, xm[j], ym[j], z0,
                                xo_mid_j, yo_mid_j, z1,
                                xo_mid_i, yo_mid_i, z1);
        }
        #endif

        /* flat section for build-plate */
        if( edgeWidth > 0.0 ) {
            #if 0
            trig_list_add(list, xm[i], ym[i], z0,
                                xi[i], yi[i], z0,
                                xi[j], yi[j], z0);
            #else
            double edgeSize=(11-2);    // TODO: this needs to be dynamic
            for(int k=0; k<edgeSize; k+=1) {
                double xi_diff = xi[j]-xi[i];
                double yi_diff = yi[j]-yi[i];
                double xi_mid_i = k*xi_diff/edgeSize+xi[i];
                double yi_mid_i = k*yi_diff/edgeSize+yi[i];
                double xi_mid_j = (k+1)*xi_diff/edgeSize+xi[i];
                double yi_mid_j = (k+1)*yi_diff/edgeSize+yi[i];
                trig_list_add(list, xm[i], ym[i], z0,
                                    xi_mid_i, yi_mid_i, z0,
                                    xi_mid_j, yi_mid_j, z0);
            }
            #endif // 0
            trig_list_add(list, xm[i], ym[i], z0,
                                xi[j], yi[j], z0,
                                xm[j], ym[j], z0);
        }
    }

    return 0;
}


int maze_add_maze_printable_face(maze_t *maze, trig_list_t *list, maze_output_opts_t *opts, int face) {

    double edgeWidth = opts->edgeWidth;
    double scale = opts->scale;
    int realFace = face;
    if( face < 0 ) face = 0;

    int rows = maze->faces[face].rows;
    int cols = maze->faces[face].cols;

    /* add face */
    trig_list_t faceTrigs0, faceTrigs1, faceTrigs2;
    trig_list_init(&faceTrigs0);
    trig_list_init(&faceTrigs1);
    trig_list_init(&faceTrigs2);
    maze_add_maze_face(maze, face, &faceTrigs0, opts);

    /* add edges to faces */
    trig_list_t border0;
    trig_list_init(&border0);

    /* add perimeters */
    maze_add_flat_border(&border0, -0.5, -0.5, cols, rows, edgeWidth/scale, 1.0);
    trig_list_concatenate(&faceTrigs0, &border0);
    trig_list_free(&border0);

    double minZ = 0.0;
    if( edgeWidth > 0.0 ) {
        minZ = -0.5 + edgeWidth/scale;
        trig_list_set_minimum(&faceTrigs0, minZ, 2);

        /* translate face to account for edgeWidth */
        /* i.e., move minZ to z=0.0 */
        if( face == 0 ) {
            trig_list_move(&faceTrigs0, 0.0, 0.0, -minZ);
        }
    }

    trig_list_concatenate(&faceTrigs1, &faceTrigs0);
    trig_list_concatenate(&faceTrigs2, &faceTrigs0);
    trig_list_free(&faceTrigs0);

    /* translate face */
    if( realFace == -1 ) {
        trig_list_scale(&faceTrigs2, -1.0, 1.0, 1.0);
        trig_list_move(&faceTrigs2, 0.0, 0.0, 0.5);
        double deltaZ = maze->startPos[2]-(maze->dimensions[2]-0.75)-0.5;
        printf("deltaZ = %g\n", deltaZ);
        trig_list_move(&faceTrigs2, 2.5*maze->dimensions[0], 1.5*maze->dimensions[1], deltaZ);
    } if( realFace == 0 ) {
        trig_list_move(&faceTrigs1, 0.0, 0.0, 0.5);
        double deltaZ = maze->startPos[2]-(maze->dimensions[2]-0.75)-0.5;
        trig_list_move(&faceTrigs1, 2.5*maze->dimensions[0], -0.5*maze->dimensions[1], deltaZ);
    } else if( face == 1 ) {
        trig_list_scale(&faceTrigs2, 1.0, 1.0, -1.0);
        trig_list_rotate_axial(&faceTrigs1, 0, M_PI/2.0);
        trig_list_rotate_axial(&faceTrigs2, 0, M_PI/2.0);
        trig_list_move(&faceTrigs1, 0.0, 0.0, 0.0);
        trig_list_move(&faceTrigs2, 0.0, maze->dimensions[1]-1.0, 0.0);
    } else if( face == 2 ) {
        trig_list_scale(&faceTrigs2, 1.0, 1.0, -1.0);
        trig_list_rotate_axial(&faceTrigs1, 2, M_PI/2.0);
        trig_list_rotate_axial(&faceTrigs2, 2, M_PI/2.0);
        trig_list_rotate_axial(&faceTrigs1, 1, -M_PI/2.0);
        trig_list_rotate_axial(&faceTrigs2, 1, -M_PI/2.0);
        trig_list_move(&faceTrigs2, 0.0, 0.0, 0.0);
        trig_list_move(&faceTrigs1, maze->dimensions[0]-1.0, 0.0, 0.0);
    }

    /* add face to maze list */
    if( realFace >= 0 )
        trig_list_concatenate(list, &faceTrigs1);
    if( realFace != 0 )
        trig_list_concatenate(list, &faceTrigs2);

    trig_list_free(&faceTrigs1);
    trig_list_free(&faceTrigs2);

    trig_list_move(list, 0.5, 0.5, 0.5);

    return 0;
}

static int maze_add_maze_slider(maze_t *maze, trig_list_t *list) {
    double scale = 1.0;

    /* add slider */
    double startX = maze->startPos[0] + 0.5;
    double startY = maze->startPos[1] + 0.5;
    double startZ = maze->startPos[2] + 0.5;
    double sizeX = 2*maze->dimensions[0] - 1.5;
    double sizeY = 2*maze->dimensions[1] - 1.5;
    double sizeZ = 2*maze->dimensions[2] - 1.5;
    double scale2 = scale*0.85;
    maze_add_cube(list, startX, startY, startZ, 0, sizeX*scale, scale2, scale2);
    maze_add_cube(list, startX, startY, startZ, 0, scale2, sizeY*scale, scale2);
    maze_add_cube(list, startX, startY, startZ, 0, scale2, scale2, sizeZ*scale);

    trig_list_move(list, 0.5, 0.5, 0.5);

    return 0;
}


static int maze_add_flat_slider(trig_list_t *list, double xOffset, double yOffset, double xSize, double ySize, double scale) {

    double size = 0.9*scale;
    double size2 = size/2.0;
    double xPos = xOffset*scale;
    double yPos = yOffset*scale;
    double zPos = (size-scale)/2.0;

    /* double sizes, so it doesn't fall out */
    xSize *= 2.0;
    ySize *= 2.0;

    /* 'y' direction */
    int dir = 1;
    if( ySize < 0.0 )
        dir = -1;
    ySize = fabs(ySize);
    maze_add_cube(list, xPos, yPos+dir*(scale*ySize/4.0+size2), zPos,
                            0, size, scale*ySize/2.0, size);
    maze_add_cube(list, xPos, yPos+dir*size2/2.0, zPos,
                            0, size, size2, size);

    /* 'x' directions */
    maze_add_cube(list, xPos+scale*(xSize-1)/4.0+size2, yPos, zPos,
                            0, scale*(xSize-1)/2, size, size);
    maze_add_cube(list, xPos-scale*(xSize-1)/4.0-size2, yPos, zPos,
                            0, scale*(xSize-1)/2, size, size);

    return 0;
}


int maze_add_maze_flat(maze_t *maze, trig_list_t *list, maze_output_opts_t *opts) {

    double scale = opts->scale;
    double edgeWidth = opts->edgeWidth;
    char stlName[PATH_MAX];

    if( opts->flatPackDirName ) {
        if( mkdir(opts->flatPackDirName, 0755) ) {
            fprintf(stderr, "mkdir: %s\n", strerror(errno));
        }
    }

    /* write faces */
    /* for each face */
    for(int face=0; face<maze->numFaces; ++face) {
        int rows = maze->faces[face].rows;
        int cols = maze->faces[face].cols;

        trig_list_t faceTrigs1, faceTrigs2;
        trig_list_init(&faceTrigs1);
        trig_list_init(&faceTrigs2);
        maze_add_maze_face(maze, face, &faceTrigs1, opts);
        maze_add_maze_face(maze, face, &faceTrigs2, opts);

        /* scale face */
        trig_list_scale(&faceTrigs1, scale, scale, scale);
        trig_list_scale(&faceTrigs2, -scale, scale, scale);

        trig_list_t border1, border2;
        trig_list_init(&border1);
        trig_list_init(&border2);

        /* add perimeters */
        maze_add_flat_border(&border1, -0.5, -0.5, cols, rows, edgeWidth, scale);
        maze_add_flat_border(&border2, -0.5, -0.5, cols, rows, edgeWidth, scale);

        /* move border2 into proper location */
        trig_list_move(&border2, (-cols+1) * scale, 0.0, 0.0);

        trig_list_concatenate(&faceTrigs1, &border1);
        trig_list_concatenate(&faceTrigs2, &border2);

        /* compute base offsets for face */
        double xBase1 = -cols-1.0;
        double yBase1 = 0.0;
        for(int f2 = 0; f2<face; ++f2) {
            yBase1 += maze->faces[f2].rows+3;
        }
        double xBase2 = cols+1;
        double yBase2 = yBase1;

        /* translate face */
        trig_list_move(&faceTrigs1, xBase1*scale, yBase1*scale, -edgeWidth);
        trig_list_move(&faceTrigs2, xBase2*scale, yBase2*scale, -edgeWidth);

        double minZ = -0.5*scale + edgeWidth;
        trig_list_set_minimum(&faceTrigs1, minZ, 2);
        trig_list_set_minimum(&faceTrigs2, minZ, 2);
        
        if( opts->flatPackDirName != NULL ) {
            /* write to files */
            snprintf(stlName, sizeof(stlName), "%s/Face_%da.stl",
                opts->flatPackDirName, face+1);
            trig_list_write_stl(&faceTrigs1, stlName, "MazeCubeFace");
            snprintf(stlName, sizeof(stlName), "%s/Face_%db.stl",
                opts->flatPackDirName, face+1);
            trig_list_write_stl(&faceTrigs2, stlName, "MazeCubeFace");
        }

        /* add face to maze list */
        trig_list_concatenate(list, &faceTrigs1);
        trig_list_concatenate(list, &faceTrigs2);

        trig_list_free(&faceTrigs1);
        trig_list_free(&faceTrigs2);
        trig_list_free(&border1);
        trig_list_free(&border2);
    }

    /* write slider pieces */
    trig_list_t slider;
    trig_list_init(&slider);
    int xSize, ySize, zSize;
    xSize = maze->dimensions[0];
    ySize = maze->dimensions[1];
    zSize = maze->dimensions[2];
    maze_add_flat_slider(&slider, 0.0, ySize+1.0, xSize+1, -(ySize+1), scale);
    maze_add_flat_slider(&slider, 0.0, ySize+zSize+4.0, zSize+1, ySize+1, scale);

    /* move slider to account for edgeWidth */
    trig_list_move(&slider, 0.0, 0.0, edgeWidth);
    if( opts->flatPackDirName != NULL ) {
        /* write to file */
        snprintf(stlName, sizeof(stlName), "%s/slider.stl",
                opts->flatPackDirName);
        trig_list_write_stl(&slider, stlName, "MazeCubeSlider");
    }
    trig_list_concatenate(list, &slider);
    trig_list_free(&slider);

    trig_list_move(list, 0.5, 0.5, 0.5);

    return 0;
}

int maze_add_solution(maze_t *maze, trig_list_t *list) {
    double scale = 1.0;
    for(int i=0; i<maze->solution.num; ++i) {
        double x, y, z;
        x = maze->solution.positions[i][0]+0.5*scale;
        y = maze->solution.positions[i][1]+0.5*scale;
        z = maze->solution.positions[i][2]+0.5*scale;
        maze_add_cube(list, x, y, z, 0, scale, scale, scale);
    }

    return 1;
}


int maze_export_stl(maze_t *maze, char *filename, maze_output_opts_t *opts) {
    if( maze->numDimensions != 3 ) {
        fprintf(stderr,"%s: STL export is only supported for 3D mazes.\n", __FUNCTION__);
        return -1;
    }

    double scale = opts->scale;

    trig_list_t trigs;
    trig_list_init(&trigs);
    maze_add_maze(maze, &trigs, opts);
    maze_add_maze_slider(maze, &trigs);

    /* scale output */
    trig_list_scale(&trigs, scale, scale, scale);

    /* write to file */
    trig_list_write_stl(&trigs, filename, "MazeCube");

    /* free triangle list */
    trig_list_free(&trigs);

    return 0;
}


int maze_export_stl_printable(maze_t *maze, char *dirname, maze_output_opts_t *opts) {
    if( maze->numDimensions != 3 ) {
        fprintf(stderr,"%s: STL export is only supported for 3D mazes.\n", __FUNCTION__);
        return -1;
    }

    double scale = opts->scale;

    /* construct filename */
    if( dirname == NULL )   dirname = ".";
    else                    mkdir(dirname, 0700);

    trig_list_t trigs;
    trig_list_init(&trigs);
    char filename[PATH_MAX];

    /* write slider */
    snprintf(filename, sizeof(filename), "%s/slider.stl", dirname);
    maze_add_maze_slider(maze, &trigs);
    trig_list_move(&trigs, -0.5, -0.5, -0.5);
    trig_list_scale(&trigs, scale, scale, scale);
    trig_list_write_stl(&trigs, filename, "Slider");
    trig_list_free(&trigs);

    char name[64];
    char name1[64];
    char name2[64];
    char filename1[PATH_MAX];
    char filename2[PATH_MAX];
    for(int face=-1; face<maze->numFaces; ++face) {
        trig_list_init(&trigs);

        if( face > 0 ) {
            snprintf(filename, sizeof(filename), "%s/face_%i.stl", dirname, face);
            snprintf(filename1, sizeof(filename1), "%s/face_%i_marker_1.stl", dirname, face);
            snprintf(filename2, sizeof(filename2), "%s/face_%i_marker_2.stl", dirname, face);
            snprintf(name, sizeof(name), "Face%i", face);
            snprintf(name1, sizeof(name1), "Face%iMarker1", face);
            snprintf(name2, sizeof(name2), "Face%iMarker2", face);
        } else if( face == 0 ) {
            // Face -1 (i.e. end)
            snprintf(filename, sizeof(filename), "%s/end1.stl", dirname);
            snprintf(filename1, sizeof(filename1), "%s/end1_marker_1.stl", dirname);
            snprintf(filename2, sizeof(filename2), "%s/end1_marker_2.stl", dirname);
            snprintf(name, sizeof(name), "End1");
            snprintf(name1, sizeof(name1), "End1Marker1");
            snprintf(name2, sizeof(name2), "End1Marker2");
        } else if( face < 0 ) {
            // Face -1 (i.e. end)
            snprintf(filename, sizeof(filename), "%s/end2.stl", dirname);
            snprintf(filename1, sizeof(filename1), "%s/end2_marker_1.stl", dirname);
            snprintf(filename2, sizeof(filename2), "%s/end2_marker_2.stl", dirname);
            snprintf(name, sizeof(name), "End2");
            snprintf(name1, sizeof(name1), "End2Marker1");
            snprintf(name2, sizeof(name2), "End2Marker2");
        }
        maze_add_maze_printable_face(maze, &trigs, opts, face);

        /* scale output */
        trig_list_scale(&trigs, scale, scale, scale);

        /* initialize per group id trig lists */
        trig_list_t open_trigs, marker1, marker2, face_trigs;
        trig_list_init(&open_trigs);
        trig_list_init(&face_trigs);
        trig_list_init(&marker1);
        trig_list_init(&marker2);

    int num_open = find_open_edges(&trigs);
    printf("%i open edges in model.\n", num_open);

        /* split by group id into separate files. */
        for(int i=0; i<trigs.num; ++i) {
            trig_t *trig_i = &trigs.trig[i];
            trig_get_normal(trig_i);
            face_group_t groupId = trig_i->groupId;
            switch( groupId ) {
                case FACE_OPEN:
                    trig_list_append(&open_trigs, trig_i);
                    break;
                case FACE_MARKER_1:
                    trig_list_append(&marker1, trig_i);
                    break;
                case FACE_MARKER_2:
                    trig_list_append(&marker2, trig_i);
                    break;
                default:
                    trig_list_append(&face_trigs, trig_i);
                    break;
            }
        }

    num_open = find_open_edges(&face_trigs);
    printf("\t%i open edges in faces model.\n", num_open);

    num_open = find_open_edges(&marker1);
    printf("\t%i open edges in marker1 model.\n", num_open);

    num_open = find_open_edges(&marker2);
    printf("\t%i open edges in marker2 model.\n", num_open);

        /* write to files */
        trig_list_write_stl(&face_trigs, filename, name);
        trig_list_write_stl(&marker1, filename1, name1);
        trig_list_write_stl(&marker2, filename2, name2);
        if( open_trigs.num > 0 ) {
            char name[64];
            char filename[PATH_MAX];
            snprintf(filename, sizeof(filename), "errors_%i.stl", face);
            snprintf(name, sizeof(name), "errors_%i", face);
            trig_list_write_stl(&open_trigs, filename, name);
        }

        /* free triangle list */
        trig_list_free(&marker2);
        trig_list_free(&marker1);
        trig_list_free(&face_trigs);
        trig_list_free(&open_trigs);
        trig_list_free(&trigs);
    }

    return 0;
}


int maze_export_stl_flat(maze_t *maze, char *filename, maze_output_opts_t *opts) {

    if( maze->numDimensions != 3 ) {
        fprintf(stderr,"%s: STL export is only supported for 3D mazes.\n", __FUNCTION__);
        return -1;
    }

    trig_list_t trigs;
    trig_list_init(&trigs);
    maze_add_maze_flat(maze, &trigs, opts);

    if( filename != NULL ) {
        /* write to file */
        trig_list_write_stl(&trigs, filename, "MazeCubeFlat");
    }

    /* free triangle list */
    trig_list_free(&trigs);

    return 0;
}


int maze_export_stl_solution(maze_t *maze, char *filename, maze_output_opts_t *opts) {
    if( maze->numDimensions != 3 ) {
        fprintf(stderr,"%s: STL export is only supported for 3D mazes.\n", __FUNCTION__);
        return -1;
    }

    double scale = opts->scale;

    trig_list_t trigs;
    trig_list_init(&trigs);
    maze_add_solution(maze, &trigs);

    /* scale output */
    trig_list_scale(&trigs, scale, scale, scale);

    /* write to file */
    trig_list_write_stl(&trigs, filename, "MazeCubeSolution");

    /* free triangle list */
    trig_list_free(&trigs);

    return 1;
}


/* Allocate opaque triangle list to be handed to the caller */
void *maze_export_trig_list(maze_t *maze) {

    trig_list_t* trigs = calloc(1, sizeof(trig_list_t));
    trig_list_init(trigs);

    /* set default values */
    maze_output_opts_t maze_opts;
    maze_opts.scale = 1.0;
    maze_opts.edgeWidth = 0.4;
    maze_opts.separateMarkers = 1;
    maze_add_maze(maze, trigs, &maze_opts);

    return ((void*)trigs);
}


/* Free opaque triangle list */
void maze_export_trig_list_free(void *list) {
    if( list == NULL )  return;
    trig_list_free((trig_list_t*)list);
    return;
}


/* Get number of spatial dimensions, which will always be 3, for now. */
int maze_export_num_dims(void *list) {
    if( list == NULL )  return -1;
    return 3;   /* hard coded for now, may change later */
}


/* Return the number of triangles in a list. */
int maze_export_num_trigs(void *list) {
    if( list == NULL )  return 0;
    return ((trig_list_t*)list)->num;
}


/* Return a single component of a triangle vertex. */
float maze_export_vertex_dim(void *list, int trig, int vertex, int dim) {
    if( list == NULL )  return -1.0;
    if( trig < 0 || trig > ((trig_list_t*)list)->num )  return -2.0;
    if( vertex < 0 || vertex >= 3 )  return -3.0;
    if( dim < 0 || dim >= 3 )  return -4.0;

    trig_t* trig_ptr = &((trig_list_t*)list)->trig[trig];
    switch(dim) {
    case 0:
        return trig_ptr->x[vertex];
        break;
    case 1:
        return trig_ptr->y[vertex];
        break;
    case 2:
        return trig_ptr->z[vertex];
        break;
    default:
        return -5.0;
        break;
    }
}

/* Return a single component of a triangle vertex normal. */
float maze_export_normal_dim(void *list, int trig, int vertex, int dim) {
    if( list == NULL )  return -1.0;
    if( trig < 0 || trig > ((trig_list_t*)list)->num )  return -2.0;
    if( vertex < 0 || vertex >= 3 )  return -3.0;
    if( dim < 0 || dim >= 3 )  return -4.0;

    trig_t* trig_ptr = &((trig_list_t*)list)->trig[trig];
    switch(dim) {
    case 0:
        return trig_ptr->nx[vertex];
        break;
    case 1:
        return trig_ptr->ny[vertex];
        break;
    case 2:
        return trig_ptr->nz[vertex];
        break;
    default:
        return -5.0;
        break;
    }
}


/* return the group id of a triangle */
int maze_export_groupid(void *list, int trig) {
    if( list == NULL )  return -1.0;
    if( trig < 0 || trig > ((trig_list_t*)list)->num )  return -2.0;
    
    return ((trig_list_t*)list)->trig[trig].groupId;
}
