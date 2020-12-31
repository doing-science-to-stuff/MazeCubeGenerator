#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "maze.h"

static void maze_export_stl_triangle(FILE *fp,
                                    double x1, double y1, double z1,
                                    double x2, double y2, double z2,
                                    double x3, double y3, double z3,
                                    double nx, double ny, double nz,
                                    int rev) {
    if( rev <= 0 ) {
        fprintf(fp, "facet normal %g %g %g\n", nx, ny, nz);
        fprintf(fp, "  outer loop\n");
        fprintf(fp, "    vertex %g %g %g\n", x1, y1, z1);
        fprintf(fp, "    vertex %g %g %g\n", x2, y2, z2);
        fprintf(fp, "    vertex %g %g %g\n", x3, y3, z3);
    } else {
        fprintf(fp, "facet normal %g %g %g\n", -nx, -ny, -nz);
        fprintf(fp, "  outer loop\n");
        fprintf(fp, "    vertex %g %g %g\n", x1, y1, z1);
        fprintf(fp, "    vertex %g %g %g\n", x3, y3, z3);
        fprintf(fp, "    vertex %g %g %g\n", x2, y2, z2);
    }
    fprintf(fp, "  endloop\n");
    fprintf(fp, "endfacet\n");
}

static void maze_export_stl_transform(maze_t *maze, int face, double scale, int dir, double *x, double *y, double *z) {
    int d1 = maze->faces[face].d1;
    int d2 = maze->faces[face].d2;

    double tx=-1.0, ty=-1.0, tz=-1.0;
    if( d1 == 0 && d2 == 1) {
        tx = *x * scale;
        ty = *y * scale;
        if( dir < 0 ) {
            tz = -(*z) - 0.5*scale;
        } else {
            tz = *z + (maze->dimensions[2]-0.5) * scale;
        }
    } else if( d1 == 0 && d2 == 2 ) {
        tx = *x * scale;
        tz = *y * scale;
        if( dir < 0 ) {
            ty = -(*z) - 0.5*scale;
        } else {
            ty = *z + (maze->dimensions[1]-0.5) * scale;
        }
    } else if( d1 == 1 && d2 == 2 ) {
        ty = *x * scale;
        tz = *y * scale;
        if( dir < 0 ) {
            tx = -(*z) - 0.5*scale;
        } else {
            tx = *z + (maze->dimensions[0]-0.5) * scale;
        }
    } else {
        fprintf(stderr, "Unhandled dimension combination d1=%i,d2=%i\n", d1, d2);
    }

    /* write back to original variables */
    *x  = tx;
    *y  = ty;
    *z  = tz;
}

static void maze_export_get_normal(double x1, double y1, double z1,
                                   double x2, double y2, double z2,
                                   double x3, double y3, double z3,
                                   double *nx, double *ny, double *nz) {
    /* get two edge vectors */
    double ux = x2-x1;
    double uy = y2-y1;
    double uz = z2-z1;
    double vx = x3-x1;
    double vy = y3-y1;
    double vz = z3-z1;

    /* compute normal for triangle defined by coordinates
     * see: https://mathworld.wolfram.com/CrossProduct.html
     * Equation 2 */
    *nx = uy*vz - uz*vy;
    *ny = uz*vx - ux*vz;
    *nz = ux*vy - uy*vx;
}

/* circular marker */
static void maze_export_stl_marker1(FILE *fp, maze_t *maze, int face, position_t pos, double radius, double scale, int dir) {
    int d1 = maze->faces[face].d1;
    int d2 = maze->faces[face].d2;

    /* get row and column for marker */
    int c = pos[d2];
    int r = pos[d1];

    /* some faces need normals inverted */
    int reverse = 0;
    if( face == 0 && dir > 0 )
        reverse = 1;
    if( face == 1 && dir < 0 )
        reverse = 1;
    if( face == 2 && dir > 0 )
        reverse = 1;

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

        double x10 = 0.0, y10 = 0.0, z10 = 0.0;
        double x20 = 0.0, y20 = 0.0, z20 = 0.0;
        double nx, ny, nz;
        for(int j=0; j<numSegsJ; ++j) {
            double thetaJ1 =  M_PI * j / numSegsJ;
            double thetaJ2 =  M_PI * (j+1) / numSegsJ;

            /* compute raw torus coordinates */
            double x11 = cos(thetaI1) * (1+radius*cos(thetaJ1))+r;
            double y11 = sin(thetaI1) * (1+radius*cos(thetaJ1))+c;
            double z11 = radius*sin(thetaJ1);

            double x12 = cos(thetaI1) * (1+radius*cos(thetaJ2))+r;
            double y12 = sin(thetaI1) * (1+radius*cos(thetaJ2))+c;
            double z12 = radius*sin(thetaJ2);

            double x21 = cos(thetaI2) * (1+radius*cos(thetaJ1))+r;
            double y21 = sin(thetaI2) * (1+radius*cos(thetaJ1))+c;
            double z21 = radius*sin(thetaJ1);

            double x22 = cos(thetaI2) * (1+radius*cos(thetaJ2))+r;
            double y22 = sin(thetaI2) * (1+radius*cos(thetaJ2))+c;
            double z22 = radius*sin(thetaJ2);

            /* TODO: When crossing a cell boundary, end cap points should be
             * interpolated to exactly match the boundary. */

            /* transform into model space */
            maze_export_stl_transform(maze, face, scale, dir, &x11, &y11, &z11);
            maze_export_stl_transform(maze, face, scale, dir, &x12, &y12, &z12);
            maze_export_stl_transform(maze, face, scale, dir, &x21, &y21, &z21);
            maze_export_stl_transform(maze, face, scale, dir, &x22, &y22, &z22);

            /* record first point along surface for end caps */
            if( j == 0 && (!cell1 || !cell2) ) {
                x10 = x11;
                y10 = y11;
                z10 = z11;
                x20 = x21;
                y20 = y21;
                z20 = z21;
            }

            /* output transformed triangles */
            if( cell1 && cell2 ) {
                /* curved surface of marker */
                maze_export_get_normal(x11, y11, z11,
                                       x12, y12, z12,
                                       x22, y22, z22,
                                       &nx, &ny, &nz);
                maze_export_stl_triangle(fp, x11, y11, z11,
                        x12, y12, z12,
                        x22, y22, z22,
                        nx, ny, nz, reverse);
                maze_export_stl_triangle(fp, x11, y11, z11,
                        x22, y22, z22,
                        x21, y21, z21,
                        nx, ny, nz, reverse);
            } else if( cell1 && !cell2 && j>0 ) {
                /* cap one end */
                maze_export_get_normal(x10, y10, z10,
                                       x11, y11, z11,
                                       x12, y12, z12,
                                       &nx, &ny, &nz);
                maze_export_stl_triangle(fp, x10, y10, z10,
                                             x11, y11, z11,
                                             x12, y12, z12,
                                             nx, ny, nz, reverse);
            } else if( !cell1 && cell2 && j>0 ) {
                /* cap other end */
                maze_export_get_normal(x20, y20, z20,
                                       x22, y22, z22,
                                       x21, y21, z21,
                                       &nx, &ny, &nz);
                maze_export_stl_triangle(fp, x20, y20, z20,
                                             x22, y22, z22,
                                             x21, y21, z21,
                                             nx, ny, nz, reverse);
            }
        }
    }
}

/* square marker */
static void maze_export_stl_corner(FILE *fp, maze_t *maze, int r, int c, int dr, int dc, int face, double radius, double scale, int dir, int cCap, int rCap) {

    /* some faces need normals inverted */
    int reverse = 0;
    if( face == 0 && dir > 0 )
        reverse = 1;
    if( face == 1 && dir < 0 )
        reverse = 1;
    if( face == 2 && dir > 0 )
        reverse = 1;

    if( dc*dr < 0 )
        reverse ^= 1;

    double x10 = 0.0, y10 = 0.0, z10 = 0.0;
    double x30 = 0.0, y30 = 0.0, z30 = 0.0;
    double nx, ny, nz;

    int numSegs = 32;
    for(int i=0; i<numSegs; ++i) {
        double thetaI1 = M_PI * i / numSegs;
        double thetaI2 = M_PI * (i+1) / numSegs;

        /* compute raw segment coordinates */
        double x11 = r+dr+dr*radius*cos(thetaI1);
        double y11 = c+0.5*dc;
        double z11 = radius*sin(thetaI1);

        double x12 = r+dr+dr*radius*cos(thetaI2);
        double y12 = c+0.5*dc;
        double z12 = radius*sin(thetaI2);

        double x21 = r+dr+dr*radius*cos(thetaI1);
        double y21 = c+dc+dc*radius*cos(thetaI1);
        double z21 = radius*sin(thetaI1);

        double x22 = r+dr+dr*radius*cos(thetaI2);
        double y22 = c+dc+dc*radius*cos(thetaI2);
        double z22 = radius*sin(thetaI2);

        double x31 = r+0.5*dr;
        double y31 = c+dc+dc*radius*cos(thetaI1);
        double z31 = radius*sin(thetaI1);

        double x32 = r+0.5*dr;
        double y32 = c+dc+dc*radius*cos(thetaI2);
        double z32 = radius*sin(thetaI2);

        /* transform points */
        maze_export_stl_transform(maze, face, scale, dir, &x11, &y11, &z11);
        maze_export_stl_transform(maze, face, scale, dir, &x12, &y12, &z12);
        maze_export_stl_transform(maze, face, scale, dir, &x21, &y21, &z21);
        maze_export_stl_transform(maze, face, scale, dir, &x22, &y22, &z22);
        maze_export_stl_transform(maze, face, scale, dir, &x31, &y31, &z31);
        maze_export_stl_transform(maze, face, scale, dir, &x32, &y32, &z32);

        /* record first point on surface for end caps */
        if( i == 0 ) {
            x10 = x11;
            y10 = y11;
            z10 = z11;
            x30 = x31;
            y30 = y31;
            z30 = z31;
        }

        /* curved surface of marker */
        maze_export_get_normal(x11, y11, z11,
                x12, y12, z12,
                x22, y22, z22,
                &nx, &ny, &nz);
        maze_export_stl_triangle(fp, x11, y11, z11,
                x12, y12, z12,
                x22, y22, z22,
                nx, ny, nz, reverse);
        maze_export_stl_triangle(fp, x11, y11, z11,
                x22, y22, z22,
                x21, y21, z21,
                nx, ny, nz, reverse);

        maze_export_get_normal(x31, y31, z31,
                x22, y22, z22,
                x32, y32, z32,
                &nx, &ny, &nz);
        maze_export_stl_triangle(fp, x31, y31, z31,
                x22, y22, z22,
                x32, y32, z32,
                nx, ny, nz, reverse);
        maze_export_stl_triangle(fp, x31, y31, z31,
                x21, y21, z21,
                x22, y22, z22,
                nx, ny, nz, reverse);

        /* add end caps */
        if( rCap==0 ) {
            maze_export_get_normal(x11, y11, z11,
                    x10, y10, z10,
                    x12, y12, z12,
                    &nx, &ny, &nz);
            maze_export_stl_triangle(fp, x11, y11, z11,
                    x10, y10, z10,
                    x12, y12, z12,
                    nx, ny, nz, reverse);
        }

        if( cCap==0 ) {
            maze_export_get_normal(x31, y31, z31,
                    x32, y32, z32,
                    x30, y30, z30,
                    &nx, &ny, &nz);
            maze_export_stl_triangle(fp, x31, y31, z31,
                    x32, y32, z32,
                    x30, y30, z30,
                    nx, ny, nz, reverse);
        }
    }
}

static void maze_export_stl_edge1(FILE *fp, maze_t *maze, int r, int c, int dr, int face, double radius, double scale, int dir) {

    /* some faces need normals inverted */
    int reverse = 0;
    if( face == 0 && dir > 0 )
        reverse = 1;
    if( face == 1 && dir < 0 )
        reverse = 1;
    if( face == 2 && dir > 0 )
        reverse = 1;

    if( dr < 0 )
        reverse ^= 1;

    double nx, ny, nz;

    int numSegs = 32;
    for(int i=0; i<numSegs; ++i) {
        double thetaI1 = M_PI * i / numSegs;
        double thetaI2 = M_PI * (i+1) / numSegs;

        /* compute raw segment coordinates */
        double x11 = r+dr+dr*radius*cos(thetaI1);
        double y11 = c+0.5;
        double z11 = radius*sin(thetaI1);

        double x12 = r+dr+dr*radius*cos(thetaI2);
        double y12 = c+0.5;
        double z12 = radius*sin(thetaI2);

        double x21 = r+dr+dr*radius*cos(thetaI1);
        double y21 = c-0.5;
        double z21 = radius*sin(thetaI1);

        double x22 = r+dr+dr*radius*cos(thetaI2);
        double y22 = c-0.5;
        double z22 = radius*sin(thetaI2);

        /* transform points */
        maze_export_stl_transform(maze, face, scale, dir, &x11, &y11, &z11);
        maze_export_stl_transform(maze, face, scale, dir, &x12, &y12, &z12);
        maze_export_stl_transform(maze, face, scale, dir, &x21, &y21, &z21);
        maze_export_stl_transform(maze, face, scale, dir, &x22, &y22, &z22);

        /* outer shell of marker */
        maze_export_get_normal(x11, y11, z11,
                x22, y22, z22,
                x12, y12, z12,
                &nx, &ny, &nz);
        maze_export_stl_triangle(fp, x11, y11, z11,
                x22, y22, z22,
                x12, y12, z12,
                nx, ny, nz, reverse);
        maze_export_stl_triangle(fp, x11, y11, z11,
                x21, y21, z21,
                x22, y22, z22,
                nx, ny, nz, reverse);
    }
}

static void maze_export_stl_edge2(FILE *fp, maze_t *maze, int r, int c, int dc, int face, double radius, double scale, int dir) {

    /* some faces need normals inverted */
    int reverse = 0;
    if( face == 0 && dir > 0 )
        reverse = 1;
    if( face == 1 && dir < 0 )
        reverse = 1;
    if( face == 2 && dir > 0 )
        reverse = 1;

    if( dc < 0 )
        reverse ^= 1;

    double nx, ny, nz;

    int numSegs = 32;
    for(int i=0; i<numSegs; ++i) {
        double thetaI1 = M_PI * i / numSegs;
        double thetaI2 = M_PI * (i+1) / numSegs;

        /* compute raw segment coordinates */
        double x11 = r+0.5;
        double y11 = c+dc+dc*radius*cos(thetaI1);
        double z11 = radius*sin(thetaI1);

        double x12 = r+0.5;
        double y12 = c+dc+dc*radius*cos(thetaI2);
        double z12 = radius*sin(thetaI2);

        double x21 = r-0.5;
        double y21 = c+dc+dc*radius*cos(thetaI1);
        double z21 = radius*sin(thetaI1);

        double x22 = r-0.5;
        double y22 = c+dc+dc*radius*cos(thetaI2);
        double z22 = radius*sin(thetaI2);

        /* transform points */
        maze_export_stl_transform(maze, face, scale, dir, &x11, &y11, &z11);
        maze_export_stl_transform(maze, face, scale, dir, &x12, &y12, &z12);
        maze_export_stl_transform(maze, face, scale, dir, &x21, &y21, &z21);
        maze_export_stl_transform(maze, face, scale, dir, &x22, &y22, &z22);

        /* outer shell of marker */
        maze_export_get_normal(x11, y11, z11,
                x12, y12, z12,
                x22, y22, z22,
                &nx, &ny, &nz);
        maze_export_stl_triangle(fp, x11, y11, z11,
                x12, y12, z12,
                x22, y22, z22,
                nx, ny, nz, reverse);
        maze_export_stl_triangle(fp, x11, y11, z11,
                x22, y22, z22,
                x21, y21, z21,
                nx, ny, nz, reverse);
    }
}

static void maze_export_stl_marker2(FILE *fp, maze_t *maze, int face, position_t pos, double radius, double scale, int dir) {
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

    /* add corners */
    maze_export_stl_corner(fp, maze, r, c, -1, -1, face, radius, scale, dir, tSide, lSide);
    maze_export_stl_corner(fp, maze, r, c, -1, 1, face, radius, scale, dir, bSide, lSide);
    maze_export_stl_corner(fp, maze, r, c, 1, -1, face, radius, scale, dir, tSide, rSide);
    maze_export_stl_corner(fp, maze, r, c, 1, 1, face, radius, scale, dir, bSide, rSide);

    /* add straight segments */
    if( lSide != 0 )
        maze_export_stl_edge1(fp, maze, r, c, -1, face, radius, scale, dir);
    if( rSide != 0 )
        maze_export_stl_edge1(fp, maze, r, c, 1, face, radius, scale, dir);
    if( tSide != 0 )
        maze_export_stl_edge2(fp, maze, r, c, -1, face, radius, scale, dir);
    if( bSide != 0 )
        maze_export_stl_edge2(fp, maze, r, c, 1, face, radius, scale, dir);
}

static void maze_export_stl_cube(FILE *fp, int x, int y, int z, char face_mask, double scaleX, double scaleY, double scaleZ) {
    double dx = scaleX/2.0;
    double dy = scaleY/2.0;
    double dz = scaleZ/2.0;

    if( (face_mask & (1<<0)) == 0) {
        /* left (-x) */
        maze_export_stl_triangle(fp, x-dx, y-dy, z-dz,
                x-dx, y-dy, z+dz,
                x-dx, y+dy, z+dz,
                -1, 0, 0, 0);
        maze_export_stl_triangle(fp, x-dx, y-dy, z-dz,
                x-dx, y+dy, z+dz,
                x-dx, y+dy, z-dz,
                -1, 0, 0, 0);
    }
    if( (face_mask & (1<<1)) == 0) {
        /* right (+x) */
        maze_export_stl_triangle(fp, x+dx, y-dy, z-dz,
                x+dx, y+dy, z+dz,
                x+dx, y-dy, z+dz,
                1, 0, 0, 0);
        maze_export_stl_triangle(fp, x+dx, y-dy, z-dz,
                x+dx, y+dy, z-dz,
                x+dx, y+dy, z+dz,
                1, 0, 0, 0);
    }

    if( (face_mask & (1<<2)) == 0) {
        /* front (-y) */
        maze_export_stl_triangle(fp, x-dx, y-dy, z-dz,
                x+dx, y-dy, z+dz,
                x-dx, y-dy, z+dz,
                0, -1, 0, 0);
        maze_export_stl_triangle(fp, x-dx, y-dy, z-dz,
                x+dx, y-dy, z-dz,
                x+dx, y-dy, z+dz,
                0, -1, 0, 0);
    }
    if( (face_mask & (1<<3)) == 0) {
        /* back (+y) */
        maze_export_stl_triangle(fp, x-dx, y+dy, z-dz,
                x-dx, y+dy, z+dz,
                x+dx, y+dy, z+dz,
                0, 1, 0, 0);
        maze_export_stl_triangle(fp, x-dx, y+dy, z-dz,
                x+dx, y+dy, z+dz,
                x+dx, y+dy, z-dz,
                0, 1, 0, 0);
    }

    if( (face_mask & (1<<4)) == 0) {
        /* bottom (-z) */
        maze_export_stl_triangle(fp, x-dx, y-dy, z-dz,
                x+dx, y+dy, z-dz,
                x+dx, y-dy, z-dz,
                0, 0, -1, 0);
        maze_export_stl_triangle(fp, x-dx, y+dy, z-dz,
                x+dx, y+dy, z-dz,
                x-dx, y-dy, z-dz,
                0, 0, -1, 0);
    }
    if( (face_mask & (1<<5)) == 0) {
        /* top (+z face) */
        maze_export_stl_triangle(fp, x-dx, y-dy, z+dz,
                x+dx, y-dy, z+dz,
                x+dx, y+dy, z+dz,
                0, 0, 1, 0);
        maze_export_stl_triangle(fp, x-dx, y+dy, z+dz,
                x-dx, y-dy, z+dz,
                x+dx, y+dy, z+dz,
                0, 0, 1, 0);
    }
}


int maze_export_stl(maze_t *maze, char *filename) {
    if( maze->numDimensions != 3 ) {
        fprintf(stderr,"%s: STL export is only supported for 3D mazes.\n", __FUNCTION__);
        return -1;
    }

    /* open file */
    FILE *fp = fopen(filename,"w");
    if( fp == NULL ) {
        perror("fopen");
        return -1;
    }

    /* start solid */
    fprintf(fp,"solid MazeCube\n");

    double scale = 1.0;
    /* for each face */
    for(int face=0; face<maze->numFaces; ++face) {
        int d1 = maze->faces[face].d1;
        int d2 = maze->faces[face].d2;
        /* for each cell */
        for(int row=1; row<maze->faces[face].rows-1; ++row) {
            for(int col=1; col<maze->faces[face].cols-1; ++col) {
                if( face_get_cell(&maze->faces[face], row, col)!=0 ) {
                    /* output small cube for high and low faces */
                    int x1=0,y1=0,z1=0;
                    int x2=0,y2=0,z2=0;
                    if( d1 == 0 && d2 == 1 ) {
                        x1 = x2 = row;
                        y1 = y2 = col;
                        z1 = 0;
                        z2 = maze->dimensions[2]-1;
                    } else if( d1 == 0 && d2 == 2 ) {
                        x1 = x2 = row;
                        y1 = 0;
                        y2 = maze->dimensions[1]-1;
                        z1 = z2 = col;
                    } else if( d1 == 1 && d2 == 2 ) {
                        x1 = 0;
                        x2 = maze->dimensions[0]-1;
                        y1 = y2 = row;
                        z1 = z2 = col;
                    } else {
                        fprintf(stderr, "Unhandled dimension combination d1=%i,d2=%i\n", d1, d2);
                    }

                    /* compute face mask for cube */
                    char mask1 = 0;
                    int d1 = maze->faces[face].d1;
                    int d2 = maze->faces[face].d2;
                    int rows = maze->faces[face].rows;
                    int cols = maze->faces[face].cols;
                    if( face_get_cell(&maze->faces[face], row-1, col) !=0 ) {
                        mask1 |= 1<<(2*d1);
                    }
                    if( face_get_cell(&maze->faces[face], row+1, col) !=0 ) {
                        mask1 |= 1<<(2*d1+1);
                    }
                    if( face_get_cell(&maze->faces[face], row, col-1) !=0 ) {
                        mask1 |= 1<<(2*d2);
                    }
                    if( face_get_cell(&maze->faces[face], row, col+1) !=0 ) {
                        mask1 |= 1<<(2*d2+1);
                    }
                    char mask2 = mask1;

                    if( row==0 || col==0
                        || row == rows-1 || col == cols-1 ) {
                        /* is a border cell */
                        if( d1 == 0 && d2 == 1 ) {
                            /* mask z */
                            mask1 |= (1<<4);
                            mask2 |= (1<<5);
                        } else if( d1 == 0 && d2 == 2 ) {
                            /* mask y */
                            mask1 |= (1<<2);
                            mask2 |= (1<<3);
                        } else if( d1 == 1 && d2 == 2 ) {
                            /* mask x */
                            mask1 |= (1<<0);
                            mask2 |= (1<<1);
                        } else {
                            fprintf(stderr, "Unhandled dimension combo! (d1=%i,d2=%i)\n",d1,d2);
                        }
                    }

                    /* export cubes */
                    maze_export_stl_cube(fp, x1, y1, z1, mask1, scale, scale, scale);
                    maze_export_stl_cube(fp, x2, y2, z2, mask2, scale, scale, scale);
                }
            }
        }

        /* add end markers */
        double markerRadius = 0.1;
        maze_export_stl_marker2(fp, maze, face, maze->startPos, markerRadius, scale, -1);
        maze_export_stl_marker2(fp, maze, face, maze->startPos, markerRadius, scale, 1);
        maze_export_stl_marker1(fp, maze, face, maze->endPos, markerRadius, scale, -1);
        maze_export_stl_marker1(fp, maze, face, maze->endPos, markerRadius, scale, 1);
    }

    /* add frame corners */
    int xEnd = maze->dimensions[0]-1;
    int yEnd = maze->dimensions[1]-1;
    int zEnd = maze->dimensions[2]-1;
    maze_export_stl_cube(fp,    0,    0,    0, (1<<1)|(1<<3)|(1<<5), scale, scale, scale);
    maze_export_stl_cube(fp, xEnd,    0,    0, (1<<0)|(1<<3)|(1<<5), scale, scale, scale);
    maze_export_stl_cube(fp,    0, yEnd,    0, (1<<1)|(1<<2)|(1<<5), scale, scale, scale);
    maze_export_stl_cube(fp, xEnd, yEnd,    0, (1<<0)|(1<<2)|(1<<5), scale, scale, scale);
    maze_export_stl_cube(fp,    0,    0, zEnd, (1<<1)|(1<<3)|(1<<4), scale, scale, scale);
    maze_export_stl_cube(fp, xEnd,    0, zEnd, (1<<0)|(1<<3)|(1<<4), scale, scale, scale);
    maze_export_stl_cube(fp,    0, yEnd, zEnd, (1<<1)|(1<<2)|(1<<4), scale, scale, scale);
    maze_export_stl_cube(fp, xEnd, yEnd, zEnd, (1<<0)|(1<<2)|(1<<4), scale, scale, scale);

    /* add frame edges */
    /* x */
    maze_export_stl_cube(fp, xEnd/2,    0,    0, (1<<0)|(1<<1), (xEnd-1)*scale, scale, scale);
    maze_export_stl_cube(fp, xEnd/2, yEnd,    0, (1<<0)|(1<<1), (xEnd-1)*scale, scale, scale);
    maze_export_stl_cube(fp, xEnd/2,    0, zEnd, (1<<0)|(1<<1), (xEnd-1)*scale, scale, scale);
    maze_export_stl_cube(fp, xEnd/2, yEnd, zEnd, (1<<0)|(1<<1), (xEnd-1)*scale, scale, scale);
    /* y */
    maze_export_stl_cube(fp,    0, yEnd/2,    0, (1<<2)|(1<<3), scale, (yEnd-1)*scale, scale);
    maze_export_stl_cube(fp, xEnd, yEnd/2,    0, (1<<2)|(1<<3), scale, (yEnd-1)*scale, scale);
    maze_export_stl_cube(fp,    0, yEnd/2, zEnd, (1<<2)|(1<<3), scale, (yEnd-1)*scale, scale);
    maze_export_stl_cube(fp, xEnd, yEnd/2, zEnd, (1<<2)|(1<<3), scale, (yEnd-1)*scale, scale);
    /* z */
    maze_export_stl_cube(fp,    0,    0, zEnd/2, (1<<4)|(1<<5), scale, scale, (zEnd-1)*scale);
    maze_export_stl_cube(fp, xEnd,    0, zEnd/2, (1<<4)|(1<<5), scale, scale, (zEnd-1)*scale);
    maze_export_stl_cube(fp,    0, yEnd, zEnd/2, (1<<4)|(1<<5), scale, scale, (zEnd-1)*scale);
    maze_export_stl_cube(fp, xEnd, yEnd, zEnd/2, (1<<4)|(1<<5), scale, scale, (zEnd-1)*scale);

    /* add slider */
    double startX = maze->startPos[0];
    double startY = maze->startPos[1];
    double startZ = maze->startPos[2];
    double sizeX = 2*maze->dimensions[0]-0.5;
    double sizeY = 2*maze->dimensions[1]-0.5;
    double sizeZ = 2*maze->dimensions[2]-0.5;
    double scale2 = scale*0.95;
    maze_export_stl_cube(fp, startX, startY, startZ, 0, sizeX*scale, scale2, scale2);
    maze_export_stl_cube(fp, startX, startY, startZ, 0, scale, sizeY*scale2, scale2);
    maze_export_stl_cube(fp, startX, startY, startZ, 0, scale, scale2, sizeZ*scale2);

    /* close solid */
    fprintf(fp,"endsolid MazeCube\n");

    /* close file */
    fclose(fp); fp=NULL;

    return 0;
}

int maze_export_stl_solution(maze_t *maze, char *filename) {
    if( maze->numDimensions != 3 ) {
        fprintf(stderr,"%s: STL export is only supported for 3D mazes.\n", __FUNCTION__);
        return -1;
    }

    /* open file */
    FILE *fp = fopen(filename,"w");
    if( fp == NULL ) {
        perror("fopen");
        return -1;
    }

    /* open solid */
    fprintf(fp,"solid MazeCubeSolution\n");

    double scale = 1.0;
    for(int i=0; i<maze->solution.posListNum; ++i) {
        double x, y, z;
        x = maze->solution.positions[i][0];
        y = maze->solution.positions[i][1];
        z = maze->solution.positions[i][2];
        maze_export_stl_cube(fp, x, y, z, 0, scale, scale, scale);
    }

    /* close solid */
    fprintf(fp,"endsolid MazeCubeSolution\n");

    /* close file */
    fclose(fp); fp=NULL;

    return 0;
}
