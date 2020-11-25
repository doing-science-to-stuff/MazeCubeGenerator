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

    /* compute normal for triangle defines by coordinates
     * see: https://mathworld.wolfram.com/CrossProduct.html
     * Equation 2 */
    *nx = uy*vz - uz*vy;
    *ny = uz*vx - ux*vz;
    *nz = ux*vy - uy*vx;
}

static void maze_export_stl_marker1(FILE *fp, maze_t *maze, int face, position_t pos, double radius, double scale, int dir) {
    int d1 = maze->faces[face].d1;
    int d2 = maze->faces[face].d2;

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

    /* circle marker */
    position_t pos1 = calloc(maze->numDimensions,sizeof(int));
    position_t pos2 = calloc(maze->numDimensions,sizeof(int));
    const int numSegsI = 64;
    const int numSegsJ = 8;
    for(int i=0; i<numSegsI; ++i) {
        /* compute points along circle */
        double thetaIm1 =  2.0 * M_PI * (i-1) / numSegsI;
        double xm1 = cos(thetaIm1);
        double ym1 = sin(thetaIm1);
        double thetaI1 =  2.0 * M_PI * i / numSegsI;
        double x1 = cos(thetaI1)*(1-radius);
        double y1 = sin(thetaI1)*(1-radius);
        double thetaI2 =  2.0 * M_PI * (i+1) / numSegsI;
        double x2 = cos(thetaI2)*(1-radius);
        double y2 = sin(thetaI2)*(1-radius);

        /* position cicle points onto face */
        int col1 = y1+c+0.5;
        int col2 = y2+c+0.5;
        int row1 = x1+r+0.5;
        int row2 = x2+r+0.5;
        int rowm1 = xm1+r+0.5;
        int colm1 = ym1+c+0.5;

        /* check of cell under each point is set */
        int cell1 = face_get_cell(&maze->faces[face], row1, col1);
        int cell2 = face_get_cell(&maze->faces[face], row2, col2);
        int cellm1 = face_get_cell(&maze->faces[face], rowm1, colm1);

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

            /* transform into model space */
            maze_export_stl_transform(maze, face, scale, dir, &x11, &y11, &z11);
            maze_export_stl_transform(maze, face, scale, dir, &x12, &y12, &z12);
            maze_export_stl_transform(maze, face, scale, dir, &x21, &y21, &z21);
            maze_export_stl_transform(maze, face, scale, dir, &x22, &y22, &z22);

            /* record first point on surface for end caps */
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

static void maze_export_stl_marker2(FILE *fp, maze_t *maze, int face, double radius, double scale, int dir) {
    /* square marker */
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
    if( (face_mask & (1<<5)) == 0) {
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
}


int maze_export_stl(maze_t *maze, char *filename) {
    if( maze->numDimensions != 3 ) {
        fprintf(stderr,"%s: STL export is only supported for 3D mazes.\n", __FUNCTION__);
        return -1;
    }

    /* open file */
    FILE *fp = fopen(filename,"w");

    /* open solid */
    fprintf(fp,"solid puzzle\n");

    double scale = 1.0;
    /* for each face */
    for(int face=0; face<maze->numFaces; ++face) {
        int d1 = maze->faces[face].d1;
        int d2 = maze->faces[face].d2;
        /* for each cell */
        for(int row=0; row<maze->faces[face].rows; ++row) {
            for(int col=0; col<maze->faces[face].cols; ++col) {
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
                    char mask2 = 0;
#if 0
                    int d1 = maze->faces[face].d1;
                    int d2 = maze->faces[face].d2;
                    int rows = maze->faces[face].rows;
                    int cols = maze->faces[face].cols;
                    if( row > 0
                        && face_get_cell(&maze->faces[face], row-1, col) !=0 ) {
                        mask1 |= 1<<(2*d1);
                    }
                    if( row < rows-1
                        && face_get_cell(&maze->faces[face], row+1, col) !=0 ) {
                        mask1 |= 1<<(2*d1+1);
                    }
                    if( col > 0
                        && face_get_cell(&maze->faces[face], row, col-1) !=0 ) {
                        mask1 |= 1<<(2*d2);
                    }
                    if( col < cols-1
                        && face_get_cell(&maze->faces[face], row, col+1) !=0 ) {
                        mask1 |= 1<<(2*d2);
                    }
                    mask2 = mask1;
#endif /* 1 */
#if 0
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
#endif /* 1 */

                    #if 1
                    /* export cubes */
                    maze_export_stl_cube(fp, x1, y1, z1, mask1, scale, scale, scale);
                    maze_export_stl_cube(fp, x2, y2, z2, mask2, scale, scale, scale);
                    #endif /* 1 */
                }
            }
        }

        /* add end markers */
        double markerRadius = 0.1;
        maze_export_stl_marker1(fp, maze, face, maze->startPos, markerRadius, scale, -1);
        maze_export_stl_marker1(fp, maze, face, maze->startPos, markerRadius, scale, 1);
        maze_export_stl_marker1(fp, maze, face, maze->endPos, markerRadius, scale, -1);
        maze_export_stl_marker1(fp, maze, face, maze->endPos, markerRadius, scale, 1);
    }

    /* add movable piece */
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
    fprintf(fp,"endsolid puzzle\n");

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

    /* open solid */
    fprintf(fp,"solid puzzle\n");

    double scale = 1.0;
    for(int i=0; i<maze->solution.posListNum; ++i) {
        double x, y, z;
        x = maze->solution.positions[i][0];
        y = maze->solution.positions[i][1];
        z = maze->solution.positions[i][2];
        maze_export_stl_cube(fp, x, y, z, 0, scale, scale, scale);
    }

    /* close solid */
    fprintf(fp,"endsolid puzzle\n");

    /* close file */
    fclose(fp); fp=NULL;

    return 0;
}
