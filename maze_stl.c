#include <stdio.h>
#include "maze.h"

static void maze_export_stl_triangle(FILE *fp,
                                    double x1, double y1, double z1,
                                    double x2, double y2, double z2,
                                    double x3, double y3, double z3,
                                    double nx, double ny, double nz) {
    fprintf(fp, "facet normal %g %g %g\n", nx, ny, nz);
    fprintf(fp, "  outer loop\n");
    fprintf(fp, "    vertex %g %g %g\n", x1, y1, z1);
    fprintf(fp, "    vertex %g %g %g\n", x2, y2, z2);
    fprintf(fp, "    vertex %g %g %g\n", x3, y3, z3);
    fprintf(fp, "  endloop\n");
    fprintf(fp, "endfacet\n");
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
                -1, 0, 0);
        maze_export_stl_triangle(fp, x-dx, y-dy, z-dz,
                x-dx, y+dy, z+dz,
                x-dx, y+dy, z-dz,
                -1, 0, 0);
    }
    if( (face_mask & (1<<1)) == 0) {
        /* right (+x) */
        maze_export_stl_triangle(fp, x+dx, y-dy, z-dz,
                x+dx, y+dy, z+dz,
                x+dx, y-dy, z+dz,
                1, 0, 0);
        maze_export_stl_triangle(fp, x+dx, y-dy, z-dz,
                x+dx, y+dy, z-dz,
                x+dx, y+dy, z+dz,
                1, 0, 0);
    }

    if( (face_mask & (1<<2)) == 0) {
        /* front (-y) */
        maze_export_stl_triangle(fp, x-dx, y-dy, z-dz,
                x+dx, y-dy, z+dz,
                x-dx, y-dy, z+dz,
                0, -1, 0);
        maze_export_stl_triangle(fp, x-dx, y-dy, z-dz,
                x+dx, y-dy, z-dz,
                x+dx, y-dy, z+dz,
                0, -1, 0);
    }
    if( (face_mask & (1<<3)) == 0) {
        /* back (+y) */
        maze_export_stl_triangle(fp, x-dx, y+dy, z-dz,
                x-dx, y+dy, z+dz,
                x+dx, y+dy, z+dz,
                0, 1, 0);
        maze_export_stl_triangle(fp, x-dx, y+dy, z-dz,
                x+dx, y+dy, z+dz,
                x+dx, y+dy, z-dz,
                0, 1, 0);
    }

    if( (face_mask & (1<<4)) == 0) {
        /* top (+z face) */
        maze_export_stl_triangle(fp, x-dx, y-dy, z+dz,
                x+dx, y-dy, z+dz,
                x+dx, y+dy, z+dz,
                0, 0, 1);
        maze_export_stl_triangle(fp, x-dx, y+dy, z+dz,
                x-dx, y-dy, z+dz,
                x+dx, y+dy, z+dz,
                0, 0, 1);
    }
    if( (face_mask & (1<<5)) == 0) {
        /* bottom (-z) */
        maze_export_stl_triangle(fp, x-dx, y-dy, z-dz,
                x+dx, y+dy, z-dz,
                x+dx, y-dy, z-dz,
                0, 0, -1);
        maze_export_stl_triangle(fp, x-dx, y+dy, z-dz,
                x+dx, y+dy, z-dz,
                x-dx, y-dy, z-dz,
                0, 0, -1);
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

                    /* export cubes */
                    maze_export_stl_cube(fp, x1, y1, z1, mask1, scale, scale, scale);
                    maze_export_stl_cube(fp, x2, y2, z2, mask2, scale, scale, scale);
                }
            }
        }
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
