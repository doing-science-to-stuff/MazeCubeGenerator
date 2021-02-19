#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maze.h"

typedef struct trig {
    double x[3], y[3], z[3];
} trig_t;

/* move individual triangle */
static void trig_move(trig_t *trig, double dx, double dy, double dz) {
    for(int i=0; i<3; ++i) {
        trig->x[i] += dx;
        trig->y[i] += dy;
        trig->z[i] += dz;
    }
}

/* rescale triangle */
static void trig_scale(trig_t *trig, double sx, double sy, double sz) {
    for(int i=0; i<3; ++i) {
        trig->x[i] *= sx;
        trig->y[i] *= sy;
        trig->z[i] *= sz;
    }

    if( sx*sy*sz < 0.0 ) {
        /* normal will be reversed */
        double temp;
        temp = trig->x[1]; trig->x[1] = trig->x[2]; trig->x[2] = temp;
        temp = trig->y[1]; trig->y[1] = trig->y[2]; trig->y[2] = temp;
        temp = trig->z[1]; trig->z[1] = trig->z[2]; trig->z[2] = temp;
    }
}

/* rotate triangle */
static void trig_rotate_axial(trig_t *trig, int axis, double rad) {
    /* rotate triangle rad radians around spcified axis */
    for(int i=0; i<3; ++i) {
        double x0, y0, x1, y1;

        /* select coordinates to rotate */
        switch(axis) {
            case 0:
                x0 = trig->y[i];
                y0 = trig->z[i];
                break;
            case 1:
                x0 = trig->x[i];
                y0 = trig->z[i];
                break;
            case 2:
            default:
                x0 = trig->x[i];
                y0 = trig->y[i];
                break;
        }

        /* rotate x0,y0 by r radians to get x1,y1 */
        /* see: https://en.wikipedia.org/wiki/Rotation_matrix#In_two_dimensions*/
        x1 = x0 * cos(rad) - y0 * sin(rad);
        y1 = x0 * sin(rad) + y0 * cos(rad);

        /* update appropriate coordinates */
        switch(axis) {
            case 0:
                trig->y[i] = x1;
                trig->z[i] = y1;
                break;
            case 1:
                trig->x[i] = x1;
                trig->z[i] = y1;
                break;
            case 2:
            default:
                trig->x[i] = x1;
                trig->y[i] = y1;
                break;
        }
    }
}

/* rotate triangle around point */
static void trig_rotate_axial_around(trig_t *trig, int axis, double rad, double cx, double cy, double cz) {
    trig_move(trig, -cx, -cy, -cz);
    trig_rotate_axial(trig, axis, rad);
    trig_move(trig, cx, cy, cz);
}

/* export single triangle as STL */
static void trig_get_normal(trig_t *trig,
                            double *nx, double *ny, double *nz) {
    /* get two edge vectors */
    double ux = trig->x[1]-trig->x[0];
    double uy = trig->y[1]-trig->y[0];
    double uz = trig->z[1]-trig->z[0];
    double vx = trig->x[2]-trig->x[0];
    double vy = trig->y[2]-trig->y[0];
    double vz = trig->z[2]-trig->z[0];

    /* compute normal for triangle defined by coordinates
     * see: https://mathworld.wolfram.com/CrossProduct.html
     * Equation 2 */
    *nx = uy*vz - uz*vy;
    *ny = uz*vx - ux*vz;
    *nz = ux*vy - uy*vx;
}

static void trig_export_stl(FILE *fp, trig_t *trig) {

    double nx, ny, nz;
    trig_get_normal(trig, &nx, &ny, &nz);
    fprintf(fp, "facet normal %g %g %g\n", nx, ny, nz);
    fprintf(fp, "  outer loop\n");
    fprintf(fp, "    vertex %g %g %g\n", trig->x[0], trig->y[0], trig->z[0]);
    fprintf(fp, "    vertex %g %g %g\n", trig->x[1], trig->y[1], trig->z[1]);
    fprintf(fp, "    vertex %g %g %g\n", trig->x[2], trig->y[2], trig->z[2]);
    fprintf(fp, "  endloop\n");
    fprintf(fp, "endfacet\n");
}


typedef struct trig_list {
    int num;
    int cap;
    trig_t *trig;
} trig_list_t;

/* initialize empty triangle list */
static int trig_list_init(trig_list_t *list) {
    memset(list,'\0',sizeof(*list));
    int initial_cap = 10;
    list->trig = calloc(initial_cap, sizeof(trig_t));
    list->cap = initial_cap;

    return 1;
}

/* free triangle list */
static void trig_list_free(trig_list_t *list) {
    free(list->trig); list->trig=NULL;
    memset(list,'\0',sizeof(*list));
}

/* add triangle to list */
static int trig_list_add(trig_list_t *list,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3, int rev) {

    /* reallocate list, if needed */
    if( list->num == list->cap ) {
        int new_cap = (list->cap*2) + 1;

        trig_t *new_buf = calloc(new_cap, sizeof(*list->trig));
        if( !new_buf ) { return 0; }

        memcpy(new_buf, list->trig, list->num*sizeof(*list->trig));
        free(list->trig); list->trig=NULL;

        list->trig = new_buf;
        list->cap = new_cap;
    }

    int pos = list->num;
    list->trig[pos].x[0] = x1;
    list->trig[pos].y[0] = y1;
    list->trig[pos].z[0] = z1;

    if( rev <= 0 ) {
        list->trig[pos].x[1] = x2;
        list->trig[pos].y[1] = y2;
        list->trig[pos].z[1] = z2;

        list->trig[pos].x[2] = x3;
        list->trig[pos].y[2] = y3;
        list->trig[pos].z[2] = z3;
    } else {
        list->trig[pos].x[1] = x3;
        list->trig[pos].y[1] = y3;
        list->trig[pos].z[1] = z3;

        list->trig[pos].x[2] = x2;
        list->trig[pos].y[2] = y2;
        list->trig[pos].z[2] = z2;
    }
    ++list->num;

    return 0;
}

/* copy one list into another */
static void trig_list_copy(trig_list_t *dst, trig_list_t *src) {
    for(int i=0; i<src->num; ++i) {
        trig_t *t = &src->trig[i];
        trig_list_add(dst, t->x[0], t->y[0], t->z[0],
            t->x[1], t->y[1], t->z[1],
            t->x[2], t->y[2], t->z[2], 0);
    }
}

/* move all triangles in list */
static void trig_list_move(trig_list_t *list, double dx, double dy, double dz) {
    for(int i=0; i<list->num; ++i) {
        trig_move(&list->trig[i], dx, dy, dz);
    }
}

/* scale all triangles in list */
static void trig_list_scale(trig_list_t *list, double sx, double sy, double sz) {
    for(int i=0; i<list->num; ++i) {
        trig_scale(&list->trig[i], sx, sy, sz);
    }
}

/* rotate triangles in list */
static void trig_list_rotate_axial(trig_list_t *list, int axis, double rad) {
    for(int i=0; i<list->num; ++i) {
        trig_rotate_axial(&list->trig[i], axis, rad);
    }
}

/* rotate triangles in list around specified point */
static void trig_list_rotate_axial_around(trig_list_t *list, int axis, double rad, double cx, double cy, double cz) {
    for(int i=0; i<list->num; ++i) {
        trig_rotate_axial_around(&list->trig[i], axis, rad, cx, cy, cz);
    }
}

/* export list of triangles as STL */
static void trig_list_export_stl(FILE *fp, trig_list_t *list) {
    for(int i=0; i<list->num; ++i) {
        trig_export_stl(fp, &list->trig[i]);
    }
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
static void maze_add_marker1(trig_list_t *list, maze_t *maze, int face, position_t pos, double radius, double scale, int dir) {
    int d1 = maze->faces[face].d1;
    int d2 = maze->faces[face].d2;

    /* get row and column for marker */
    int c = pos[d2];
    int r = pos[d1];

    int reverse = 0;

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
        for(int j=0; j<numSegsJ; ++j) {
            double thetaJ1 =  M_PI * j / numSegsJ;
            double thetaJ2 =  M_PI * (j+1) / numSegsJ;

            /* compute raw torus coordinates */
            double x11 = cos(thetaI1) * (1+radius*cos(thetaJ1))+r;
            double y11 = sin(thetaI1) * (1+radius*cos(thetaJ1))+c;
            double z11 = radius*sin(thetaJ1)+0.5;

            double x12 = cos(thetaI1) * (1+radius*cos(thetaJ2))+r;
            double y12 = sin(thetaI1) * (1+radius*cos(thetaJ2))+c;
            double z12 = radius*sin(thetaJ2)+0.5;

            double x21 = cos(thetaI2) * (1+radius*cos(thetaJ1))+r;
            double y21 = sin(thetaI2) * (1+radius*cos(thetaJ1))+c;
            double z21 = radius*sin(thetaJ1)+0.5;

            double x22 = cos(thetaI2) * (1+radius*cos(thetaJ2))+r;
            double y22 = sin(thetaI2) * (1+radius*cos(thetaJ2))+c;
            double z22 = radius*sin(thetaJ2)+0.5;

            /* TODO: When crossing a cell boundary, end cap points should be
             * interpolated to exactly match the boundary. */

            /* record first point along surface for end caps */
            if( j == 0 && (!cell1 || !cell2) ) {
                x10 = x11;
                y10 = y11;
                z10 = z11;
                x20 = x21;
                y20 = y21;
                z20 = z21;
            }

            /* output triangles */
            if( cell1 && cell2 ) {
                /* curved surface of marker */
                trig_list_add(list, x11, y11, z11,
                        x22, y22, z22,
                        x12, y12, z12, reverse);
                trig_list_add(list, x11, y11, z11,
                        x21, y21, z21,
                        x22, y22, z22, reverse);
            } else if( cell1 && !cell2 && j>0 ) {
                /* cap one end */
                trig_list_add(list, x10, y10, z10,
                                    x12, y12, z12,
                                    x11, y11, z11, reverse);
            } else if( !cell1 && cell2 && j>0 ) {
                /* cap other end */
                trig_list_add(list, x20, y20, z20,
                                    x21, y21, z21,
                                    x22, y22, z22, reverse);
            }
        }
    }
}

/* square marker */
static void maze_add_corner(trig_list_t *list, maze_t *maze, int r, int c, int dr, int dc, int face, double radius, double scale, int dir, int rCap, int cCap) {

    /* local triangle list */
    trig_list_t corner;
    trig_list_init(&corner);

    /* some faces need normals inverted */
    int reverse = 0;

    double x10 = 0.0, y10 = 0.0, z10 = 0.0;
    double x30 = 0.0, y30 = 0.0, z30 = 0.0;

    int numSegs = 32;
    for(int i=0; i<numSegs; ++i) {
        double thetaI1 = M_PI * i / numSegs;
        double thetaI2 = M_PI * (i+1) / numSegs;

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

        /* curved surface of marker */
        trig_list_add(&corner, x11, y11, z11,
                x22, y22, z22,
                x12, y12, z12, reverse);
        trig_list_add(&corner, x11, y11, z11,
                x21, y21, z21,
                x22, y22, z22, reverse);

        trig_list_add(&corner, x31, y31, z31,
                x32, y32, z32,
                x22, y22, z22, reverse);
        trig_list_add(&corner, x31, y31, z31,
                x22, y22, z22,
                x21, y21, z21, reverse);

        /* add end caps */
        if( rCap==0 ) {
            trig_list_add(&corner, x11, y11, z11,
                    x12, y12, z12,
                    x10, y10, z10, reverse);
        }

        if( cCap==0 ) {
            trig_list_add(&corner, x31, y31, z31,
                    x30, y30, z30,
                    x32, y32, z32, reverse);
        }
    }

    trig_list_scale(&corner, dr, dc, 1.0);
    trig_list_move(&corner, r, c, 0.0);

    trig_list_copy(list, &corner);
    trig_list_free(&corner);
}

static void maze_add_edge(trig_list_t *list, maze_t *maze, int r, int c, int face, double radius, double scale, int rotated, double rOffset, double cOffset) {

    /* local triangle list */
    trig_list_t edge;
    trig_list_init(&edge);

    /* some faces need normals inverted */
    int reverse = 0;

    #if 0
    if( dr < 0 )
        reverse ^= 1;
    #endif

    int numSegs = 32;
    for(int i=0; i<numSegs; ++i) {
        double thetaI1 = M_PI * i / numSegs;
        double thetaI2 = M_PI * (i+1) / numSegs;

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

        /* outer shell of marker */
        trig_list_add(&edge, x11, y11, z11,
                x12, y12, z12,
                x22, y22, z22,
                reverse);
        trig_list_add(&edge, x11, y11, z11,
                x22, y22, z22,
                x21, y21, z21,
                reverse);
    }

    if( rotated != 0 )
        trig_list_rotate_axial(&edge, 2, M_PI/2.0);
    trig_list_move(&edge, r+rOffset, c+cOffset, 0.0);

    trig_list_copy(list, &edge);
    trig_list_free(&edge);
}

static void maze_add_marker2(trig_list_t *list, maze_t *maze, int face, position_t pos, double radius, double scale, int dir) {
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
    maze_add_corner(list, maze, r, c, -1, -1, face, radius, scale, dir, lSide, tSide);
    maze_add_corner(list, maze, r, c, -1, 1, face, radius, scale, dir, lSide, bSide);
    maze_add_corner(list, maze, r, c, 1, -1, face, radius, scale, dir, rSide, tSide);
    maze_add_corner(list, maze, r, c, 1, 1, face, radius, scale, dir, rSide, bSide);

    /* add straight segments */
    if( lSide != 0 )
        maze_add_edge(list, maze, r, c, face, radius, scale, 0, -1, 0);
    if( rSide != 0 )
        maze_add_edge(list, maze, r, c, face, radius, scale, 0, 1, 0);
    if( tSide != 0 )
        maze_add_edge(list, maze, r, c, face, radius, scale, 1, 0, -1);
    if( bSide != 0 )
        maze_add_edge(list, maze, r, c, face, radius, scale, 1, 0, 1);
}

static void maze_add_cube(trig_list_t *list, double x, double y, double z, char face_mask, double scaleX, double scaleY, double scaleZ) {
    double dx = scaleX/2.0;
    double dy = scaleY/2.0;
    double dz = scaleZ/2.0;

    #if 0
    printf("%s:\n", __FUNCTION__);
    printf("cube centered at %g,%g,%g\n", x, y, z);
    printf("cube sizes are %g,%g,%g\n", scaleX, scaleY, scaleZ);
    printf("cube extents are %g,%g,%g to %g,%g,%g\n", x-dx, y-dy, z-dz, x+dx, y+dy, z+dz);
    #endif /* 0 */

    if( (face_mask & (1<<0)) == 0) {
        /* left (-x) */
        trig_list_add(list, x-dx, y-dy, z-dz,
                x-dx, y-dy, z+dz,
                x-dx, y+dy, z+dz, 0);
        trig_list_add(list, x-dx, y-dy, z-dz,
                x-dx, y+dy, z+dz,
                x-dx, y+dy, z-dz, 0);
    }
    if( (face_mask & (1<<1)) == 0) {
        /* right (+x) */
        trig_list_add(list, x+dx, y-dy, z-dz,
                x+dx, y+dy, z+dz,
                x+dx, y-dy, z+dz, 0);
        trig_list_add(list, x+dx, y-dy, z-dz,
                x+dx, y+dy, z-dz,
                x+dx, y+dy, z+dz, 0);
    }

    if( (face_mask & (1<<2)) == 0) {
        /* front (-y) */
        trig_list_add(list, x-dx, y-dy, z-dz,
                x+dx, y-dy, z+dz,
                x-dx, y-dy, z+dz, 0);
        trig_list_add(list, x-dx, y-dy, z-dz,
                x+dx, y-dy, z-dz,
                x+dx, y-dy, z+dz, 0);
    }
    if( (face_mask & (1<<3)) == 0) {
        /* back (+y) */
        trig_list_add(list, x-dx, y+dy, z-dz,
                x-dx, y+dy, z+dz,
                x+dx, y+dy, z+dz, 0);
        trig_list_add(list, x-dx, y+dy, z-dz,
                x+dx, y+dy, z+dz,
                x+dx, y+dy, z-dz, 0);
    }

    if( (face_mask & (1<<4)) == 0) {
        /* bottom (-z) */
        trig_list_add(list, x-dx, y-dy, z-dz,
                x+dx, y+dy, z-dz,
                x+dx, y-dy, z-dz, 0);
        trig_list_add(list, x-dx, y+dy, z-dz,
                x+dx, y+dy, z-dz,
                x-dx, y-dy, z-dz, 0);
    }
    if( (face_mask & (1<<5)) == 0) {
        /* top (+z face) */
        trig_list_add(list, x-dx, y-dy, z+dz,
                x+dx, y-dy, z+dz,
                x+dx, y+dy, z+dz, 0);
        trig_list_add(list, x-dx, y+dy, z+dz,
                x-dx, y-dy, z+dz,
                x+dx, y+dy, z+dz, 0);
    }
}


int maze_add_maze_face(maze_t *maze, int face, trig_list_t *list) {

    double scale = 1.0;

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

                /* export cube */
                maze_add_cube(list, x1, y1, z1, mask1, scale, scale, scale);
            }
        }
    }

    /* add end markers */
    double markerRadius = 0.1;
    maze_add_marker2(list, maze, face, maze->startPos, markerRadius, scale, 1);
    maze_add_marker1(list, maze, face, maze->endPos, markerRadius, scale, 1);

    return 1;
}


int maze_add_maze(maze_t *maze, trig_list_t *list) {

    double scale = 1.0;
    /* for each face */
    for(int face=0; face<maze->numFaces; ++face) {

        trig_list_t faceTrigs1, faceTrigs2;
        trig_list_init(&faceTrigs1);
        trig_list_init(&faceTrigs2);
        maze_add_maze_face(maze, face, &faceTrigs1);
        maze_add_maze_face(maze, face, &faceTrigs2);

        /* translate face */
        if( face == 0 ) {
            trig_list_scale(&faceTrigs2, 1.0, 1.0, -1.0);
            trig_list_move(&faceTrigs2, 0.0, 0.0, 0.5);
            trig_list_move(&faceTrigs2, 0.0, 0.0, -0.5);
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
        
        /* add face to maze list */
        trig_list_copy(list, &faceTrigs1);
        trig_list_copy(list, &faceTrigs2);

        trig_list_free(&faceTrigs1);
        trig_list_free(&faceTrigs2);
    }

    #if 0
    /* add frame corners */
    int xEnd = maze->dimensions[0]-1;
    int yEnd = maze->dimensions[1]-1;
    int zEnd = maze->dimensions[2]-1;
    maze_add_cube(list,    0,    0,    0, (1<<1)|(1<<3)|(1<<5), scale, scale, scale);
    maze_add_cube(list, xEnd,    0,    0, (1<<0)|(1<<3)|(1<<5), scale, scale, scale);
    maze_add_cube(list,    0, yEnd,    0, (1<<1)|(1<<2)|(1<<5), scale, scale, scale);
    maze_add_cube(list, xEnd, yEnd,    0, (1<<0)|(1<<2)|(1<<5), scale, scale, scale);
    maze_add_cube(list,    0,    0, zEnd, (1<<1)|(1<<3)|(1<<4), scale, scale, scale);
    maze_add_cube(list, xEnd,    0, zEnd, (1<<0)|(1<<3)|(1<<4), scale, scale, scale);
    maze_add_cube(list,    0, yEnd, zEnd, (1<<1)|(1<<2)|(1<<4), scale, scale, scale);
    maze_add_cube(list, xEnd, yEnd, zEnd, (1<<0)|(1<<2)|(1<<4), scale, scale, scale);

    /* add frame edges */
    /* x */
    maze_add_cube(list, xEnd/2,    0,    0, (1<<0)|(1<<1), (xEnd-1)*scale, scale, scale);
    maze_add_cube(list, xEnd/2, yEnd,    0, (1<<0)|(1<<1), (xEnd-1)*scale, scale, scale);
    maze_add_cube(list, xEnd/2,    0, zEnd, (1<<0)|(1<<1), (xEnd-1)*scale, scale, scale);
    maze_add_cube(list, xEnd/2, yEnd, zEnd, (1<<0)|(1<<1), (xEnd-1)*scale, scale, scale);
    /* y */
    maze_add_cube(list,    0, yEnd/2,    0, (1<<2)|(1<<3), scale, (yEnd-1)*scale, scale);
    maze_add_cube(list, xEnd, yEnd/2,    0, (1<<2)|(1<<3), scale, (yEnd-1)*scale, scale);
    maze_add_cube(list,    0, yEnd/2, zEnd, (1<<2)|(1<<3), scale, (yEnd-1)*scale, scale);
    maze_add_cube(list, xEnd, yEnd/2, zEnd, (1<<2)|(1<<3), scale, (yEnd-1)*scale, scale);
    /* z */
    maze_add_cube(list,    0,    0, zEnd/2, (1<<4)|(1<<5), scale, scale, (zEnd-1)*scale);
    maze_add_cube(list, xEnd,    0, zEnd/2, (1<<4)|(1<<5), scale, scale, (zEnd-1)*scale);
    maze_add_cube(list,    0, yEnd, zEnd/2, (1<<4)|(1<<5), scale, scale, (zEnd-1)*scale);
    maze_add_cube(list, xEnd, yEnd, zEnd/2, (1<<4)|(1<<5), scale, scale, (zEnd-1)*scale);
    #endif // 0?

    /* add slider */
    double startX = maze->startPos[0];
    double startY = maze->startPos[1];
    double startZ = maze->startPos[2];
    double sizeX = 2*maze->dimensions[0]-0.5;
    double sizeY = 2*maze->dimensions[1]-0.5;
    double sizeZ = 2*maze->dimensions[2]-0.5;
    double scale2 = scale*0.95;
    maze_add_cube(list, startX, startY, startZ, 0, sizeX*scale, scale2, scale2);
    maze_add_cube(list, startX, startY, startZ, 0, scale, sizeY*scale2, scale2);
    maze_add_cube(list, startX, startY, startZ, 0, scale, scale2, sizeZ*scale2);

    trig_list_move(list, 0.5, 0.5, 0.5);

    return 0;
}


static int maze_add_flat_border(trig_list_t *list, double xOffset, double yOffset, double xSize, double ySize, double scale) {

    #if 0
    printf("%s\n", __FUNCTION__);
    printf("x,y: %g,%g; w,h: %g,%g\n", xOffset, yOffset, xSize, ySize);
    #endif /* 0 */

    double xi[4], yi[4], xo[4], yo[4];

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

    double z0 = -0.5*scale;
    double z1 = 0.5*scale;

    /* adjust scaling */
    for(int i=0; i<4; ++i) {
        xo[i] *= scale;
        yo[i] *= scale;
        xi[i] *= scale;
        yi[i] *= scale;
    }

    for(int i=0; i<4; ++i) {
        int j = (i+1)%4;

        /* top flat face */
        trig_list_add(list, xo[i], yo[i], z1,
                                     xo[j], yo[j], z1,
                                     xi[i], yi[i], z1, 0);
        trig_list_add(list, xi[i], yi[i], z1,
                                     xo[j], yo[j], z1,
                                     xi[j], yi[j], z1, 0);

        /* outer sloped face */
        trig_list_add(list, xo[i], yo[i], z1,
                                     xi[i], yi[i], z0,
                                     xi[j], yi[j], z0, 0);
        trig_list_add(list, xo[i], yo[i], z1,
                                     xi[j], yi[j], z0,
                                     xo[j], yo[j], z1, 0);

        /* inner vertical face */
        trig_list_add(list, xi[i], yi[i], z1,
                                     xi[j], yi[j], z1,
                                     xi[j], yi[j], z0, 0);
        trig_list_add(list, xi[i], yi[i], z1,
                                     xi[j], yi[j], z0,
                                     xi[i], yi[i], z0, 0);
    }

    return 0;
}


static int maze_add_flat_slider(trig_list_t *list, double xOffset, double yOffset, double xSize, double ySize, double scale) {

    double size = 0.9*scale;
    double size2 = size/2.0;
    double xPos = xOffset*scale;
    double yPos = yOffset*scale;
    double zPos = (size-scale)/2.0;

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


int maze_add_maze_flat(maze_t *maze, trig_list_t *list) {
    double scale = 1.0;
    /* write faces */
    /* for each face */
    for(int face=0; face<maze->numFaces; ++face) {
        int rows = maze->faces[face].rows;
        int cols = maze->faces[face].cols;

        trig_list_t faceTrigs1, faceTrigs2;
        trig_list_init(&faceTrigs1);
        trig_list_init(&faceTrigs2);
        maze_add_maze_face(maze, face, &faceTrigs1);
        maze_add_maze_face(maze, face, &faceTrigs2);

        /* add perimeters */
        maze_add_flat_border(&faceTrigs1, -0.5, -0.5, cols, rows, scale);
        maze_add_flat_border(&faceTrigs2, -0.5, -0.5, cols, rows, scale);

        /* compute base offsets for face */
        double xBase1 = -cols-1.0;
        double yBase1 = 0.0;
        for(int f2 = 0; f2<face; ++f2) {
            yBase1 += maze->faces[f2].rows+3;
        }
        double xBase2 = cols+1;
        double yBase2 = yBase1;

        /* translate face */
        trig_list_scale(&faceTrigs2, -1.0, 1.0, 1.0);
        trig_list_move(&faceTrigs1, xBase1, yBase1, 0.0);
        trig_list_move(&faceTrigs2, xBase2, yBase2, 0.0);
        
        /* add face to maze list */
        trig_list_copy(list, &faceTrigs1);
        trig_list_copy(list, &faceTrigs2);

        trig_list_free(&faceTrigs1);
        trig_list_free(&faceTrigs2);
    }

    /* write slider pieces */
    int xSize, ySize, zSize;
    xSize = maze->dimensions[0];
    ySize = maze->dimensions[1];
    zSize = maze->dimensions[2];
    maze_add_flat_slider(list, 0.0, ySize+1.0, xSize+1, -(ySize+1), scale);
    maze_add_flat_slider(list, 0.0, ySize+zSize+4.0, zSize+1, ySize+1, scale);

    trig_list_move(list, 0.5, 0.5, 0.5);

    return 0;
}

int maze_export_stl_flat(maze_t *maze, char *filename) {

    if( maze->numDimensions != 3 ) {
        fprintf(stderr,"%s: STL export is only supported for 3D mazes.\n", __FUNCTION__);
        return -1;
    }

    trig_list_t trigs;
    trig_list_init(&trigs);
    maze_add_maze_flat(maze, &trigs);

    /* open file */
    FILE *fp = fopen(filename,"w");
    if( fp == NULL ) {
        perror("fopen");
        return -1;
    }

    /* start solid */
    fprintf(fp,"solid MazeCubeFaces\n");

    /* write triangles to output file */
    trig_list_export_stl(fp, &trigs);

    /* free triangle list */
    trig_list_free(&trigs);

    /* close solid */
    fprintf(fp,"endsolid MazeCubeFaces\n");

    /* close file */
    fclose(fp); fp=NULL;

    return 0;
}


int maze_add_solution(maze_t *maze, trig_list_t *list) {
    double scale = 1.0;
    for(int i=0; i<maze->solution.posListNum; ++i) {
        double x, y, z;
        x = maze->solution.positions[i][0];
        y = maze->solution.positions[i][1];
        z = maze->solution.positions[i][2];
        maze_add_cube(list, x, y, z, 0, scale, scale, scale);
    }

    return 1;
}

int maze_export_stl_solution(maze_t *maze, char *filename) {
    if( maze->numDimensions != 3 ) {
        fprintf(stderr,"%s: STL export is only supported for 3D mazes.\n", __FUNCTION__);
        return -1;
    }

    trig_list_t trigs;
    trig_list_init(&trigs);
    maze_add_solution(maze, &trigs);

    /* open file */
    FILE *fp = fopen(filename,"w");
    if( fp == NULL ) {
        perror("fopen");
        return -1;
    }

    /* open solid */
    fprintf(fp,"solid MazeCubeSolution\n");

    /* write triangles to output file */
    trig_list_export_stl(fp, &trigs);

    /* free triangle list */
    trig_list_free(&trigs);

    /* open solid */
    fprintf(fp,"endsolid MazeCubeSolution\n");

    fclose(fp); fp=NULL;

    return 1;
}

int maze_export_stl(maze_t *maze, char *filename) {
    if( maze->numDimensions != 3 ) {
        fprintf(stderr,"%s: STL export is only supported for 3D mazes.\n", __FUNCTION__);
        return -1;
    }

    trig_list_t trigs;
    trig_list_init(&trigs);
    maze_add_maze(maze, &trigs);

    /* open file */
    FILE *fp = fopen(filename,"w");
    if( fp == NULL ) {
        perror("fopen");
        return -1;
    }

    /* start solid */
    fprintf(fp,"solid MazeCube\n");

    /* write triangles to output file */
    trig_list_export_stl(fp, &trigs);

    /* free triangle list */
    trig_list_free(&trigs);

    /* close solid */
    fprintf(fp,"endsolid MazeCube\n");

    /* close file */
    fclose(fp); fp=NULL;

    return 0;
}
