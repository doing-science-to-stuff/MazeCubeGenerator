/*
 * trig-list.h
 * MazeCubeGen: maze cube generator
 *
 * Copyright (c) 2026 Bryan Franklin. All rights reserved.
 */
#ifndef TRIG_LIST_H
#define TRIG_LIST_H

const double epsilon = 1e-8;

typedef enum face_group {
    FACE_OPEN = -2,
    FACE_NONE = -1,
    FACE_0A,
    FACE_0B,
    FACE_1A,
    FACE_1B,
    FACE_2A,
    FACE_2B,
    FACE_MARKER_1,
    FACE_MARKER_2
} face_group_t;


typedef struct trig {
    double x[3], y[3], z[3];    /* vertex coordinates */
    double nx[3], ny[3], nz[3];   /* vertex normals */
    face_group_t groupId;
} trig_t;


/* adjust lengths of normals to be unit length */
static void trig_unitize_normals(trig_t *trig) {
    for(int i=0; i<3; ++i) {
        /* get lengths of the normal vectors */
        double len = sqrt(pow(trig->nx[i],2.0)
                        + pow(trig->ny[i],2.0)
                        + pow(trig->nz[i],2.0));

        /* "normalize" the normals */
        if( fabs(len) > epsilon ) {
            double invLen = 1.0/len;
            trig->nx[i] *= invLen;
            trig->ny[i] *= invLen;
            trig->nz[i] *= invLen;
        }
    }
}


/* set the normal vector for all vertices */
static void trig_set_normals(trig_t *trig,
                             double x1, double y1, double z1,
                             double x2, double y2, double z2,
                             double x3, double y3, double z3) {
    /* update triangle */
    trig->nx[0] = x1; trig->nx[1] = x2; trig->nx[2] = x3;
    trig->ny[0] = y1; trig->ny[1] = y2; trig->ny[2] = y3;
    trig->nz[0] = z1; trig->nz[1] = z2; trig->nz[2] = z3;

    /* "normalize" the normals */
    trig_unitize_normals(trig);
}


/* compute the normal of flat triangle using the cross product of two edge vectors */
static void trig_get_normal(trig_t *trig) {
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
    double nx = uy*vz - uz*vy;
    double ny = uz*vx - ux*vz;
    double nz = ux*vy - uy*vx;

    /* round to epsilon precision */
    nx = epsilon*round(nx/epsilon);
    ny = epsilon*round(ny/epsilon);
    nz = epsilon*round(nz/epsilon);

    trig_set_normals(trig,
                     nx, ny, nz,
                     nx, ny, nz,
                     nx, ny, nz);
}


/* initialize a trig */
static void trig_init(trig_t *trig) {
    if( !trig ) return;
    memset(trig, '\0', sizeof(*trig));
    trig->groupId = FACE_NONE;
}


/* fill in a triangle */
static void trig_fill(trig_t *trig,
                      double x1, double y1, double z1,
                      double x2, double y2, double z2,
                      double x3, double y3, double z3) {
    if( !trig ) return;
    trig_init(trig);
    trig->x[0] = x1; trig->x[1] = x2; trig->x[2] = x3;
    trig->y[0] = y1; trig->y[1] = y2; trig->y[2] = y3;
    trig->z[0] = z1; trig->z[1] = z2; trig->z[2] = z3;
    trig_get_normal(trig);
}


/* set a grouping id to aid in color assignment */
static void trig_set_group(trig_t *trig, face_group_t id) {
    if( !trig ) return;
    trig->groupId = id;
}


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
    #if 0
    printf("Scaling trig:\t<%g, %g, %g>,\t<%g, %g, %g>,\t<%g, %g, %g> by %g, %g, %g\n",
        trig->x[0], trig->y[0], trig->z[0],
        trig->x[1], trig->y[1], trig->z[1],
        trig->x[2], trig->y[2], trig->z[2],
        sx, sy, sz);
    #if 0
    printf("normals: \t<%g, %g, %g>,\t<%g, %g, %g>,\t<%g, %g, %g>\n",
            trig->nx[0], trig->ny[0], trig->nz[0],
            trig->nx[1], trig->ny[1], trig->nz[1],
            trig->nx[2], trig->ny[2], trig->nz[2]);
    #endif /* 0 */
    #endif /* 0 */

    for(int i=0; i<3; ++i) {
        trig->x[i] *= sx;
        trig->y[i] *= sy;
        trig->z[i] *= sz;
        /* see: https://paroj.github.io/gltut/Illumination/Tut09%20Normal%20Transformation.html */
        if( fabs(sx) > epsilon ) trig->nx[i] /= sx;
        if( fabs(sy) > epsilon ) trig->ny[i] /= sy;
        if( fabs(sz) > epsilon ) trig->nz[i] /= sz;
    }

    if( sx*sy*sz < 0.0 ) {
        #if 0
        printf("swapping vertices due to normal.\n");
        #endif /* 0 */
        /* normal will be reversed, so vertex order needs to reverse as well */
        double temp;
        temp = trig->x[1]; trig->x[1] = trig->x[2]; trig->x[2] = temp;
        temp = trig->y[1]; trig->y[1] = trig->y[2]; trig->y[2] = temp;
        temp = trig->z[1]; trig->z[1] = trig->z[2]; trig->z[2] = temp;
        /* normals need to be swapped along with vertices */
        temp = trig->nx[1]; trig->nx[1] = trig->nx[2]; trig->nx[2] = temp;
        temp = trig->ny[1]; trig->ny[1] = trig->ny[2]; trig->ny[2] = temp;
        temp = trig->nz[1]; trig->nz[1] = trig->nz[2]; trig->nz[2] = temp;
    }

    /* "normalize" the normals */
    trig_unitize_normals(trig);

    #if 0
    printf("\t->\t<%g, %g, %g>,\t<%g, %g, %g>,\t<%g, %g, %g>\n",
        trig->x[0], trig->y[0], trig->z[0],
        trig->x[1], trig->y[1], trig->z[1],
        trig->x[2], trig->y[2], trig->z[2]);
    #if 0
    printf("normals: \t<%g, %g, %g>,\t<%g, %g, %g>,\t<%g, %g, %g>\n",
            trig->nx[0], trig->ny[0], trig->nz[0],
            trig->nx[1], trig->ny[1], trig->nz[1],
            trig->nx[2], trig->ny[2], trig->nz[2]);
    #endif /* 0 */
    #endif /* 0 */
}


static void trig_set_minimum(trig_t *trig, double min, int dim) {
    for(int i=0; i<3; ++i) {
        switch(dim) {
        case 0:
            if( trig->x[i] < min ) trig->x[i] = min;
            break;
        case 1:
            if( trig->y[i] < min ) trig->y[i] = min;
            break;
        case 2:
            if( trig->z[i] < min ) trig->z[i] = min;
            break;
        }
    }
    trig_get_normal(trig);
}


/* rotate triangle */
static void trig_rotate_axial(trig_t *trig, int axis, double rad) {
    /* rotate triangle rad radians around spcified axis */
    for(int i=0; i<3; ++i) {
        double x0, y0, x1, y1;
        double nx0, ny0, nx1, ny1;

        /* select coordinates to rotate */
        switch(axis) {
            case 0:
                x0 = trig->y[i];
                y0 = trig->z[i];
                nx0 = trig->ny[i];
                ny0 = trig->nz[i];
                break;
            case 1:
                x0 = trig->x[i];
                y0 = trig->z[i];
                nx0 = trig->nx[i];
                ny0 = trig->nz[i];
                break;
            case 2:
            default:
                x0 = trig->x[i];
                y0 = trig->y[i];
                nx0 = trig->nx[i];
                ny0 = trig->ny[i];
                break;
        }

        /* rotate x0,y0 and nx0,ny0 by r radians to get x1,y1 and nx1,ny1 */
        /* see: https://en.wikipedia.org/wiki/Rotation_matrix#In_two_dimensions*/
        x1 = x0 * cos(rad) - y0 * sin(rad);
        y1 = x0 * sin(rad) + y0 * cos(rad);
        nx1 = nx0 * cos(rad) - ny0 * sin(rad);
        ny1 = nx0 * sin(rad) + ny0 * cos(rad);

        /* update appropriate coordinates */
        switch(axis) {
            case 0:
                trig->y[i] = x1;
                trig->z[i] = y1;
                trig->ny[i] = nx1;
                trig->nz[i] = ny1;
                break;
            case 1:
                trig->x[i] = x1;
                trig->z[i] = y1;
                trig->nx[i] = nx1;
                trig->nz[i] = ny1;
                break;
            case 2:
            default:
                trig->x[i] = x1;
                trig->y[i] = y1;
                trig->nx[i] = nx1;
                trig->ny[i] = ny1;
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
static void trig_export_stl(FILE *fp, trig_t *trig) {

    /* get normal for triangle */
    trig_get_normal(trig);

    /* round all values to remove noise */
    for(int i=0; i<3; ++i) {
        trig->x[i] = epsilon*round(trig->x[i]/epsilon);
        trig->y[i] = epsilon*round(trig->y[i]/epsilon);
        trig->z[i] = epsilon*round(trig->z[i]/epsilon);
        trig->nx[i] = epsilon*round(trig->nx[i]/epsilon);
        trig->ny[i] = epsilon*round(trig->ny[i]/epsilon);
        trig->nz[i] = epsilon*round(trig->nz[i]/epsilon);
    }

    /* output triangle to fp */
    /* Note: since STL only has one normal per facet
             and trig_get_normal sets all three to be the same,
             just use first vertex's normal. */
    fprintf(fp, "facet normal %g %g %g\n", trig->nx[0], trig->ny[0], trig->nz[0]);
    fprintf(fp, "  outer loop\n");
    fprintf(fp, "    vertex %g %g %g\n", trig->x[0], trig->y[0], trig->z[0]);
    fprintf(fp, "    vertex %g %g %g\n", trig->x[1], trig->y[1], trig->z[1]);
    fprintf(fp, "    vertex %g %g %g\n", trig->x[2], trig->y[2], trig->z[2]);
    fprintf(fp, "  endloop\n");
    fprintf(fp, "endfacet\n");
}


typedef struct trig_list {
    int num;    /* number of triangles in list */
    int cap;    /* allocated capacity of list */
    trig_t *trig;   /* list buffer */
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


/* reallocate list, if needed */
static void trig_list_resize(trig_list_t *list) {
    if( !list ) { return; }
    
    if( list->num == list->cap ) {
        int new_cap = (list->cap*2) + 1;

        trig_t *new_buf = calloc(new_cap, sizeof(*list->trig));
        if( !new_buf ) { return; }

        memcpy(new_buf, list->trig, list->num*sizeof(*list->trig));
        free(list->trig); list->trig=NULL;

        list->trig = new_buf;
        list->cap = new_cap;
    }
}


/* add triangle to list */
static int trig_list_add(trig_list_t *list,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double x3, double y3, double z3) {

    #if 0
    printf("Adding trig: <%g, %g, %g>, <%g, %g, %g>, <%g, %g, %g>\n",
        x1, y1, z1, x2, y2, z2, x3, y3, z3);
    #endif /* 0 */

    /* reallocate list, if needed */
    trig_list_resize(list);

    int pos = list->num;
    trig_init(&list->trig[pos]);
    list->trig[pos].x[0] = x1;
    list->trig[pos].y[0] = y1;
    list->trig[pos].z[0] = z1;

    list->trig[pos].x[1] = x2;
    list->trig[pos].y[1] = y2;
    list->trig[pos].z[1] = z2;

    list->trig[pos].x[2] = x3;
    list->trig[pos].y[2] = y3;
    list->trig[pos].z[2] = z3;

    trig_get_normal(&list->trig[pos]);

    ++list->num;

    return 0;
}

static int trig_list_append(trig_list_t *list, trig_t *t) {

    /* reallocate list, if needed */
    trig_list_resize(list);

    /* copy trig into list */
    trig_t *new_pos = &list->trig[list->num];
    memcpy(new_pos, t, sizeof(*t));
    ++list->num;

    return 1;
}

/* copy one list onto end of another list */
static void trig_list_concatenate(trig_list_t *dst, trig_list_t *src) {
    for(int i=0; i<src->num; ++i) {
        trig_list_append(dst, &src->trig[i]);
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


/* set minimum value in specified dimension for all points */
static void trig_list_set_minimum(trig_list_t *list, double min, int dim) {
    for(int i=0; i<list->num; ++i) {
        trig_set_minimum(&list->trig[i], min, dim);
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


/* set group id for all triangles in a list */
static void trig_list_set_groupid(trig_list_t *list, face_group_t id) {
    if( !list ) return;
    for(int i=0; i<list->num; ++i) {
        trig_set_group(&list->trig[i], id);
    }
}


/* replace group ids for all triangles with the target group id */
static void trig_list_replace_groupid(trig_list_t *list, face_group_t id, face_group_t target_id) {
    if( !list ) return;
    for(int i=0; i<list->num; ++i) {
        if( list->trig[i].groupId == target_id ) {
            trig_set_group(&list->trig[i], id);
        }
    }
}


static int trig_list_write_stl(trig_list_t *trigs, char *filename, char *name) {

    /* remove destination if list is empty. */
    if( trigs-> num <= 0 ) {
        //printf("  skipping empty '%s'.\n", name);
        return unlink(filename);
    }
    printf("  writing '%s' to '%s'.\n", name, filename);

    /* open file */
    FILE *fp = fopen(filename,"w");
    if( fp == NULL ) {
        perror("fopen");
        return -1;
    }

    /* start solid */
    fprintf(fp,"solid %s\n", name);

    /* write triangles to output file */
    trig_list_export_stl(fp, trigs);

    /* close solid */
    fprintf(fp,"endsolid %s\n", name);

    /* close file */
    fclose(fp); fp=NULL;

    return 0;
}


int is_same_point(double x1, double y1, double z1,
                  double x2, double y2, double z2) {
    if( fabs(x2-x1) < epsilon 
        && fabs(y2-y1) < epsilon
        && fabs(z2-z1) < epsilon)
            return 1;
    return 0;
}


int trig_has_vertex(trig_t *trig, double x, double y, double z) {
    for(int i=0; i<3; ++i) {
        if( is_same_point(trig->x[i], trig->y[i], trig->z[i], x, y, z) )
            return 1;
    }

    return 0;
}

int is_same_trig(trig_t *trig1, trig_t *trig2) {
    for(int i=0; i<3; ++i) {
        int found = 0;
        for(int k=0; k<3; ++k) {
            if( is_same_point(trig1->x[i], trig1->y[i], trig1->z[i],
                              trig2->x[k], trig2->y[k], trig2->z[k]) )
                found = 1;
        }
        if( !found )
            return 0;
    }

    return 1;   // all points in trig1 found in trig2
}


int trig_shares_edge_num(trig_t *trig1, int edge1, trig_t *trig2) {
    unsigned int i=edge1%3;
    unsigned int j=(i+1)%3;
    if( trig_has_vertex(trig2, trig1->x[i], trig1->y[i], trig1->z[i])
            && trig_has_vertex(trig2, trig1->x[j], trig1->y[j], trig1->z[j]) )
            return 1;

    return 0;
}


int trig_shares_edge(trig_t *trig1, trig_t *trig2) {
    for(int i=0; i<3; ++i) {
        if( trig_shares_edge_num(trig1, i, trig2) )
            return 1;
    }

    return 0;
}


int trig_has_open_edge(trig_t *trig, trig_list_t *trigs){
    // for each edge in trig
    for(int j=0; j<3; ++j) {
        int found=0;
        for(int i=0; i<trigs->num; ++i) {
            // check is trigs[i] shares edge j with another trig
            if( trig_shares_edge_num(trig, j, &trigs->trig[i])
                    && !is_same_trig(trig, &trigs->trig[i]) )
                found = 1;
        }
        if( !found )
            return 1;
    }

    return 0;
}


int find_open_edges(trig_list_t *trigs) {
    // for each trig
    int count=0;
    for(int i=0; i<trigs->num; ++i) {
        if( trig_has_open_edge(&trigs->trig[i], trigs) ) {
            trigs->trig[i].groupId = FACE_OPEN;
            ++count;
        }
    }

    return count;
}

#endif // TRIG_LIST_H
