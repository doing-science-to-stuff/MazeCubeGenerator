/*
 * maze-cube.c
 *
 * Copyright (c) 2019-2020 Bryan Franklin. All rights reserved.
 */
#include <stdio.h>
#include "../scene.h"

#include "../../PuzzleMaze/maze.c"

static maze_t maze;
/* there are 90 moves in the original puzzle each direction */
#if 0
static int framesPerMove = 4;  /* 30s@24fps */
#else
static int framesPerMove = 10;  /* 30s@60fps */
#endif
static double edge_size = 30;
static int framesPerSpin = 240; /* 4s@60fps */
static double marker_radius = 0.05;

/* scene_frames is optional, but gives the total number of frames to render
 * for an animated scene. */
int scene_frames(int dimensions, char *config) {
    maze_load(&maze,config);
    int numFrames = (maze.solution.posListNum-1)*framesPerMove + 2*framesPerSpin;
    return numFrames;
}

static int add_mirror(scene *scn, int dimensions, int which, double mirror_dist) {

    vectNd offset, normal;
    vectNd_calloc(&offset,dimensions);
    vectNd_calloc(&normal,dimensions);

    /* create a mirrored hplane */
    object *mirror=NULL;
    scene_alloc_object(scn,dimensions,&mirror,"hplane");

    mirror->red = 0.1;
    mirror->green = 0.1;
    mirror->blue = 0.1;
    mirror->red_r = 0.95;
    mirror->green_r = 0.95;
    mirror->blue_r = 0.95;

    /* set normal */
    vectNd_reset(&offset);
    vectNd_set(&normal, which, -mirror_dist);
    object_add_dir(mirror, &normal);

    /* set position */
    vectNd_reset(&offset);
    vectNd_set(&offset,which,mirror_dist);
    object_add_pos(mirror,&offset);

    return 0;
}

static void set_face_color(object *obj, int face) {
    switch(face) {
        case 0:
            obj->red = 1.0;
            obj->green = 1.0;
            obj->blue = 1.0;
            break;
        case 1:
            obj->red = 0.0;
            obj->green = 1.0;
            obj->blue = 0.0;
            break;
        case 2:
            obj->red = 0.0;
            obj->green = 0.0;
            obj->blue = 1.0;
            break;
        case 3:
            obj->red = 1.0;
            obj->green = 1.0;
            obj->blue = 0.0;
            break;
        case 4:
            obj->red = 0.0;
            obj->green = 1.0;
            obj->blue = 1.0;
            break;
        case 5:
            obj->red = 0.5;
            obj->green = 1.0;
            obj->blue = 0.0;
            break;
        case 6:
            obj->red = 0.0;
            obj->green = 0.5;
            obj->blue = 1.0;
            break;
        case 7:
            obj->red = 1.0;
            obj->green = 0.0;
            obj->blue = 1.0;
            break;
        default:
            obj->red = 1.0;
            obj->green = 0.0;
            obj->blue = 0.0;
            break;
    }
}

static void add_maze_marker_segment(object *marker, int face,
    double x1, double y1, double x2, double y2,
    double scale, double radius) {

    vectNd pos1, pos2;
    vectNd_calloc(&pos1, maze.numDimensions);
    vectNd_calloc(&pos2, maze.numDimensions);

    int d1 = maze.faces[face].d1;
    int d2 = maze.faces[face].d2;
    vectNd_reset(&pos1);
    vectNd_reset(&pos2);
    vectNd_set(&pos1, d1, x1*scale);
    vectNd_set(&pos1, d2, y1*scale);
    vectNd_set(&pos2, d1, x2*scale);
    vectNd_set(&pos2, d2, y2*scale);

    /* add an object from pos1 to pos2 */
    object *cyl = object_alloc(maze.numDimensions, "hcylinder", "marker piece");
    object_add_obj(marker, cyl);
    object_add_pos(cyl, &pos1);
    object_add_pos(cyl, &pos2);
    vectNd extra;
    vectNd_calloc(&extra, maze.numDimensions);
    for(int i=0; i<maze.numDimensions-3; ++i) {
        #if 1
        vectNd_copy(&extra,&pos1);
        for(int j=0; j<maze.numDimensions; ++j) {
            if( j == d1 || j == d2 ) continue;
            vectNd_set(&extra, j, scale*radius);
        }
        object_add_pos(cyl, &extra);
        #else
        object_add_pos(cyl, &pos2);
        #endif // 0
    }
    vectNd_free(&extra);
    object_add_size(cyl, radius*scale);
    object_add_flag(cyl, 0);
    set_face_color(cyl, face);

    /* add caps at end */
    vectNd norm;
    vectNd_calloc(&norm, maze.numDimensions);
    vectNd_sub(&pos1, &pos2, &norm);
    object *end = object_alloc(maze.numDimensions, "hdisk", "marker piece end 1");
    object_add_obj(marker, end);
    object_add_pos(end, &pos1);
    object_add_dir(end, &norm);
    object_add_size(end, radius*scale);
    set_face_color(end, face);
    end = object_alloc(maze.numDimensions, "hdisk", "marker piece end 2");
    object_add_obj(marker, end);
    object_add_pos(end, &pos2);
    object_add_dir(end, &norm);
    object_add_size(end, radius*scale);
    set_face_color(end, face);

}

static object *make_maze_marker1(int face, double scale, int r, int c) {

    int numSegs = 32;
    double radius = marker_radius;

    object *marker = object_alloc(maze.numDimensions, "cluster", "maze marker 1");
    object_add_flag(marker,4);

    int d1 = maze.faces[face].d1;
    int d2 = maze.faces[face].d2;
    #if 1
    vectNd pos1, pos2, posm1;
    vectNd_calloc(&pos1, maze.numDimensions);
    vectNd_calloc(&pos2, maze.numDimensions);
    vectNd_calloc(&posm1, maze.numDimensions);
    for(int i=0; i<numSegs; ++i) {
        double thetam1 =  2.0 * M_PI * (i-1) / numSegs;
        double xm1 = cos(thetam1);
        double ym1 = sin(thetam1);
        double theta1 =  2.0 * M_PI * i / numSegs;
        double x1 = cos(theta1);
        double y1 = sin(theta1);
        double theta2 =  2.0 * M_PI * (i+1) / numSegs;
        double x2 = cos(theta2);
        double y2 = sin(theta2);

        vectNd_reset(&pos1);
        vectNd_reset(&pos2);
        vectNd_reset(&posm1);
        vectNd_set(&pos1, d1, x1*scale);
        vectNd_set(&pos1, d2, y1*scale);
        vectNd_set(&pos2, d1, x2*scale);
        vectNd_set(&pos2, d2, y2*scale);
        vectNd_set(&posm1, d1, xm1*scale);
        vectNd_set(&posm1, d2, ym1*scale);

        /* determine which sections should be drawn */
        int col1 = y1+c+0.5;
        int col2 = y2+c+0.5;
        int row1 = x1+r+0.5;
        int row2 = x2+r+0.5;
        int rowm1 = xm1+r+0.5;
        int colm1 = ym1+c+0.5;
        //printf("1) %g,%g -> %i,%i\n", x1, y1, row1, col1);
        //printf("2) %g,%g -> %i,%i\n", x2, y2, row2, col2);
        int cell1 = face_get_cell(&maze.faces[face], row1, col1);
        int cell2 = face_get_cell(&maze.faces[face], row2, col2);
        int cellm1 = face_get_cell(&maze.faces[face], rowm1, colm1);
        //printf("cell1=%i, cell2=%i, cellm1=%i\n", cell1, cell2, cellm1);


        /* add a joint between cylinders */
        if( cell2 != 0 ) {
            object *sph = object_alloc(maze.numDimensions, "sphere", "joint");
            snprintf(sph->name, sizeof(sph->name), "joint %i", i);
            object_add_obj(marker, sph);
            object_add_pos(sph, &pos2);
            object_add_size(sph, radius*scale);
            set_face_color(sph, face);
        } 
        if( cell1!=0 && cellm1 == 0 ) {
            object *sph = object_alloc(maze.numDimensions, "sphere", "joint");
            snprintf(sph->name, sizeof(sph->name), "joint %i-1", i);
            object_add_obj(marker, sph);
            object_add_pos(sph, &pos2);
            object_add_size(sph, radius*scale);
            set_face_color(sph, face);
        }

        if( cell1 == 0 || cell2 == 0 )
            continue;

        /* add an object from pos1 to pos2 */
        #if 0
        object *cyl = object_alloc(maze.numDimensions, "hcylinder", "marker piece");
        object_add_obj(marker, cyl);
        object_add_pos(cyl, &pos1);
        object_add_pos(cyl, &pos2);
        object_add_size(cyl, radius*scale);
        object_add_flag(cyl, 0);
        set_face_color(cyl, face);
        #else
        add_maze_marker_segment(marker, face, x1, y1, x2, y2, scale, radius);
        #endif
    }
    #else
    /* add central marker */
    vectNd centerPos;
    vectNd_calloc(&centerPos,maze.numDimensions);
    object *sph = object_alloc(maze.numDimensions, "sphere", "marker center");
    object_add_obj(marker, sph);
    //vectNd_set(&centerPos, d1, c*scale);
    //vectNd_set(&centerPos, d2, r*scale);
    object_add_pos(sph, &centerPos);
    object_add_size(sph, 3.0*radius*scale);
    set_face_color(sph, face);
    #endif

    return marker;
}

static void add_maze_marker_corner(object *marker, int face,
    double x1, double y1,
    double x2, double y2,
    double x3, double y3,
    double scale, double radius) {

    add_maze_marker_segment(marker, face, x1, y1, x2, y2, scale, radius);
    add_maze_marker_segment(marker, face, x2, y2, x3, y3, scale, radius);

    int d1 = maze.faces[face].d1;
    int d2 = maze.faces[face].d2;
    vectNd pos;
    vectNd_calloc(&pos, maze.numDimensions);
    vectNd_reset(&pos);
    vectNd_set(&pos, d1, x2*scale);
    vectNd_set(&pos, d2, y2*scale);

    object *sph = object_alloc(maze.numDimensions, "sphere", "joint");
    object_add_obj(marker, sph);
    object_add_pos(sph, &pos);
    object_add_size(sph, radius*scale);
    set_face_color(sph, face);
}

static object *make_maze_marker2(int face, double scale, int r, int c) {

    double radius = marker_radius;

    object *marker = object_alloc(maze.numDimensions, "cluster", "maze marker 2");
    object_add_flag(marker,4);

    #if 1
    double d=0.5;
    /* add corners */
    if( face_get_cell(&maze.faces[face], r-1, c-1) != 0 ) {
        add_maze_marker_corner(marker, face, -1, -d,
                                             -1, -1,
                                             -d, -1, scale, radius);
    }
    if( face_get_cell(&maze.faces[face], r-1, c+1) != 0 ) {
        add_maze_marker_corner(marker, face, -d, 1,
                                             -1, 1,
                                             -1, d, scale, radius);
    }
    if( face_get_cell(&maze.faces[face], r+1, c-1) != 0 ) {
        add_maze_marker_corner(marker, face, 1, -d,
                                             1, -1,
                                             d, -1, scale, radius);
    }
    if( face_get_cell(&maze.faces[face], r+1, c+1) != 0 ) {
        add_maze_marker_corner(marker, face, 1, d,
                                             1, 1,
                                             d, 1, scale, radius);
    }

    /* add sides of square */
    if( face_get_cell(&maze.faces[face], r, c-1) != 0 ) {
        add_maze_marker_segment(marker, face, -d,-1,
                                             d, -1, scale, radius);
    }
    if( face_get_cell(&maze.faces[face], r-1, c) != 0 ) {
        add_maze_marker_segment(marker, face, -1, -d,
                                             -1, d, scale, radius);
    }
    if( face_get_cell(&maze.faces[face], r, c+1) != 0 ) {
        add_maze_marker_segment(marker, face, -d, 1,
                                             d, 1, scale, radius);
    }
    if( face_get_cell(&maze.faces[face], r+1, c) != 0 ) {
        add_maze_marker_segment(marker, face, 1, -d,
                                             1, d, scale, radius);
    }
    #else
    /* add central marker */
    vectNd centerPos;
    vectNd_calloc(&centerPos,maze.numDimensions);
    object *sph = object_alloc(maze.numDimensions, "sphere", "marker center");
    object_add_obj(marker, sph);
    //vectNd_set(&centerPos, d1, c*scale);
    //vectNd_set(&centerPos, d2, r*scale);
    object_add_pos(sph, &centerPos);
    object_add_size(sph, 3.0*radius*scale);
    set_face_color(sph, face);
    #endif

    return marker;
}

static void add_maze_faces(object *puzzle, maze_t *maze, double edge_size) {

    printf("%s\n", __FUNCTION__);
    int dim = puzzle->dimensions;
    double scale = edge_size/maze->faces[0].rows;
    double epsilon = 0.01;

    /* for each face in maze */
    for(int face=0; face < maze->numFaces; ++face) {
        int d1 = maze->faces[face].d1;
        int d2 = maze->faces[face].d2;
        int rows = maze->faces[face].rows;
        int cols = maze->faces[face].cols;
        /* for high and low values in all dimensions not in face's plane */
        vectNd offset;
        vectNd_calloc(&offset,dim);
        vectNd counter;
        vectNd_calloc(&counter,dim-2);
        vectNd markerOffset;
        vectNd_calloc(&markerOffset, dim);
        int done = 0;
        //printf("face %i: d1=%i; d2=%i\n", face, d1, d2);
        for(int i=0; i<dim; ++i)
            vectNd_set(&offset,i,-1.0);
        vectNd_set(&offset,d1,0.0);
        vectNd_set(&offset,d2,0.0);
        while( !done ) {
            //vectNd_print(&counter,"\tcounter");

            /* create cluster */
            char faceName[64];
            snprintf(faceName,sizeof(faceName),"face %i,%i", d1, d2);
            object *faceCluster = object_alloc(dim, "cluster", faceName);
            object_add_flag(faceCluster,4);

            #if 1
            /* fill in face with hcube for each cell */
            for(int row=0; row<rows; ++row) {
                for(int col=0; col<cols; ++col) {
                    int cell = face_get_cell(&maze->faces[face], row, col);
                    if(cell==0)
                        continue;

                    /* create hcube for face cell */
                    char cellName[64];
                    snprintf(cellName,sizeof(cellName),"cell %i,%i for face %i", col, row, face);
                    object *box = object_alloc(dim, "hcube", cellName);

                    /* configure hcube */
                    vectNd cellPos;
                    vectNd_calloc(&cellPos,dim);
                    vectNd_reset(&cellPos);
                    vectNd_set(&cellPos, d1, row*scale);
                    vectNd_set(&cellPos, d2, col*scale);
                    object_add_pos(box, &cellPos);
                    //vectNd_print(&cellPos, "cellPos");

                    /* set up basis for hcube */
                    vectNd dir;
                    vectNd_calloc(&dir, dim);
                    for(int i=0; i<dim; ++i) {
                        vectNd_reset(&dir);
                        vectNd_set(&dir, i, 1);
                        object_add_dir(box, &dir);
                        object_add_size(box, scale);
                    }

                    set_face_color(box, face);

                    /* add cell to face */
                    object_add_obj(faceCluster, box);

                    vectNd_free(&dir);
                    vectNd_free(&cellPos);
                }
            }
            #endif // 1

            /* convert counter into offset vectors, skipping d1 and d2 */
            int k = 0, j = 0;
            vectNd_reset(&offset);
            vectNd_reset(&markerOffset);
            while(k < dim-2 && j < dim) {
                while( j == d1 || j == d2 ) {
                    ++j;
                    continue;
                }
                /* subtract 1 here to account for thickness of face */
                int dimK = maze->dimensions[k]-1;
                /* subtract 0.5 from counter to center at 0 */
                /* subtract 0.5 to shift both faces by half of their thickness */
                int countK = counter.v[k];
                double dist = dimK*countK+(countK?1:(-1))*epsilon*(maze->numFaces-face);
                double markerDist = (countK?1:-1)*(0.5+marker_radius+epsilon)*scale;
                vectNd_set(&offset, j, dist);
                vectNd_set(&markerOffset, j, markerDist);
                ++j;
                ++k;
            }

            #if 1
            /* add start/end markers */
            object* startMarker = make_maze_marker2(face, scale, maze->startPos[d1], maze->startPos[d2]);
            object_add_obj(faceCluster, startMarker);
            vectNd_set(&markerOffset, d1, maze->startPos[d1]*scale);
            vectNd_set(&markerOffset, d2, maze->startPos[d2]*scale);
            object_move(startMarker, &markerOffset);

            object* endMarker = make_maze_marker1(face, scale, maze->endPos[d1], maze->endPos[d2]);
            object_add_obj(faceCluster, endMarker);
            vectNd_set(&markerOffset, d1, maze->endPos[d1]*scale);
            vectNd_set(&markerOffset, d2, maze->endPos[d2]*scale);
            object_move(endMarker, &markerOffset);
            #endif // 1

            /* move face into position */
            vectNd scaledOffset;
            vectNd_calloc(&scaledOffset,dim);
            vectNd_scale(&offset, scale, &scaledOffset);
            object_move(faceCluster, &scaledOffset);
            //vectNd_print(&scaledOffset,"\tscaledOffset");

            /* update counter */
            j = 0;
            while(j < dim-2 && counter.v[j] == 1 ) {
                vectNd_set(&counter,j++,0);
            }
            if( j < dim-2 )
                vectNd_set(&counter,j,1);
            else
                done = 1;

            /* add face cluster to puzzle */
            object_add_obj(puzzle, faceCluster);
            faceCluster = NULL;
        }
        vectNd_free(&counter);
        vectNd_free(&offset);
    }
}

static void add_slider(object *puzzle, maze_t *maze, double edge_size, int frame, int frames) {
    
    double scale = edge_size / maze->faces[0].rows;

    int dimensions = maze->numDimensions;
    /* create cluster to add slider sub-pieces to */
    object *slider = object_alloc(dimensions, "cluster", "slider");
    object_add_flag(slider, 1);
    object_add_obj(puzzle, slider);

    /* for each dimension, add an hcube lengthened in that dimension */
    vectNd hcubeDir;
    vectNd_calloc(&hcubeDir, dimensions);
    vectNd offset;
    vectNd_calloc(&offset, dimensions);
    for(int d1=0; d1<dimensions; ++d1) {
        for(int d2=d1+1; d2<dimensions; ++d2) {
            object *obj = object_alloc(dimensions, "hcube", "movable slider part");
            snprintf(obj->name, sizeof(obj->name), "slider for %i,%i faces", d1, d2);
            object_add_obj(slider, obj);

            for(int i=0; i<dimensions; ++i) {
                if( i==d1 || i==d2 ) 
                    object_add_size(obj, scale);
                else
                    object_add_size(obj, 2.0*edge_size);
            
                vectNd_reset(&hcubeDir);
                vectNd_set(&hcubeDir, i, 1.0);
                object_add_dir(obj, &hcubeDir);
            }

            vectNd_reset(&offset);
    #if 0
            for(int i=0; i<dimensions; ++i)
                vectNd_set(&offset, i, -0.5*edge_size);
    #endif // 0
            object_add_pos(obj,&offset);
            obj->red = 0.8;
            obj->blue = 0.8;
            obj->green = 0.8;
        }
    }

    /* get location of center */
    vectNd pos;
    vectNd_calloc(&pos, dimensions);
    int pos1, pos2;
    double posW;
    /* perform moves */
    pos1 = frame / framesPerMove;
    pos2 = pos1+1;
    posW = (frame % framesPerMove) / (double)framesPerMove;
    if( pos1 < 0 ) {
        pos1 = 0;
        pos2 = 0;
        posW = 0.0;
    }
    if( pos2 >= maze->solution.posListNum ) {
        /* pause for final spin */
        pos1 = maze->solution.posListNum-1;
        pos2 = maze->solution.posListNum-1;
        posW = 0.0;
    }
    for(int i=0; i<dimensions; ++i) {
        vectNd_set(&pos, i,
            ((1.0-posW) * maze->solution.positions[pos1][i]
            + posW * maze->solution.positions[pos2][i]));
    }
    //vectNd_print(&pos, "slider pos is");
    printf("pos1: %i, pos2: %i, posW: %g\n", pos1, pos2, posW);
    vectNd_scale(&pos, scale, &pos);
    vectNd_print(&pos, "slider at");

    /* move cluster to correct location */
    object_move(slider,&pos);
}

int scene_setup(scene *scn, int dimensions, int frame, int frames, char *config)
{
    double t = frame/(double)frames;
    scene_init(scn, "maze-cube", dimensions);

    printf("Generating frame %i of %i scene '%s' (%.2f%% through animation).\n",
            frame, frames, scn->name, 100.0*t);

    /* basic setup */
    #if 1
    scn->bg_red = 0.09;
    scn->bg_green = 0.25;
    scn->bg_blue = 0.64;
    #endif /* 0 */
    scn->ambient.red = 0.0;
    scn->ambient.green = 0.0;
    scn->ambient.blue = 0.0;

    /* setup lighting */
    light *lgt=NULL;
    scene_alloc_light(scn,&lgt);
    lgt->type = LIGHT_AMBIENT;
    lgt->red = 0.3;
    lgt->green = 0.3;
    lgt->blue = 0.3;

    scene_alloc_light(scn,&lgt);
    #if 0
    lgt->type = LIGHT_POINT;
    vectNd_calloc(&lgt->pos,dimensions);
    vectNd_setStr(&lgt->pos,"0,40,0,-40");
    lgt->red = 300;
    lgt->green = 300;
    lgt->blue = 300;
    #else
    lgt->type = LIGHT_DIRECTIONAL;
    vectNd_calloc(&lgt->dir,dimensions);
    vectNd_setStr(&lgt->dir,"0,-1,0,-1");
    lgt->red = 0.4;
    lgt->green = 0.4;
    lgt->blue = 0.4;
    #endif

    #if 1
    scene_alloc_light(scn,&lgt);
    lgt->type = LIGHT_DIRECTIONAL;
    vectNd_calloc(&lgt->dir,dimensions);
    vectNd_setStr(&lgt->dir,"-180,-40,0,0");
    lgt->red = 0.15;
    lgt->green = 0.15;
    lgt->blue = 0.15;

    scene_alloc_light(scn,&lgt);
    lgt->type = LIGHT_DIRECTIONAL;
    vectNd_calloc(&lgt->dir,dimensions);
    vectNd_setStr(&lgt->dir,"0,-40,140,0");
    lgt->red = 0.15;
    lgt->green = 0.15;
    lgt->blue = 0.15;
    #endif // 1

    vectNd temp;
    vectNd_calloc(&temp,dimensions);

    #if 1
    /* add background */
    object *ground=NULL;
    scene_alloc_object(scn,dimensions,&ground,"hplane");
    vectNd_reset(&temp);
    vectNd_set(&temp,1,-45);
    object_add_pos(ground,&temp);  /* position */
    vectNd_reset(&temp);
    vectNd_set(&temp,1,1);
    object_add_dir(ground,&temp);  /* normal */
    double grassScaling = 0.75;
    ground->red = 0.0225 * grassScaling;
    ground->green = 1.0 * grassScaling;
    ground->blue = 0.04 * grassScaling;
    #endif /* 0 */

    #if 1
    /* add mirrors */
    double mirror_dist = 66;
    /* positive z */
    add_mirror(scn, dimensions, 2, mirror_dist);
    /* negative x */
    add_mirror(scn, dimensions, 0, -mirror_dist);
    #if 0
    if( dimensions > 3 ) {
        /* positive w */
        add_mirror(scn, dimensions, 3, mirror_dist);
    }
    #endif /* 1? */
    #endif /* 0 */

    /* cluster will contain puzzle and be used to rotate it */
    object *clstr = NULL;
    scene_alloc_object(scn,dimensions,&clstr,"cluster");
    object_add_flag(clstr, 8);  /* set number of clusters */

    /* create objects */
    #if 0
    /* add placeholder cube */
    /* use object_alloc instead of scene_alloc_object when the object will be
     * added to a cluster */
    vectNd hcubeDir;
    vectNd_calloc(&hcubeDir, dimensions);
    object *obj = object_alloc(dimensions, "hcube", "placeholder cube");
    for(int i=0; i<dimensions; ++i)
        object_add_size(obj, edge_size);
    for(int i=0; i<dimensions; ++i) {
        vectNd_reset(&hcubeDir);
        vectNd_set(&hcubeDir, i, 1.0);
        object_add_dir(obj, &hcubeDir);
    }
    vectNd offset;
    vectNd_calloc(&offset, dimensions);
    #if 0
    for(int i=0; i<dimensions; ++i)
        vectNd_set(&offset, i, -edge_size);
    #endif // 0
    object_add_pos(obj,&offset);
    obj->red = 0.1;
    obj->blue = 0.8;
    obj->green = 0.1;
    object_add_obj(clstr, obj);
    #else
    add_maze_faces(clstr, &maze, edge_size);
    #endif
    #if 1
    add_slider(clstr, &maze, edge_size, frame-framesPerSpin, frames);
    #endif // 1

    #if 1
    /* move puzzle to be centered at origin */
    double scale = edge_size/maze.faces[0].rows;
    vectNd centeringOffset;
    vectNd_calloc(&centeringOffset, dimensions);
    for(int i=0; i<dimensions; ++i) {
        vectNd_set(&centeringOffset, i, -0.5*scale*maze.dimensions[i]);
    }
    object_move(clstr, &centeringOffset);
    vectNd_print(&centeringOffset, "centeringOffset");
    #endif // 1

    #if 1
    /* rotate puzzle into view orientation */
    vectNd rotateCenter, rotV1, rotV2;
    vectNd_calloc(&rotateCenter, dimensions);
    vectNd_calloc(&rotV1, dimensions);
    vectNd_calloc(&rotV2, dimensions);
    for(int i=0; i<3; ++i) {
        vectNd_set(&rotV1,i,1.0);
        vectNd_set(&rotV2,i,1.0);
    }
    vectNd_set(&rotV1,1,0.0);
    double angle = (M_PI/2.0) - atan(1.0/sqrt(dimensions-1));
    object_rotate2(clstr, &rotateCenter, &rotV1, &rotV2, angle);
    angle = 35*M_PI/180.0;
    //angle = 90*M_PI/180.0;
    int endFrame = frame - framesPerSpin - (maze.solution.posListNum-1)*framesPerMove;
    if( frame < framesPerSpin )
        angle += frame*2.0*M_PI/framesPerSpin;
    else if( endFrame >= 0 )
        angle += endFrame*2.0*M_PI/framesPerSpin;
    printf("rotating %g degrees in x,z plane (frame=%i).\n", angle*180/M_PI, frame);
    object_rotate(clstr, &rotateCenter, 0, 2, angle);
    #endif // 1
    
    /* zero out camera */
    camera_reset(&scn->cam);

    /* move camera into position */
    vectNd viewPoint;
    vectNd viewTarget;
    vectNd up_vect;
    vectNd lookVec;
    vectNd camCentering;
    vectNd_calloc(&viewPoint,dimensions);
    vectNd_calloc(&viewTarget,dimensions);
    vectNd_calloc(&up_vect,dimensions);
    vectNd_calloc(&lookVec, dimensions);
    vectNd_calloc(&camCentering, dimensions);

    vectNd_setStr(&viewTarget,"-5,-5,20,0");
    vectNd_setStr(&viewPoint,"160,45,-120,0");
    //vectNd_setStr(&viewTarget,"-5,-5,20,0.0");
    //vectNd_setStr(&viewPoint,"160,45,-120,0.0");
    vectNd_copy(&camCentering, &centeringOffset);
    for(int i=0; i<3; ++i)
        vectNd_set(&camCentering, i, 0.0);
    vectNd_add(&viewTarget, &camCentering, &viewTarget);
    vectNd_add(&viewPoint, &camCentering, &viewPoint);
    vectNd_set(&up_vect,1,1);  /* 0,1,0,0... */
    vectNd_sub(&viewPoint, &viewTarget, &lookVec);
    //vectNd_rotate2(&viewPoint, &viewTarget, &lookVec, &up_vect, 10.0*M_PI/180.0, &viewPoint);
    camera_set_aim(&scn->cam, &viewPoint, &viewTarget, &up_vect, 0.0);
    camera_set_flip(&scn->cam, 1, 0);
    vectNd_free(&up_vect);
    vectNd_free(&viewPoint);
    vectNd_free(&viewTarget);

    return 1;
}

int scene_cleanup() {
    /* If any persistent resources were allocated,
     * they should be freed here. */
    maze_free(&maze);
    return 0;
}
