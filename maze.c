/*
 * maze.c
 * PuzzleMaze: A 3D puzzle model generator
 *
 * Copyright (c) 2020 Bryan Franklin. All rights reserved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "maze.h"

int face_init(face_t *face, int *sizes, int d1, int d2) {
    face->d1 = d1;
    face->d2 = d2;
    face->rows = sizes[d1];
    face->cols = sizes[d2];
    face->cells = calloc(face->rows*face->cols,sizeof(char));
    memset(face->cells, 1, face->rows*face->cols);

    return 0;
}


int face_free(face_t *face) {
    free(face->cells); face->cells = NULL;
    memset(face, '\0', sizeof(*face));

    return 0;
}

static int face_set_cell(face_t *face, int row, int col, int value) {
    int pos = row*face->cols + col;
    int old = face->cells[pos];
    face->cells[pos] = value;
    return old;
}


int face_get_cell(face_t *face, int row, int col) {
    if( row < 0 || row>=face->rows
        || col < 0 || col >=face->cols ) {
        fprintf(stderr, "Request for out of bounds face cell %i,%i (size: %ix%i)\n",
            row, col, face->rows, face->cols);
        return 0;
    }
    int pos = row*face->cols + col;
    return face->cells[pos];
}


static int pos_list_init(position_list_t *list, int numDimensions) {
    list->numDimensions = numDimensions;
    list->posListCap = 10;
    list->posListNum = 0;
    list->positions = calloc(list->posListCap,sizeof(position_t));
    return (list->positions!=NULL);
}


static int pos_list_free(position_list_t *list) {
    free(list->positions); list->positions=NULL;
    memset(list,'\0',sizeof(*list));
    return 0;
}


static int pos_list_push(position_list_t *list, position_t pos) {
    /* add to list of restart locations */
    if( list->posListNum == list->posListCap ) {
        int newCap = list->posListCap*2+1;
        void *tmp = realloc(list->positions,newCap*sizeof(int*));
        if( tmp==NULL ) {
            perror("realloc");
            return 0;
        }
        list->positions = tmp;
        list->posListCap = newCap;
    }
    list->positions[list->posListNum] = calloc(list->numDimensions,sizeof(int));
    memcpy(list->positions[list->posListNum], pos, list->numDimensions*sizeof(int));
    ++list->posListNum;
    return 1;
}


static int pos_list_pop(position_list_t *list, position_t pos) {
    if( list->posListNum < 0 )
        return 0;
    memcpy(pos,list->positions[list->posListNum-1],list->numDimensions*sizeof(int));
    free(list->positions[list->posListNum-1]); list->positions[list->posListNum-1]=NULL;
    --list->posListNum;
    return 1;
}


static int pos_list_random(position_list_t *list, position_t pos) {
    if( list->posListNum <= 0 )
        return 0;
    int which = rand()%list->posListNum;
    memcpy(pos, list->positions[which], list->numDimensions*sizeof(int));
    return 1;
}


static int pos_list_rfind(position_list_t *list, position_t pos) {
    for(int i=list->posListNum-1; i>=0; --i) {
        if( memcmp(list->positions[i], pos, list->numDimensions*sizeof(int)) == 0 )
        return i;
    }
    return -1;
}


int position_increment(maze_t *maze, position_t pos) {
    int done = 0;
    int j=0;
    while(j<maze->numDimensions && pos[j]==maze->dimensions[j]-1) {
        pos[j++] = 1;
    }
    if( j < maze->numDimensions )
        ++pos[j];
    else
        done = 1;
    return done;
}


int maze_init(maze_t *maze, int numDimensions, int *sizes) {

    /* initialize maze structure */
    printf("Initializing %iD maze.\n", numDimensions);
    maze->numDimensions = numDimensions;
    maze->numFaces = (numDimensions*(numDimensions-1))/2;
    maze->dimensions = calloc(numDimensions,sizeof(int));
    for(int i=0; i<numDimensions; ++i) {
        if( (sizes[i]%2) == 0 ) {
            fprintf(stderr, "Warning: dimension sizes must be odd, adjusting dimension %i from %i to %i.\n", i, sizes[i], sizes[i]+1);
            sizes[i] += 1;
        }
        maze->dimensions[i] = sizes[i];
    }

    /* allocate start and end positions */
    maze->startPos = calloc(numDimensions,sizeof(int));
    maze->endPos = calloc(numDimensions,sizeof(int));

    /* allocate and initialize faces */
    maze->faces = calloc(maze->numFaces,sizeof(face_t));
    int face=0;
    for(int d1 = 0; d1 < numDimensions; ++d1) {
        for(int d2 = d1+1; d2 < numDimensions; ++d2) {
            printf("  face %i -> dimensions %i,%i\n", face, d1, d2);
            face_init(&maze->faces[face++], sizes, d1, d2);
        }
    }

    return 0;
}


int maze_free(maze_t *maze) {
    for(int face=0; face<maze->numFaces; ++face) {
        face_free(&maze->faces[face]);
    }
    free(maze->startPos); maze->startPos=NULL;
    free(maze->endPos); maze->endPos=NULL;
    free(maze->faces); maze->faces=NULL;
    memset(maze,'\0', sizeof(*maze));

    return 0;
}


static int maze_clear_cell(maze_t *maze, int *pos) {
    int ret = 0;

    for(int face = 0; face < maze->numFaces; ++face) {
        int row = pos[maze->faces[face].d1];
        int col = pos[maze->faces[face].d2];
        ret |= face_set_cell(&maze->faces[face], row, col, 0);
    }

    return ret;
}


static int maze_allow_clear(maze_t *maze, int *pos, int move) {
    /* check for border violations */
    if( move < 0 ) {
        if( pos[-move-1] <= 2 )
            return 0;
    } else {
        if( pos[move-1] >= maze->dimensions[move-1]-2)
            return 0;
    }

    /* determine where move would end up */
    int *midPos = calloc(maze->numDimensions, sizeof(int));
    int *nextPos = calloc(maze->numDimensions, sizeof(int));
    memcpy(midPos, pos, maze->numDimensions*sizeof(int));
    memcpy(nextPos, pos, maze->numDimensions*sizeof(int));
    if( move >= 0 )
        midPos[move-1] += 1;
    else
        midPos[-move-1] -= 1;
    if( move >= 0 )
        nextPos[move-1] += 2;
    else
        nextPos[-move-1] -= 2;

    /* check to see if that position can be safely cleared */
    int allowed = 1;
    int breaksWall = 0;
    for(int face = 0; allowed && face < maze->numFaces; ++face) {
        if( abs(move)-1 != maze->faces[face].d1
            && abs(move)-1 != maze->faces[face].d2 )
            continue;
        int row = nextPos[maze->faces[face].d1];
        int col = nextPos[maze->faces[face].d2];
        if( row >= maze->faces[face].rows-1 
            || col >= maze->faces[face].cols-1 )
            continue;
        int filledCells = 0;
        for(int i=-1; i<=1; ++i) {
            for(int j=-1; j<=1; ++j) {
                if( face_get_cell(&maze->faces[face], row+i, col+j)!=0 )
                    ++filledCells;
            }
        }
        int destCell = face_get_cell(&maze->faces[face], row, col);
        int midRow = midPos[maze->faces[face].d1];
        int midCol = midPos[maze->faces[face].d2];
        int midCell = face_get_cell(&maze->faces[face], midRow, midCol);
        if( !(filledCells == 9 || (midCell==0 && destCell==0)) )
            allowed = 0;
        if( filledCells == 9 )
            breaksWall = 1;
    }

    return (allowed&&breaksWall);
}


static int maze_gen_step(maze_t *maze, int *pos) {

    /* clear cell at position pos */
    maze_clear_cell(maze,pos);

    /* enumerate neighbors that can be moved to */
    int *validMoves = calloc(maze->numDimensions*2, sizeof(int));
    int numMoves = 0;
    for(int i=0; i<maze->numDimensions; ++i) {
        if( maze_allow_clear(maze,pos,i+1) ) {
            validMoves[numMoves++] = i+1;
        }
        if( maze_allow_clear(maze,pos,-(i+1)) ) {
            validMoves[numMoves++] = -(i+1);
        }
    }

    /* if valid move to neighbor found */
    if( numMoves > 0) {
        /* pick move randomly */
        int move = validMoves[rand()%numMoves];

        /* recurse to new position */
        int *nextPos = calloc(maze->numDimensions,sizeof(int));
        memcpy(nextPos,pos,maze->numDimensions*sizeof(int));
        if( move > 0 )
            nextPos[move-1] += 1;
        else
            nextPos[-move-1] -= 1;
        maze_clear_cell(maze,nextPos);
        if( move > 0 )
            nextPos[move-1] += 1;
        else
            nextPos[-move-1] -= 1;
        maze_gen_step(maze,nextPos);
        free(nextPos); nextPos=NULL;

        /* recurse to current position, to check for remaining valid moves */
        maze_gen_step(maze,pos);
    }

    return 0;
}


int maze_unfinished(maze_t *maze) {

    /* check all faces for 2x2 regions of all uncleared cells */
    for(int face = 0; face < maze->numFaces; ++face) {
        int rows = maze->faces[face].rows;
        int cols = maze->faces[face].cols;
        for(int row = 1; row < rows-1; ++row) {
            for(int col = 1; col < cols-1; ++col) {
                if( face_get_cell(&maze->faces[face],row,col)
                    && face_get_cell(&maze->faces[face],row+1,col)
                    && face_get_cell(&maze->faces[face],row,col+1)
                    && face_get_cell(&maze->faces[face],row+1,col+1) ) {
                    /* unfinished section found */
                    //printf("position %i,%i of face %i is unfinished.\n", row, col, face);
                    return 1;
                }
            }
        }
    }

    /* no unfinished sections found */
    return 0;
}


int maze_position_clear(maze_t *maze, int *pos) {
    int isClear = 1;
    for(int face = 0; isClear && face < maze->numFaces; ++face) {
        int row = pos[maze->faces[face].d1];
        int col = pos[maze->faces[face].d2];
        if( face_get_cell(&maze->faces[face], row, col)!=0 )
            isClear = 0;
    }

    return isClear;
}


int maze_get_restart_location(maze_t *maze, int *pos) {

    /* start position counter at all 1s */
    for(int i = 0; i<maze->numDimensions; ++i) {
        pos[i] = 1;
    }

    /* make list of valid moves */
    int *moves = calloc(2*maze->numDimensions, sizeof(int));
    for(int i=0; i<2*maze->numDimensions; i+=2) {
        int dir = i/2;
        moves[i] = dir+1;
        moves[i+1] = -(dir+1);
        //printf("moves[%i & %i] = %i & %i\n", i, i+1, moves[i], moves[i+1]);
    }

    /* create list of possible locations */
    position_list_t posList;
    pos_list_init(&posList, maze->numDimensions);

    /* for every cleared cell */
    int done = 0;
    while( !done ) {

        if( maze_position_clear(maze, pos) ) {

            #if 0
            printf("clear position at: ");
            for(int i=0; i<maze->numDimensions; ++i)
                printf("%i ", pos[i]);
            printf("\n");
            #endif // 0

            /* check all moves from current pos */
            for(int m = 0; m<2*maze->numDimensions; ++m) {

                //printf("\ttrying move %i\n", moves[m]);
                /* if a valid move from it exists */
                if( maze_allow_clear(maze, pos, moves[m]) ) {

                    /* add to list of restart locations */
                    pos_list_push(&posList, pos);
                }
            }
        }

        /* update pos */
        done = position_increment(maze, pos);
    }

    if( !pos_list_random(&posList, pos)) {
        fprintf(stderr, "No restart locations found!\n");
    }

    /* free list of possible locations */
    free(moves); moves = NULL;
    pos_list_free(&posList);

    return 0;
}


int maze_pick_goals(maze_t *maze) {

    #if 1
    /* start position counter at all 1s */
    int *pos = calloc(maze->numDimensions,sizeof(int));
    for(int i = 0; i<maze->numDimensions; ++i) {
        pos[i] = 1;
    }

    /* create list of possible locations */
    position_list_t posList;
    pos_list_init(&posList, maze->numDimensions);

    /* for every cleared cell */
    int done = 0;
    while( !done ) {

        if( maze_position_clear(maze, pos) ) {

            printf("Adding: ");
            for(int i=0; i<maze->numDimensions; ++i) {
                printf("%i ", pos[i]);
            }
            printf("\n");

            pos_list_push(&posList, pos);
        }

        /* update pos */
        done = position_increment(maze, pos);
    }

    #if 1
    if( posList.posListNum > 2 ) {
        printf("posListNum=%i\n", posList.posListNum);
        /* pick random start location */
        pos_list_random(&posList, maze->startPos);

        /* pick random end location */
        pos_list_random(&posList, maze->endPos);
    } else {
        printf("Insufficient free spaces to pick start/end locations.\n");
    }
    #endif // 1
    pos_list_free(&posList);
    free(pos); pos = NULL;
    #endif // 1

    return 0;
}


int maze_solve(maze_t *maze) {
    return 0;
}


int maze_generate(maze_t *maze) {
    printf("Generating %iD maze.\n", maze->numDimensions);

    /* set starting point */
    int *start = calloc(maze->numDimensions,sizeof(int));
    for(int i=0; i<maze->numDimensions; ++i) {
        start[i] = 1;
    }

    /* recursively clear cells */
    int ret = maze_gen_step(maze,start);
    free(start); start=NULL;

    int done = 0;
    int retries = 10;
    while(!done && --retries) {
        /* pick start and end locations */
        maze_pick_goals(maze);

        printf("Start: ");
        for(int i=0; i<maze->numDimensions; ++i) {
            printf("%i ", maze->startPos[i]);
        }
        printf("\nEnd:   ");
        for(int i=0; i<maze->numDimensions; ++i) {
            printf("%i ", maze->endPos[i]);
        }
        printf("\n");

        /* find and record solution */
        maze_solve(maze);

        done = 1;
    }

    #if 1
    /* while not completely full (i.e., any 2x2 region is full) */
    int *restartPos = calloc(maze->numDimensions,sizeof(int));
    retries = 10;
    while(maze_unfinished(maze) && --retries) {
        maze_get_restart_location(maze, restartPos);

        /* try using as an unreachable starting point */
        ret = maze_gen_step(maze, restartPos);
    }
    #endif 
    free(start); start=NULL;

    return ret;
}


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

static void maze_export_stl_cube(FILE *fp, int x, int y, int z, char face_mask, double scale) {
    double d = scale/2.0;

    if( (face_mask & (1<<0)) == 0) {
        /* left (-x) */
        maze_export_stl_triangle(fp, x-d, y-d, z-d,
                x-d, y-d, z+d,
                x-d, y+d, z+d,
                -1, 0, 0);
        maze_export_stl_triangle(fp, x-d, y-d, z-d,
                x-d, y+d, z+d,
                x-d, y+d, z-d,
                -1, 0, 0);
    }
    if( (face_mask & (1<<1)) == 0) {
        /* right (+x) */
        maze_export_stl_triangle(fp, x+d, y-d, z-d,
                x+d, y+d, z+d,
                x+d, y-d, z+d,
                1, 0, 0);
        maze_export_stl_triangle(fp, x+d, y-d, z-d,
                x+d, y+d, z-d,
                x+d, y+d, z+d,
                1, 0, 0);
    }

    if( (face_mask & (1<<2)) == 0) {
        /* front (-y) */
        maze_export_stl_triangle(fp, x-d, y-d, z-d,
                x+d, y-d, z+d,
                x-d, y-d, z+d,
                0, -1, 0);
        maze_export_stl_triangle(fp, x-d, y-d, z-d,
                x+d, y-d, z-d,
                x+d, y-d, z+d,
                0, -1, 0);
    }
    if( (face_mask & (1<<3)) == 0) {
        /* back (+y) */
        maze_export_stl_triangle(fp, x-d, y+d, z-d,
                x-d, y+d, z+d,
                x+d, y+d, z+d,
                0, 1, 0);
        maze_export_stl_triangle(fp, x-d, y+d, z-d,
                x+d, y+d, z+d,
                x+d, y+d, z-d,
                0, 1, 0);
    }

    if( (face_mask & (1<<4)) == 0) {
        /* top (+z face) */
        maze_export_stl_triangle(fp, x-d, y-d, z+d,
                x+d, y-d, z+d,
                x+d, y+d, z+d,
                0, 0, 1);
        maze_export_stl_triangle(fp, x-d, y+d, z+d,
                x-d, y-d, z+d,
                x+d, y+d, z+d,
                0, 0, 1);
    }
    if( (face_mask & (1<<5)) == 0) {
        /* bottom (-z) */
        maze_export_stl_triangle(fp, x-d, y-d, z-d,
                x+d, y+d, z-d,
                x+d, y-d, z-d,
                0, 0, -1);
        maze_export_stl_triangle(fp, x-d, y+d, z-d,
                x+d, y+d, z-d,
                x-d, y-d, z-d,
                0, 0, -1);
    }
}


int maze_export_stl(maze_t *maze, char *filename) {
    if( maze->numDimensions != 3 ) {
        fprintf(stderr,"%s: STL export is only supported for 3D mazes.\n", __FUNCTION__);
        return -1;
    }
    printf("Exporting %iD maze to STL files `%s`.\n", maze->numDimensions, filename);

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
                    maze_export_stl_cube(fp, x1, y1, z1, mask1, scale);
                    maze_export_stl_cube(fp, x2, y2, z2, mask2, scale);
                }
            }
        }
    }

    /* close solid */
    fprintf(fp,"endsolid puzzle\n");

    /* close file */
    fclose(fp); fp=NULL;

    return 0;
}
