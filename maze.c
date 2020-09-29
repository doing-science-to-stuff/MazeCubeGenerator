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
    int pos = row*face->cols + col;
    return face->cells[pos];
}


int maze_init(maze_t *maze, int numDimensions, int *sizes) {

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

    /* check to see if that position can be savely cleared */
    int allowed = 1;
    int breaksWall = 0;
    for(int face = 0; allowed && face < maze->numFaces; ++face) {
        if( abs(move)-1 != maze->faces[face].d1
            && abs(move)-1 != maze->faces[face].d2 )
            continue;
        int row = nextPos[maze->faces[face].d1];
        int col = nextPos[maze->faces[face].d2];
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

    /* enumerate neighbor that can be moved to */
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
        for(int row = 0; row < rows-1; ++row) {
            for(int col = 0; col < cols-1; ++col) {
                if( face_get_cell(&maze->faces[face],row,col)
                    && face_get_cell(&maze->faces[face],row+1,col)
                    && face_get_cell(&maze->faces[face],row,col+1)
                    && face_get_cell(&maze->faces[face],row+1,col+1) ) {
                    /* unfinished section found */
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
        moves[i] = i+1;
        moves[i] = -(i+1);
    }

    /* create list of possible locations */
    int posListCap = 10;
    int posListNum = 0;
    int **posList = calloc(posListCap,sizeof(int*));

    /* for every cleared cell */
    int done = 0;
    while( !done ) {

        if( maze_position_clear(maze, pos) ) {

            /* check all moves from current pos */
            for(int m = 0; m<2*maze->numDimensions; ++m) {

                /* if a valid move from it exists */
                if( maze_allow_clear(maze, pos, moves[m]) ) {

                    /* add to list of restart locations */
                    if( posListNum == posListCap ) {
                        int newCap = posListCap*2+1;
                        void *tmp = realloc(posList,newCap*sizeof(int*));
                        if( tmp==NULL ) {
                            perror("calloc");
                            exit(1);
                        }
                        posList = tmp;
                        posListCap = newCap;
                    }
                    posList[posListNum] = calloc(maze->numDimensions,sizeof(int));
                    memcpy(posList[posListNum], pos, maze->numDimensions*sizeof(int));
                    ++posListNum;
                }
            }
        }

        /* update pos */
        int j=0;
        while(pos[j]==maze->dimensions[j]-1) {
            pos[j++] = 1;
        }
        if( j < maze->numDimensions )
            ++pos[j];
        else
            done = 1;
    }

    /* pick random restart location */
    int which = rand()%posListNum;
    memcpy(pos, posList[which], maze->numDimensions*sizeof(int));

    /* free list of possible locations */
    free(moves); moves = NULL;
    for(int i=0; i<posListNum; ++i) {
        free(posList[i]); posList[i] = NULL;
    }
    free(posList); posList = NULL;

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

    /* pick start and end locations */

    /* find and record solution */

    /* while not complely full (i.e., any 2x2 region is full) */
    int *restartPos = calloc(maze->numDimensions,sizeof(int));
    while(maze_unfinished(maze)) {
        maze_get_restart_location(maze, restartPos);

        /* try using as an unreachable starting point */
        ret = maze_gen_step(maze, restartPos);
    }
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

static void maze_export_stl_cube(FILE *fp, int x, int y, int z, double scale) {
    double d = scale/2.0;
#if 1
    /* top (+z face) */
    maze_export_stl_triangle(fp, x-d, y-d, z+d,
                                 x+d, y-d, z+d,
                                 x+d, y+d, z+d,
                                 0, 0, 1);
    maze_export_stl_triangle(fp, x-d, y+d, z+d,
                                 x-d, y-d, z+d,
                                 x+d, y+d, z+d,
                                 0, 0, 1);
    /* bottom (-z) */
    maze_export_stl_triangle(fp, x-d, y-d, z-d,
                                 x+d, y+d, z-d,
                                 x+d, y-d, z-d,
                                 0, 0, -1);
    maze_export_stl_triangle(fp, x-d, y+d, z-d,
                                 x+d, y+d, z-d,
                                 x-d, y-d, z-d,
                                 0, 0, -1);
#endif /* 1 */

#if 1
    /* front (-y) */
    maze_export_stl_triangle(fp, x-d, y-d, z-d,
                                 x+d, y-d, z+d,
                                 x-d, y-d, z+d,
                                 0, -1, 0);
    maze_export_stl_triangle(fp, x-d, y-d, z-d,
                                 x+d, y-d, z-d,
                                 x+d, y-d, z+d,
                                 0, -1, 0);
    /* back (+y) */
    maze_export_stl_triangle(fp, x-d, y+d, z-d,
                                 x-d, y+d, z+d,
                                 x+d, y+d, z+d,
                                 0, 1, 0);
    maze_export_stl_triangle(fp, x-d, y+d, z-d,
                                 x+d, y+d, z+d,
                                 x+d, y+d, z-d,
                                 0, 1, 0);
#endif /* 1 */

#if 1
    /* left (-x) */
    maze_export_stl_triangle(fp, x-d, y-d, z-d,
                                 x-d, y-d, z+d,
                                 x-d, y+d, z+d,
                                 -1, 0, 0);
    maze_export_stl_triangle(fp, x-d, y-d, z-d,
                                 x-d, y+d, z+d,
                                 x-d, y+d, z-d,
                                 -1, 0, 0);
    /* right (+x) */
    maze_export_stl_triangle(fp, x+d, y-d, z-d,
                                 x+d, y+d, z+d,
                                 x+d, y-d, z+d,
                                 1, 0, 0);
    maze_export_stl_triangle(fp, x+d, y-d, z-d,
                                 x+d, y+d, z-d,
                                 x+d, y+d, z+d,
                                 1, 0, 0);
#endif /* 1 */
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
                    maze_export_stl_cube(fp, x1, y1, z1, scale);
                    maze_export_stl_cube(fp, x2, y2, z2, scale);
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
