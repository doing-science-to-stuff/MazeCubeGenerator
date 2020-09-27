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
    memcpy(maze->dimensions,sizes,numDimensions*sizeof(int));
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


static int maze_clear_cell(maze_t *maze, int *pos) {
    int ret = 0;
    printf("  clearing position ");
    for(int i=0; i<maze->numDimensions; ++i)
        printf("%i ", pos[i]);
    printf("\n");
    for(int face = 0; face < maze->numFaces; ++face) {
        int row = pos[maze->faces[face].d1];
        int col = pos[maze->faces[face].d2];
        #if 0
        if( face_get_cell(&maze->faces[face], row, col)!=0 )
            printf("    clearing position %i,%i in face %i\n", row, col, face);
        #endif /* 0 */
        ret |= face_set_cell(&maze->faces[face], row, col, 0);
    }

    return ret;
}


static int maze_allow_clear(maze_t *maze, int *pos, int move) {
    /* check for border violations */
    if( move < 0 ) {
        if( pos[-(move-1)] <= 1 )
            return 0;
    } else {
        if( pos[move-1] >= maze->dimensions[move-1])
            return 0;
    }

    /* determine where move would end up */
    int *nextPos = calloc(maze->numDimensions, sizeof(int));
    memcpy(nextPos, pos, maze->numDimensions*sizeof(int));
    if( move >= 0 )
        nextPos[move-1] += 1;
    else
        nextPos[-(move-1)] -= 1;

    /* check to see if that position can be savely cleared */
    int allowed = 1;
    for(int face = 0; allowed && face < maze->numFaces; ++face) {
        int row = nextPos[maze->faces[face].d1];
        int col = nextPos[maze->faces[face].d2];
        int filledCells = 0;
        for(int i=-1; i<=1; ++i) {
            for(int j=-1; j<=1; ++j) {
                if( face_get_cell(&maze->faces[face], row+i, col+j)!=0 )
                    ++filledCells;
            }
        }
        if( filledCells < 7 )
            allowed = 0;
    }

    return allowed;
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
        if( move >= 0 )
            nextPos[move-1] += 1;
        else
            nextPos[-(move-1)] -= 1;
        maze_gen_step(maze,nextPos);

        /* recurse to current position, to check for remaining valid moves */
        maze_gen_step(maze,pos);
    }

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
    return maze_gen_step(maze,start);
}


int maze_export_stl(maze_t *maze, char *filename) {
    if( maze->numDimensions != 3 ) {
        fprintf(stderr,"%s: STL export is only supported for 3D mazes.\n", __FUNCTION__);
        return -1;
    }
    printf("Exporting %iD maze to STL files `%s`.\n", maze->numDimensions, filename);

    fprintf(stderr,"%s: STL export not implemented yet!\n", __FUNCTION__);

    return 0;
}
