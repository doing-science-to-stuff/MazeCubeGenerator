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
    if( list->posListNum <= 0 )
        return 0;
    if( pos!=NULL )
        memcpy(pos,list->positions[list->posListNum-1],list->numDimensions*sizeof(int));
    if( list->posListNum > 0 ) {
        free(list->positions[list->posListNum-1]); list->positions[list->posListNum-1]=NULL;
        --list->posListNum;
    }
    return 1;
}


static int pos_list_clear(position_list_t *list) {
    while(pos_list_pop(list,NULL));
    return 0;
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
    for(int i=0; i<numDimensions; ++i) {
        maze->startPos[i] = -1;
        maze->endPos[i] = -1;
    }

    /* initialize solution */
    pos_list_init(&maze->reachable, numDimensions);
    pos_list_init(&maze->solution, numDimensions);

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
    pos_list_push(&maze->reachable, pos);

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
    }

    /* create list of possible locations */
    position_list_t posList;
    pos_list_init(&posList, maze->numDimensions);

    /* for every cleared cell */
    int done = 0;
    while( !done ) {

        /* check for (and reject) positions with any even coordinates */
        int allOdd = 1;
        for(int i=0; allOdd && i<maze->numDimensions; ++i) {
            if( (pos[i]%2) == 0)
                allOdd = 0;
        }

        if( allOdd && maze_position_clear(maze, pos) ) {
            /* check all moves from current pos */
            for(int m = 0; m<2*maze->numDimensions; ++m) {

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

    if( maze->reachable.posListNum < 2 )
        return 0;

    pos_list_random(&maze->reachable, maze->startPos);
    pos_list_random(&maze->reachable, maze->endPos);

    return 1;
}


static int* validMoves = NULL;
static int maze_solve_dfs(maze_t *maze, position_t pos) {

    pos_list_push(&maze->solution, pos);

    /* check for endPos */
    if( memcmp(maze->endPos, pos, maze->numDimensions*sizeof(int))==0 ) {
        return 1;
    }

    /* copy position */
    position_t newPos;
    newPos = calloc(maze->numDimensions,sizeof(int));

    /* enumerate moves */
    for(int i=0; i<maze->numDimensions*2; ++i ) {
        /* update position */
        memcpy(newPos, pos, maze->numDimensions*sizeof(int));
        int move = validMoves[i];
        if( move > 0 )
            ++newPos[move-1];
        else
            --newPos[(-move)-1];

        /* if move valid */
        if( !maze_position_clear(maze,newPos) ) {
            continue;
        }

        /* check solution for next position */
        if( pos_list_rfind(&maze->solution, newPos) >= 0 ) {
            continue;
        }

        /* recurse */
        if( maze_solve_dfs(maze, newPos) )
            return 1;
    }

    pos_list_pop(&maze->solution, NULL);
    return 0;
}


int maze_solve(maze_t *maze) {

    /* initialize validMoves */
    validMoves = calloc(2*maze->numDimensions, sizeof(int));
    for(int i=0; i<maze->numDimensions; ++i) {
        validMoves[i*2] = i+1;
        validMoves[i*2+1] = -(i+1);
    }

    /* start at startPos */
    pos_list_clear(&maze->solution);
    if( maze_solve_dfs(maze, maze->startPos) ) {
        printf("Found solution with length %i\n", maze->solution.posListNum);
        return 1;
    }

    printf("No solution found!\n");
    return 0;
}


int maze_generate(maze_t *maze) {
    printf("Generating %iD maze.\n", maze->numDimensions);

    /* set starting point for generation */
    int *start = calloc(maze->numDimensions,sizeof(int));
    for(int i=0; i<maze->numDimensions; ++i) {
        start[i] = rand()%(maze->dimensions[i]-2)+1;
        start[i] |= 1;  // force initial coordinate to be all odd
    }

    /* recursively clear cells */
    maze_gen_step(maze,start);
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
        done = maze_solve(maze);
    }

    /* while not completely full (i.e., any uncleared 2x2 region exists) */
    int *restartPos = calloc(maze->numDimensions,sizeof(int));
    retries = 250;
    int restarts = 0;
    while(maze_unfinished(maze) && --retries) {
        maze_get_restart_location(maze, restartPos);

        printf("trying restart from position: ");
        for(int i=0; i<maze->numDimensions; ++i)
            printf("%i ", restartPos[i]);
        printf("\n");

        /* try using as an unreachable starting point */
        maze_gen_step(maze, restartPos);
        ++restarts;
    }
    free(start); start=NULL;

    return restarts+1;
}


int maze_write(maze_t *maze, char *filename) {

    /* open file */
    FILE *fp = fopen(filename,"w");
    if(fp==NULL ) {
        perror("fopen");
        return -1;
    }

    /* write dimensions */
    fprintf(fp, "%i %i\n", maze->numDimensions, maze->numFaces);
    for(int i=0; i<maze->numDimensions; ++i)
        fprintf(fp, "%i ", maze->dimensions[i]);
    fprintf(fp,"\n");

    /* write maze faces */
    for(int face=0; face<maze->numFaces; ++face) {
        int rows = maze->faces[face].rows;
        int cols = maze->faces[face].cols;
        int d1 = maze->faces[face].d1;
        int d2 = maze->faces[face].d2;
        fprintf(fp, "%i %i\n", cols, rows);
        for(int row=0; row<rows; ++row) {
            for(int col=0; col<cols; ++col) {
                int cell = face_get_cell(&maze->faces[face], col, row);
                int ch = cell?'1':'0';
                if( maze->startPos[d1] == col && maze->startPos[d2] == row )
                    ch = 'S';
                if( maze->endPos[d1] == col && maze->endPos[d2] == row )
                    ch = 'E';
                fprintf(fp, "%c ", ch);
            }
            fprintf(fp, "\n");
        }
    }

    /* write solution from start to end positions */
    fprintf(fp, "%i\n", maze->solution.posListNum);
    for(int m=0; m<maze->solution.posListNum; ++m) {
        for(int i=0; i<maze->numDimensions; ++i) {
            fprintf(fp, "%i ", maze->solution.positions[m][i]);
        }
        fprintf(fp, "\n");
    }

    /* close file */
    fclose(fp); fp=NULL;

    return 0;
}


int maze_load(maze_t *maze, char *filename) {
    /* open file */
    FILE *fp = fopen(filename,"r");
    if(fp==NULL ) {
        perror("fopen");
        return -1;
    }

    /* read dimensions */
    fscanf(fp, "%i %i\n", &maze->numDimensions, &maze->numFaces);
    int *sizes = calloc(maze->numDimensions, sizeof(int));
    for(int i=0; i<maze->numDimensions; ++i)
        fscanf(fp,"%i ", &sizes[i]);   

    /* initialize the maze */
    maze_init(maze, maze->numDimensions, sizes);
    free(sizes); sizes=NULL;

    /* read maze faces */
    for(int face=0; face<maze->numFaces; ++face) {
        int rows, cols;
        int d1 = maze->faces[face].d1;
        int d2 = maze->faces[face].d2;
        fscanf(fp, "%i %i\n", &cols, &rows);
        for(int row=0; row<rows; ++row) {
            for(int col=0; col<cols; ++col) {
                char cell;
                fscanf(fp, "%c ", &cell);
                if(cell=='S') {
                    /* is start location on face */
                    if( maze->startPos[d1]!=-1 && maze->startPos[d1]!=col
                        && maze->startPos[d2]!=-1 && maze->startPos[d2]!=row )
                        fprintf(stderr, "%s: Inconsistency on face %i of maze start location.\n\tstart[%i,%i]=%i,%i and %i,%i\n", __FUNCTION__, face, d1, d2, maze->startPos[d1], maze->startPos[d2], col, row);
                    maze->startPos[d1] = col;
                    maze->startPos[d2] = row;
                    cell='0';
                }
                if(cell=='E') {
                    /* is end location on face */
                    if( maze->endPos[d1]!=-1 && maze->endPos[d1]!=col
                        && maze->endPos[d2]!=-1 && maze->endPos[d2]!=row )
                        fprintf(stderr, "%s: Inconsistency on face %i of maze end location.\n\tend[%i,%i]=%i,%i and %i,%i\n", __FUNCTION__, face, d1, d2, maze->endPos[d1], maze->endPos[d2], col, row);
                    maze->endPos[d1] = col;
                    maze->endPos[d2] = row;
                    cell='0';
                }
                face_set_cell(&maze->faces[face], col, row, cell-'0');
            }
        }
    }

    /* read solution from start to end positions */
    int numPos = 0;
    position_t pos = calloc(maze->numDimensions,sizeof(int));
    fscanf(fp, "%i\n", &numPos);
    for(int m=0; m<numPos; ++m) {
        for(int i=0; i<maze->numDimensions; ++i) {
            fscanf(fp, "%i ", &pos[i]);
        }
        pos_list_push(&maze->solution, pos);
    }
    free(pos); pos=NULL;

    /* close file */
    fclose(fp); fp=NULL;

    return 0;
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
