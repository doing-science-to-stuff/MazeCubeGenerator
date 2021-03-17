/*
 * maze.c
 * MazeCubeGen: maze cube generator
 *
 * Copyright (c) 2020-2021 Bryan Franklin. All rights reserved.
 */
#include <ctype.h>
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


static int face_reset(face_t *face) {
    memset(face->cells, 1, face->rows*face->cols);

    return 0;
}


int face_free(face_t *face) {
    free(face->cells); face->cells = NULL;
    memset(face, '\0', sizeof(*face));

    return 0;
}

static int face_set_cell(face_t *face, int row, int col, int value) {
    if( row < 0 || row>=face->rows
        || col < 0 || col >=face->cols ) {
        fprintf(stderr, "Request to set out of bounds face cell %i,%i (size: %ix%i)\n",
            row, col, face->rows, face->cols);
        return 0;
    }
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
    list->capacity = 10;
    list->num = 0;
    list->positions = calloc(list->capacity,sizeof(position_t));
    return (list->positions!=NULL);
}


static int pos_list_free(position_list_t *list) {
    for(int i=0; i<list->num; ++i) {
        free(list->positions[i]); list->positions[i] = NULL;
    }
    free(list->positions); list->positions=NULL;
    memset(list,'\0',sizeof(*list));
    return 0;
}


static int pos_list_push(position_list_t *list, position_t pos) {
    /* add to list of restart locations */
    if( list->num == list->capacity ) {
        int newCap = list->capacity*2+1;
        void *tmp = realloc(list->positions,newCap*sizeof(int*));
        if( tmp==NULL ) {
            perror("realloc");
            return 0;
        }
        list->positions = tmp;
        list->capacity = newCap;
    }
    list->positions[list->num] = calloc(list->numDimensions,sizeof(int));
    memcpy(list->positions[list->num], pos, list->numDimensions*sizeof(int));
    ++list->num;
    return 1;
}


static int pos_list_pop(position_list_t *list, position_t pos) {
    if( list->num <= 0 )
        return 0;
    if( pos!=NULL )
        memcpy(pos,list->positions[list->num-1],list->numDimensions*sizeof(int));
    if( list->num > 0 ) {
        free(list->positions[list->num-1]); list->positions[list->num-1]=NULL;
        --list->num;
    }
    return 1;
}


static int pos_list_clear(position_list_t *list) {
    while(pos_list_pop(list,NULL));
    return 0;
}


static int pos_list_random(position_list_t *list, position_t pos) {
    if( list->num <= 0 )
        return 0;
    int which = rand()%list->num;
    memcpy(pos, list->positions[which], list->numDimensions*sizeof(int));
    return 1;
}


static int pos_list_rfind(position_list_t *list, position_t pos) {
    for(int i=list->num-1; i>=0; --i) {
        if( memcmp(list->positions[i], pos, list->numDimensions*sizeof(int)) == 0 )
        return i;
    }
    return -1;
}


static int position_increment(maze_t *maze, position_t pos) {
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


static int maze_parse_options(maze_t *maze, char *options) {

    if( !options )
        return 0;

    maze->pathSelMode = tolower(options[0]);

    return 1;
}


int maze_init_str(maze_t *maze, char *cfgStr) {

    int numDimensions = 3;
    int *sizes = NULL;

    /* copy config string before destructly parsing it */
    char *cfg = strdup(cfgStr);
    char *options = NULL;
    if( strchr(cfg, ':') ) {
        /* config is more than just dimensions */
        /* skip dimensions string */
        strtok(cfg,":");
        char *optionStr = strtok(NULL,":");
        if( optionStr )
            options = strdup(optionStr);
    }

    /* count separators in dimStr */
    char *dimStr = strdup(cfg);
    char *curr = dimStr;
    int numSep = 0;
    while(*curr) {
        if( *curr=='x' || *curr==',' )
            ++numSep;
        ++curr;
    }
    numDimensions = numSep+1;

    /* allocate sizes array */
    sizes = calloc(numDimensions, sizeof(int));

    /* extract dimensions form dimStr */
    char *tok = strtok(dimStr,",x");
    sizes[0] = atoi(tok);
    for(int i=1; i<numDimensions; ++i) {
        sizes[i] = atoi(tok);
        tok = strtok(NULL,",x");
    }

    #if 1
    printf("%s -> %i\t%i", cfg, numDimensions, sizes[0]);
    for(int i=1; i<numDimensions; ++i) {
        printf(",%i", sizes[i]);
    }
    printf("\n");
    #endif /* 0 */

    /* do actual initialization */
    int ret = maze_init(maze, numDimensions, sizes, options);

    /* free buffers and extra copies of strings */
    if( options ) {
        free(options); options=NULL;
    }
    free(sizes); sizes=NULL;
    free(dimStr); dimStr=NULL;
    free(cfg); cfg=NULL;

    return ret;
}


int maze_init(maze_t *maze, int numDimensions, int *sizes, char *options) {

    /* initialize maze structure */
    printf("Initializing %iD maze.\n", numDimensions);
    memset(maze,'\0',sizeof(*maze));
    maze->numDimensions = numDimensions;
    maze->numFaces = (numDimensions*(numDimensions-1))/2;
    maze->dimensions = calloc(numDimensions,sizeof(int));
    maze->maxSegments = -1;
    maze->minPathLength = -1;
    for(int i=0; i<numDimensions; ++i) {
        if( sizes[i] < 3 ) {
            fprintf(stderr, "Warning: edge sizes must be at least 3, adjusting dimension %i from %i to %i.\n", i, sizes[i], 3);
            sizes[i] = 3;
        }
        if( (sizes[i]%2) == 0 ) {
            fprintf(stderr, "Warning: edge sizes must be odd, adjusting dimension %i from %i to %i.\n", i, sizes[i], sizes[i]+1);
            sizes[i] += 1;
        }
        maze->dimensions[i] = sizes[i];
    }
    if( options )
        maze_parse_options(maze, options);

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
            if( numDimensions > 3 )
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
    pos_list_free(&maze->reachable);
    pos_list_free(&maze->solution);
    free(maze->startPos); maze->startPos=NULL;
    free(maze->endPos); maze->endPos=NULL;
    free(maze->faces); maze->faces=NULL;
    free(maze->dimensions); maze->dimensions=NULL;
    memset(maze,'\0', sizeof(*maze));

    return 0;
}


void maze_set_segments(maze_t *maze, int segs) {
    maze->maxSegments = segs;
}


void maze_set_path_length(maze_t *maze, int len) {
    maze->minPathLength = len;
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

    /* check if those positions are both cleared but not in the reachable set */
    if( maze_position_clear(maze, midPos)
            && maze_position_clear(maze, nextPos)
            && pos_list_rfind(&maze->reachable, nextPos) < 0 ) {
        /* path already exists due to overlap with previous path building */
        free(midPos); midPos=NULL;
        free(nextPos); nextPos=NULL;
        return 1;
    }

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
    free(midPos); midPos=NULL;
    free(nextPos); nextPos=NULL;

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
        if( pos[i]<maze->dimensions[i]
            && maze_allow_clear(maze,pos,i+1) ) {
            validMoves[numMoves++] = i+1;
        }
        if( pos[i]>0
            && maze_allow_clear(maze,pos,-(i+1)) ) {
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
        if( numMoves > 1 ) {
            maze_gen_step(maze,pos);
        }
    }
    free(validMoves); validMoves=NULL;

    return 0;
}


static int maze_unfinished(maze_t *maze) {

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


static int maze_get_restart_location(maze_t *maze, int *pos) {

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


static size_t maze_cell_degree(maze_t *maze, position_t pos) {

    if( !maze_position_clear(maze, pos) )
        return 0;

    position_t neighbor;
    neighbor = calloc(maze->numDimensions,sizeof(*neighbor));
    size_t degree = 0;
    for(int i=0; i<maze->numDimensions; ++i) {
        memcpy(neighbor, pos, sizeof(*neighbor)*maze->numDimensions);
        ++neighbor[i];
        if( maze_position_clear(maze, neighbor) )
            ++degree;

        memcpy(neighbor, pos, sizeof(*neighbor)*maze->numDimensions);
        --neighbor[i];
        if( maze_position_clear(maze, neighbor) )
            ++degree;
    }
    free(neighbor); neighbor=NULL;

    return degree;
}


static int maze_pick_goals_random(maze_t *maze) {

    if( maze->reachable.num < 2 )
        return 0;

    size_t posSize = sizeof(*maze->startPos) * maze->numDimensions;
    int numRetries = 100;
    do {
        pos_list_random(&maze->reachable, maze->startPos);
        pos_list_random(&maze->reachable, maze->endPos);
    } while( !memcmp(maze->startPos, maze->endPos, posSize) && --numRetries>0);

    return 1;
}


static int* validMoves = NULL;
static int maze_find_path(maze_t *maze, position_t pos, position_list_t* path, position_t goal) {

    pos_list_push(path, pos);

    /* check for endPos */
    if( memcmp(goal, pos, maze->numDimensions*sizeof(int))==0 ) {
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
        if( pos_list_rfind(path, newPos) >= 0 ) {
            continue;
        }

        /* recurse */
        if( maze_find_path(maze, newPos, path, goal) ) {
            free(newPos); newPos=NULL;
            return 1;
        }
    }

    pos_list_pop(path, NULL);
    free(newPos); newPos=NULL;

    return 0;
}


static int maze_pick_goals_optimal(maze_t *maze, char mode) {

    /* TODO Add parameter to allow selecting of solution metric */

    /* build list of all dead ends */
    position_list_t dead_ends;
    pos_list_init(&dead_ends, maze->numDimensions);
    position_t pos = NULL;
    pos = malloc(maze->numDimensions * sizeof(*pos));
    for(int i = 0; i<maze->numDimensions; ++i) {
        pos[i] = 1;
    }

    int done = 0;
    while(!done) {
        size_t degree = maze_cell_degree(maze, pos);

        /* number of dead ends */
        if( degree == 1 )
            pos_list_push(&dead_ends, pos);

        /* update pos */
        done = position_increment(maze, pos);
    }
    printf("%s: dead ends: %i\n", __FUNCTION__, dead_ends.num);

    /* initialize validMoves */
    validMoves = calloc(2*maze->numDimensions, sizeof(int));
    for(int i=0; i<maze->numDimensions; ++i) {
        validMoves[i*2] = i+1;
        validMoves[i*2+1] = -(i+1);
    }

    /* for all pairs of dead ends */
    position_list_t path;
    pos_list_init(&path, maze->numDimensions);
    int best_value = 0;
    for(int i=0; i<dead_ends.num; ++i) {
        for(int j=i+1; j<dead_ends.num; ++j) {
            /* find solution, if possible */
            pos_list_clear(&path);
            maze_find_path(maze, dead_ends.positions[i], &path, dead_ends.positions[j]);
            if( path.num == 0 )
                continue;

            /* compute metrics for solution */
            int value=0;
            if( mode == 'l' || mode == 'o' ) {
                value = path.num;
            } else if( mode == 'b' ) {
                value = 0;
                for(int k=0; k<path.num; ++k) {
                    /* count branch points */
                    if( maze_cell_degree(maze, path.positions[k]) > 2 )
                        ++value;
                }
            } else {
                fprintf(stderr, "Unrecognized path metric mode '%c'.\n", mode);
            }

            /* pick endpoints, if better than current best */
            if( value > best_value ) {
                printf("  new best path metric value (%i).\n", value);
                if( maze->startPos==NULL )
                    maze->startPos = calloc(maze->numDimensions, sizeof(*maze->startPos));
                if( maze->endPos==NULL )
                    maze->endPos = calloc(maze->numDimensions, sizeof(*maze->endPos));

                /* randomly assign start/end to the dead ends */
                if( rand()%2 ) {
                    memcpy(maze->startPos, dead_ends.positions[i],
                            maze->numDimensions*sizeof(*maze->startPos));
                    memcpy(maze->endPos, dead_ends.positions[j],
                            maze->numDimensions*sizeof(*maze->endPos));
                } else {
                    memcpy(maze->startPos, dead_ends.positions[j],
                            maze->numDimensions*sizeof(*maze->startPos));
                    memcpy(maze->endPos, dead_ends.positions[i],
                            maze->numDimensions*sizeof(*maze->endPos));
                }

                best_value = value;
            }
        }
    }
    pos_list_free(&path);
    pos_list_free(&dead_ends);
    free(pos); pos=NULL;
    free(validMoves); validMoves=NULL;

    return 1;
}


int maze_pick_goals(maze_t *maze) {
    if( maze->pathSelMode=='o'
        || maze->pathSelMode=='l'
        || maze->pathSelMode=='b' )
        return maze_pick_goals_optimal(maze, maze->pathSelMode);

    return maze_pick_goals_random(maze);
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
    if( maze_find_path(maze, maze->startPos, &maze->solution, maze->endPos) ) {
        printf("Found solution with length %i\n", maze->solution.num);
        free(validMoves); validMoves=NULL;
        return 1;
    }
    free(validMoves); validMoves=NULL;

    printf("No solution found!\n");
    return 0;
}


static void maze_reset_faces(maze_t *maze) {
    for(int i=0; i<maze->numFaces; ++i) {
        face_reset(&maze->faces[i]);
    }
}


int maze_generate(maze_t *maze) {
    printf("Generating %iD maze.\n", maze->numDimensions);

    int restarts = 0;
    do {
        maze_reset_faces(maze);
        pos_list_clear(&maze->reachable);

        /* set starting point for generation */
        int *start = calloc(maze->numDimensions,sizeof(int));
        for(int i=0; i<maze->numDimensions; ++i) {
            start[i] = rand()%(maze->dimensions[i]-2)+1;
            start[i] |= 1;  // force initial coordinate to be all odd
        }

        /* recursively clear cells */
        maze_gen_step(maze,start);
        free(start); start=NULL;

        /* pick start and end positions */
        int done = 0;
        int retries = 100;
        do {
            /* pick start and end locations */
            maze_pick_goals(maze);

            printf("\tStart: ");
            for(int i=0; i<maze->numDimensions; ++i) {
                printf("%i ", maze->startPos[i]);
            }
            printf("\n\tEnd:   ");
            for(int i=0; i<maze->numDimensions; ++i) {
                printf("%i ", maze->endPos[i]);
            }
            printf("\n");

            /* find and record solution */
            done = maze_solve(maze);
            if( done )
                printf("\tsolution length: %i\n", maze->solution.num);
        } while( --retries
            && (!done || maze->solution.num < maze->minPathLength) );

        /* while not completely full (i.e., any uncleared 2x2 region exists) */
        restarts = 0;
        int *restartPos = calloc(maze->numDimensions,sizeof(int));
        retries = 10000;
        if( maze->maxSegments > 0 )
            retries  = maze->maxSegments-1;
        while(maze_unfinished(maze) && --retries) {
            maze_get_restart_location(maze, restartPos);

            printf("trying restart from position: ");
            for(int i=0; i<maze->numDimensions; ++i)
                printf("%i ", restartPos[i]);
            printf("\n");
            pos_list_clear(&maze->reachable);

            /* try using as an unreachable starting point */
            maze_gen_step(maze, restartPos);
            ++restarts;
        }
        free(restartPos); restartPos=NULL;
    } while( maze->maxSegments>0
        && (restarts+1)>maze->maxSegments );
    printf("\t%i segments\n", restarts+1);

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
    fprintf(fp, "%i\n", maze->solution.num);
    for(int m=0; m<maze->solution.num; ++m) {
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
    maze_init(maze, maze->numDimensions, sizes, NULL);
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

static int maze_cell_trivial(maze_t *maze, position_t pos) {

    if( !maze_position_clear(maze, pos) )
        return 0;

    /* return 1 if moves from cell are only possible in opposite directions
     * along a single axis */
    position_t posDir = NULL, negDir = NULL;
    posDir = malloc(maze->numDimensions * sizeof(*posDir));
    negDir = malloc(maze->numDimensions * sizeof(*negDir));
    int num = 0;
    for(int d = 0; d<maze->numDimensions; ++d) {
        /* copy current pos */
        memcpy(posDir, pos, maze->numDimensions * sizeof(*posDir));
        memcpy(negDir, pos, maze->numDimensions * sizeof(*negDir));

        /* update along axis d */
        ++posDir[d];
        --negDir[d];

        /* check neighboring positions */
        int posClear = maze_position_clear(maze, posDir);
        int negClear = maze_position_clear(maze, negDir);
        if( (posClear && !negClear) || (!posClear && negClear) ) {
            /* can only move one direction along axis d */
            free(posDir); posDir=NULL;
            free(negDir); negDir=NULL;
            return 0;
        }
        if( posClear && negClear )
            ++num;
    }
    free(posDir); posDir=NULL;
    free(negDir); negDir=NULL;

    if( num==1 ) {
        return 1;
    }

    return 0;
}

int maze_export_gv(maze_t *maze, char *filename) {
    /* open file */
    FILE *fp = fopen(filename,"w");
    if(fp==NULL ) {
        perror("fopen");
        return -1;
    }

    /* start position counter at all 1s */
    position_t pos = NULL;
    pos = malloc(maze->numDimensions * sizeof(*pos));
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

    fprintf(fp, "graph G {\nrankdir=TB;\n");
    fprintf(fp, "node [shape = box, color = red]; \"%i", maze->startPos[0]);
    for(int i=1; i<maze->numDimensions; ++i)
        fprintf(fp, ",%i", maze->startPos[i]);
    fprintf(fp, "\" \"%i", maze->endPos[0]);
    for(int i=1; i<maze->numDimensions; ++i)
        fprintf(fp, ",%i", maze->endPos[i]);
	fprintf(fp, "\";\n");

    #if 0
    /* highlight reachable set */
    fprintf(fp, "node [shape = circle, color = green];");
    for(int j=0; j<maze->reachable.num; ++j) {
        fprintf(fp, " \"%i", maze->reachable.positions[j][0]);
        for(int i=1; i<maze->numDimensions; ++i)
            fprintf(fp, ",%i", maze->reachable.positions[j][i]);
        fprintf(fp, "\"");
    }
	fprintf(fp, ";\n");
    #endif /* 0 */

    fprintf(fp, "node [shape = circle, color = black];\n");

    /* iterate through all possible positions */
    position_t neighbor = NULL;
    neighbor = malloc(maze->numDimensions * sizeof(*neighbor));
    int done = 0;
    while( !done ) {

        /* for each non-trivial open cell */
        if( maze_position_clear(maze, pos)
            && !maze_cell_trivial(maze, pos) ) {

            /* check neighbor in each direction */
            for(int m = 0; m<2*maze->numDimensions; ++m) {

                /* compute neighbor position */
                memcpy(neighbor, pos, sizeof(*pos)*maze->numDimensions);
                int move = moves[m];

                int length = 0;
                do {
                    if( move >= 0 )
                        neighbor[move-1] += 1;
                    else
                        neighbor[-move-1] -= 1;
                    ++length;
                } while( maze_cell_trivial(maze, neighbor) );

                /* export edge to each accessible neighbor */
                if( length>0
                    && memcmp(pos, neighbor, sizeof(*pos)*maze->numDimensions) < 0
                    && maze_position_clear(maze, neighbor) ) {
                    fprintf(fp, "\"%i", pos[0]);
                    for(int i=1; i<maze->numDimensions; ++i) {
                        fprintf(fp, ",%i", pos[i]);
                    }
                    fprintf(fp, "\" -- \"%i", neighbor[0]);
                    for(int i=1; i<maze->numDimensions; ++i) {
                        fprintf(fp, ",%i", neighbor[i]);
                    }
                    fprintf(fp, "\" [ label = \"%i\" ];\n", length);
                }
            }
        }

        /* update pos */
        done = position_increment(maze, pos);
    }

    fprintf(fp, "}\n");

    free(neighbor); neighbor=NULL;
    free(pos); pos=NULL;
    free(moves); moves=NULL;

    /* close file */
    fclose(fp); fp=NULL;

    return 0;
}


int maze_metrics(maze_t *maze) {

    /* some of these are based on:
     * https://puzzling.stackexchange.com/a/5922
     */

    /* solution metrics */
    if( maze->solution.num > 0 ) {
        /* simple solution length */
        printf("solution length: %i\n", maze->solution.num);

        int sol_dead_ends = 0;
        int sol_branches = 0;
        for(int i=0; i<maze->solution.num; ++i) {
            size_t degree = maze_cell_degree(maze, maze->solution.positions[i]);

            /* number of dead ends */
            if( degree == 1 )
                ++sol_dead_ends;

            /* count decision points along solution */
            if( degree > 2 )
                ++sol_branches;
        }
        printf("dead ends: %i (solution)\n", sol_dead_ends);
        printf("branch points: %i (solution)\n", sol_branches);
    }

    /* entire cube metrics */
    /* start position counter at all 1s */
    position_t pos = NULL;
    pos = malloc(maze->numDimensions * sizeof(*pos));
    for(int i = 0; i<maze->numDimensions; ++i) {
        pos[i] = 1;
    }

    int done = 0;
    int dead_ends = 0;
    int branches = 0;
    while(!done) {
        size_t degree = maze_cell_degree(maze, pos);

        /* number of dead ends */
        if( degree == 1 )
            ++dead_ends;

        /* count decision points along solution */
        if( degree > 2 )
            ++branches;

        /* update pos */
        done = position_increment(maze, pos);
    }
    printf("dead ends: %i\n", dead_ends);
    printf("branch points: %i\n", branches);

    /* surface cell reuse */

    /* min/max extents of solution for each dimension */


    return 0;
}
