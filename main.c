/*
 * main.c
 * PuzzleMaze: A 3D puzzle model generator
 *
 * Copyright (c) 2020 Bryan Franklin. All rights reserved.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "maze.h"

static void show_help(int argc, char **argv) {
    printf("%s [-d dims] {-g n,n,n | -i file} [-s] [-m file.stl] [-o file.txt]\n",
            argv[0]);
    exit(0);
}

int main(int argc, char **argv) {
    maze_t maze;
    int dims = 3;
    int *sizes=NULL;
    char *inputFile = NULL;
    char *outputFile = NULL;
    char *stlFile = NULL;
    char *stlSolFile = NULL;
    int genMaze = 0;
    int doSolve = 0;

    /* set defaults */
    dims = 3;
    sizes = calloc(dims,sizeof(int));
    sizes[0] = 11;
    sizes[1] = 11;
    sizes[2] = 11;
    genMaze = 1;

    char ch='\0';
    while( (ch=getopt(argc, argv, "d:g:hi:m:o:p:s"))!=-1 ) {
        switch(ch) {
            case 'd':
                dims = atoi(optarg);
                sizes = realloc(sizes,dims*sizeof(int));
                break;
            case 'g':
                /* generate maze */
                genMaze = 1;
                /* argument is dimensions */
                char *str = strtok(optarg,",");
                sizes[0] = atoi(str);
                for(int i=1; i<dims; ++i) {
                    sizes[i] = atoi(str);
                    str = strtok(NULL,",");
                }
                break;
            case 'h':
            case '?':
                show_help(argc,argv);
                break;
            case 'i':
                /* load maze from file given as argument */
                inputFile = strdup(optarg);
                genMaze = 0;
                break;
            case 'm':
                /* output STL model to file given as argument */
                stlFile = strdup(optarg);
                break;
            case 'o':
                /* output maze to file given as argument */
                outputFile = strdup(optarg);
                break;
            case 'p':
                /* output solution STL model to file given as argument */
                stlSolFile = strdup(optarg);
                break;
            case 's':
                /* generate a solution */
                doSolve = 1;
                break;
        }
    }

    /* check that sufficient settings exist */
    if( outputFile==NULL && stlFile==NULL) {
        fprintf(stderr,"\n\nNo output given, use -o and/or -m to specify an output filename.\n\n");
        show_help(argc,argv);
    }
    if( inputFile==NULL && !genMaze) {
        fprintf(stderr,"\n\nNo maze source given, use -i to specify an input filename\nor -g to generate a random maze.\n\n");
        show_help(argc,argv);
    }

    if( genMaze ) {
        maze_init(&maze, dims, sizes);
        printf("Generating maze.\n");
        maze_generate(&maze);
    }
    else if( inputFile ) {
        printf("Loading maze from '%s'.\n", inputFile);
        maze_load(&maze,inputFile);
    }

    if( doSolve )
        maze_solve(&maze);

    if( outputFile ) {
        printf("Writing %iD maze to `%s`.\n", maze.numDimensions, outputFile);
        maze_write(&maze, outputFile);
    }
    if( stlFile ) {
        printf("Exporting %iD maze to STL file `%s`.\n", maze.numDimensions, stlFile);
        maze_export_stl(&maze, stlFile);
    }
    if( stlSolFile ) {
        printf("Exporting %iD maze solution to STL file `%s`.\n", maze.numDimensions, stlFile);
        maze_export_stl_solution(&maze, stlSolFile);
    }

    maze_free(&maze);

    return 0;
}
