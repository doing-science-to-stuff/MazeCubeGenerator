/*
 * bst.h
 * MazeCubeGen: maze cube generator
 *
 * Copyright (c) 2026 Bryan Franklin. All rights reserved.
 */
#ifndef BST_H
#define BST_H

typedef struct bst_node {
    struct bst_node *left, *right;
    void *value;
} bst_node_t;

typedef struct bst {
    bst_node_t *root;
    int count;
    int (*compar)(const void*, const void*);
} bst_t;

static int bst_init(bst_t *tree) {
    if( tree==NULL ) {
        return -1;
    }

    tree->root = NULL;
    tree->count = 0;
    tree->compar = NULL;

    return 0;
}


static int bst_set_comparitor(bst_t *tree, int (*compar)(const void*, const void*)) {
    tree->compar = compar;
    return 0;
}

static int bst_clear_node(bst_node_t *node) {
    if( node==NULL ) return 0;

    // post-order traversal, freeing along the way
    bst_clear_node(node->left); node->left=NULL;
    bst_clear_node(node->right); node->right=NULL;
    free(node->value); node->value=NULL;
    free(node); node=NULL;

    return 0;
}

static int bst_clear(bst_t *tree) {
    bst_clear_node(tree->root); tree->root=NULL;
    return 0;
}


static int bst_init_node(bst_node_t *node) {
    memset(node, '\0', sizeof(*node));
    return 0;
}


static int bst_insert_node(bst_t *tree, bst_node_t *node, void *value) {
    if( node == NULL ) { return -1; }

    if( node->value == NULL ) {
        node->value = value;
        return 1;
    }
    int cmp = tree->compar(value, node->value);
    if( cmp < 0 && node->left != NULL ) {
        return bst_insert_node(tree, node->left, value);
    } else if( cmp >= 0 && node->right != NULL ) {
        return bst_insert_node(tree, node->right, value);
    }

    bst_node_t *child=calloc(1, sizeof(bst_node_t));
    bst_init_node(child);
    if( cmp < 0  ) {
        node->left = child;
        bst_insert_node(tree, child, value);
    } else if( cmp >= 0  ) {
        node->right = child;
        bst_insert_node(tree, child, value);
    } else {
        return -1;
    }

    return 1;
}

static int bst_insert(bst_t *tree, void *value) {
     if( tree->root==NULL ) {
        tree->root = calloc(1, sizeof(bst_node_t));
        bst_init_node(tree->root);
        #if 0
        if( tree->root == NULL )    return 0;
        tree->root->value = value;
        tree->count = 1;
        return 1;
        #endif
     }

     if( bst_insert_node(tree, tree->root, value) > 0 )
        ++tree->count;
    
    return tree->count;
}

static bst_node_t *bst_find(bst_t *tree, void *value) {
    //printf("%s: %p (%s)\n", __FUNCTION__, tree, (char*)value);
    bst_node_t *curr = tree->root;
    while(curr!=NULL) {
        int cmp = tree->compar(value, curr->value);
        if( cmp == 0 )
            return curr;
        //if(value != NULL ) printf("\t%s\n", (char*)curr->value);

        if( cmp < 0 )
            curr = curr->left;
        else
            curr = curr->right;
    }

    return NULL;
}


#if 0
static void bst_print_node(bst_node_t *node) {
    if(node==NULL) return;
    if(node->value!=NULL)
        printf("%s\n", (char*)node->value);
    bst_print_node(node->left);
    bst_print_node(node->right);
}


static void bst_print(bst_t *tree) {
    bst_print_node(tree->root);
}
#endif /* 0 */

#endif /* BST_H */
