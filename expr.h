// AST structure definitions & core tree operations

#ifndef EXPR_H
#define EXPR_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Node types supported by the AST */
typedef enum {
    NODE_CONST,
    NODE_VAR,
    NODE_ADD,
    NODE_SUB,
    NODE_MUL,
    NODE_DIV,
    NODE_POW,
    NODE_EXP,
    NODE_SIN,
    NODE_COS
} NodeType;

/* Expression Tree (AST) Node Structure */
typedef struct Node {
    NodeType type;
    double val;            /* Value if type == NODE_CONST */
    char var_name;         /* Variable identifier if type == NODE_VAR */
    struct Node *left;     /* Left sub-tree / operand */
    struct Node *right;    /* Right sub-tree / operand */
} Node;

/* AST Memory Management */
Node* create_node(NodeType type, double val, char var_name, Node* left, Node* right);
void free_tree(Node* node);
Node* copy_tree(Node* node);

/* Expression Operations */
double evaluate(Node* node, double x_val);
double integrate_numerical(Node* root, double a, double b, int n);
Node* differentiate(Node* node, char var);
/* Symbolic integration w.r.t a variable (e.g., 'x') */
Node* integrate(Node* node, char var);\
/* Recursively simplifies an AST by removing identity operations and folding constants */
Node* simplify_tree(Node* node);
Node* simplify_full(Node* node);
void print_inorder(Node* node);

#endif /* EXPR_H */