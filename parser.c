// Shunting-Yard implementation & string tokenization

#include "parser.h"
#include <ctype.h>
#include <string.h>

/* Internal Stacks (encapsulated as module-private) */
typedef struct {
    Node* items[STACK_MAX];
    int top;
} NodeStack;

typedef struct {
    char items[STACK_MAX];
    int top;
} OpStack;

static void push_node(NodeStack* s, Node* n) { s->items[++(s->top)] = n; }
static Node* pop_node(NodeStack* s) { return s->items[(s->top)--]; }

static void push_op(OpStack* s, char c) { s->items[++(s->top)] = c; }
static char pop_op(OpStack* s) { return s->items[(s->top)--]; }
static char peek_op(OpStack* s) { return s->top >= 0 ? s->items[s->top] : '\0'; }

static int precedence(char op) {
    switch (op) {
        case '+': case '-': return 1;
        case '*': case '/': return 2;
        case '^':          return 3;
        case 'e': case 's': case 'c': return 4; /* exp, sin, cos */
        default:           return 0;
    }
}

static int is_right_associative(char op) {
    return op == '^';
}

static NodeType char_to_nodetype(char op) {
    switch (op) {
        case '+': return NODE_ADD;
        case '-': return NODE_SUB;
        case '*': return NODE_MUL;
        case '/': return NODE_DIV;
        case '^': return NODE_POW;
        case 'e': return NODE_EXP;
        default:  return NODE_CONST;
    }
}

static void process_operator(OpStack* op_stack, NodeStack* node_stack) {
    char op = pop_op(op_stack);
    if (op == 'e' || op == 's' || op == 'c') {
        Node* arg = pop_node(node_stack);
        NodeType t = (op == 'e') ? NODE_EXP : (op == 's') ? NODE_SIN : NODE_COS;
        push_node(node_stack, create_node(t, 0, 0, arg, NULL));
    } else {
        Node* right = pop_node(node_stack);
        Node* left = pop_node(node_stack);
        push_node(node_stack, create_node(char_to_nodetype(op), 0, 0, left, right));
    }
}

Node* parse_expression(const char* expr) {
    NodeStack node_stack = {.top = -1};
    OpStack op_stack = {.top = -1};
    int i = 0;

    while (expr[i] != '\0') {
        if (isspace(expr[i])) { i++; continue; }

        if (isdigit(expr[i]) || expr[i] == '.') {
            char* endptr;
            double val = strtod(&expr[i], &endptr);
            push_node(&node_stack, create_node(NODE_CONST, val, 0, NULL, NULL));
            i = endptr - expr;
            continue;
        }

        /* Parse Trigonometric & Exponential Functions */
        if (strncmp(&expr[i], "exp", 3) == 0) { push_op(&op_stack, 'e'); i += 3; continue; }
        if (strncmp(&expr[i], "sin", 3) == 0) { push_op(&op_stack, 's'); i += 3; continue; }
        if (strncmp(&expr[i], "cos", 3) == 0) { push_op(&op_stack, 'c'); i += 3; continue; }

        if (isalpha(expr[i])) {
            push_node(&node_stack, create_node(NODE_VAR, 0, expr[i], NULL, NULL));
            i++;
            continue;
        }

        if (expr[i] == '(') { push_op(&op_stack, '('); i++; continue; }

        if (expr[i] == ')') {
            while (op_stack.top >= 0 && peek_op(&op_stack) != '(') {
                process_operator(&op_stack, &node_stack);
            }
            if (op_stack.top >= 0) pop_op(&op_stack);
            i++;
            continue;
        }

        if (strchr("+-*/^", expr[i])) {
            char current_op = expr[i];
            while (op_stack.top >= 0 && peek_op(&op_stack) != '(' &&
                   (precedence(peek_op(&op_stack)) > precedence(current_op) ||
                   (precedence(peek_op(&op_stack)) == precedence(current_op) && !is_right_associative(current_op)))) {
                process_operator(&op_stack, &node_stack);
            }
            push_op(&op_stack, current_op);
            i++;
            continue;
        }

        i++;
    }

    while (op_stack.top >= 0) {
        process_operator(&op_stack, &node_stack);
    }

    return node_stack.top >= 0 ? pop_node(&node_stack) : NULL;
}