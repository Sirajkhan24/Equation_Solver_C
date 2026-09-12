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
        case 'e':          return 4; /* exp function */
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
    if (op == 'e') {
        Node* arg = pop_node(node_stack);
        push_node(node_stack, create_node(NODE_EXP, 0, 0, arg, NULL));
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
        if (isspace(expr[i])) {
            i++;
            continue;
        }

        /* Parse Numbers (Constants) */
        if (isdigit(expr[i]) || expr[i] == '.') {
            char* endptr;
            double val = strtod(&expr[i], &endptr);
            push_node(&node_stack, create_node(NODE_CONST, val, 0, NULL, NULL));
            i = endptr - expr;
            continue;
        }

        /* Parse 'exp(' function */
        if (strncmp(&expr[i], "exp", 3) == 0) {
            push_op(&op_stack, 'e');
            i += 3;
            continue;
        }

        /* Parse Variables (e.g. 'x') */
        if (isalpha(expr[i])) {
            push_node(&node_stack, create_node(NODE_VAR, 0, expr[i], NULL, NULL));
            i++;
            continue;
        }

        /* Parentheses */
        if (expr[i] == '(') {
            push_op(&op_stack, '(');
            i++;
            continue;
        }

        if (expr[i] == ')') {
            while (op_stack.top >= 0 && peek_op(&op_stack) != '(') {
                process_operator(&op_stack, &node_stack);
            }
            if (op_stack.top >= 0) pop_op(&op_stack); /* Pop '(' */
            i++;
            continue;
        }

        /* Operators */
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