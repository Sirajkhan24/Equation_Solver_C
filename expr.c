// Differentiation, integration, evaluation, printing
#include <stdbool.h>
#include <stdlib.h>
#include "expr.h"

Node* create_node(NodeType type, double val, char var_name, Node* left, Node* right) {
    Node* n = (Node*)malloc(sizeof(Node));
    if (!n) {
        fprintf(stderr, "Error: Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    n->type = type;
    n->val = val;
    n->var_name = var_name;
    n->left = left;
    n->right = right;
    return n;
}

void free_tree(Node* node) {
    if (!node) return;
    free_tree(node->left);
    free_tree(node->right);
    free(node);
}

Node* copy_tree(Node* node) {
    if (!node) return NULL;
    return create_node(node->type, node->val, node->var_name,
                       copy_tree(node->left), copy_tree(node->right));
}

double evaluate(Node* node, double x_val) {
    if (!node) return 0.0;
    switch (node->type) {
        case NODE_CONST: return node->val;
        case NODE_VAR:   return x_val;
        case NODE_ADD:   return evaluate(node->left, x_val) + evaluate(node->right, x_val);
        case NODE_SUB:   return evaluate(node->left, x_val) - evaluate(node->right, x_val);
        case NODE_MUL:   return evaluate(node->left, x_val) * evaluate(node->right, x_val);
        case NODE_DIV:   return evaluate(node->left, x_val) / evaluate(node->right, x_val);
        case NODE_POW:   return pow(evaluate(node->left, x_val), evaluate(node->right, x_val));
        case NODE_EXP:   return exp(evaluate(node->left, x_val));
        default:         return 0.0;
    }
}

/* Definite Numerical Integration using Simpson's 1/3 Rule */
double integrate_numerical(Node* root, double a, double b, int n) {
    if (n % 2 != 0) n++; // Simpson's rule requires an even number of intervals
    
    double h = (b - a) / n;
    double sum = evaluate(root, a) + evaluate(root, b);

    for (int i = 1; i < n; i++) {
        double x = a + i * h;
        if (i % 2 == 0) {
            sum += 2.0 * evaluate(root, x);
        } else {
            sum += 4.0 * evaluate(root, x);
        }
    }

    return (h / 3.0) * sum;
}

Node* differentiate(Node* node, char var) {
    if (!node) return NULL;

    switch (node->type) {
        case NODE_CONST:
            return create_node(NODE_CONST, 0.0, 0, NULL, NULL);

        case NODE_VAR:
            return create_node(NODE_CONST, (node->var_name == var) ? 1.0 : 0.0, 0, NULL, NULL);

        case NODE_ADD:
            return create_node(NODE_ADD, 0, 0,
                               differentiate(node->left, var),
                               differentiate(node->right, var));

        case NODE_SUB:
            return create_node(NODE_SUB, 0, 0,
                               differentiate(node->left, var),
                               differentiate(node->right, var));

        case NODE_MUL:
            /* Product Rule: (u * v)' = u' * v + u * v' */
            return create_node(NODE_ADD, 0, 0,
                               create_node(NODE_MUL, 0, 0, differentiate(node->left, var), copy_tree(node->right)),
                               create_node(NODE_MUL, 0, 0, copy_tree(node->left), differentiate(node->right, var)));

        case NODE_DIV:
            /* Quotient Rule: (u / v)' = (u' * v - u * v') / (v^2) */
            {
                Node* num = create_node(NODE_SUB, 0, 0,
                                        create_node(NODE_MUL, 0, 0, differentiate(node->left, var), copy_tree(node->right)),
                                        create_node(NODE_MUL, 0, 0, copy_tree(node->left), differentiate(node->right, var)));
                Node* den = create_node(NODE_POW, 0, 0, copy_tree(node->right), create_node(NODE_CONST, 2.0, 0, NULL, NULL));
                return create_node(NODE_DIV, 0, 0, num, den);
            }

        case NODE_EXP:
            /* Chain Rule: (exp(u))' = exp(u) * u' */
            return create_node(NODE_MUL, 0, 0,
                               create_node(NODE_EXP, 0, 0, copy_tree(node->left), NULL),
                               differentiate(node->left, var));

        case NODE_POW:
            /* Power Rule for polynomial terms (u^c): (u^n)' = n * u^(n-1) * u' */
            if (node->right && node->right->type == NODE_CONST) {
                double n = node->right->val;
                Node* pow_term = create_node(NODE_POW, 0, 0,
                                             copy_tree(node->left),
                                             create_node(NODE_CONST, n - 1.0, 0, NULL, NULL));
                Node* coeff_term = create_node(NODE_MUL, 0, 0,
                                               create_node(NODE_CONST, n, 0, NULL, NULL),
                                               pow_term);
                return create_node(NODE_MUL, 0, 0, coeff_term, differentiate(node->left, var));
            }
            return create_node(NODE_CONST, 0.0, 0, NULL, NULL);

        default:
            return NULL;
    }
}

Node* integrate(Node* node, char var) {
    if (!node) return NULL;

    switch (node->type) {
        case NODE_CONST:
            // Rule: ∫ c dx = c * x
            if (node->val == 0.0) {
                return create_node(NODE_CONST, 0.0, 0, NULL, NULL);
            }
            return create_node(NODE_MUL, 0, 0,
                               copy_tree(node),
                               create_node(NODE_VAR, 0, var, NULL, NULL));

        case NODE_VAR:
            // Rule: ∫ x dx = (1/2) * x^2
            if (node->var_name == var) {
                Node* pow_term = create_node(NODE_POW, 0, 0,
                                             copy_tree(node),
                                             create_node(NODE_CONST, 2.0, 0, NULL, NULL));
                return create_node(NODE_MUL, 0, 0,
                                   create_node(NODE_CONST, 0.5, 0, NULL, NULL),
                                   pow_term);
            }
            // If it's a different variable treat it as a constant: ∫ y dx = y * x
            return create_node(NODE_MUL, 0, 0,
                               copy_tree(node),
                               create_node(NODE_VAR, 0, var, NULL, NULL));

        case NODE_ADD:
            // Sum Rule: ∫ (f + g) dx = ∫ f dx + ∫ g dx
            return create_node(NODE_ADD, 0, 0,
                               integrate(node->left, var),
                               integrate(node->right, var));

        case NODE_SUB:
            // Difference Rule: ∫ (f - g) dx = ∫ f dx - ∫ g dx
            return create_node(NODE_SUB, 0, 0,
                               integrate(node->left, var),
                               integrate(node->right, var));

        case NODE_MUL:
            // Rule: ∫ (c * f) dx = c * ∫ f dx  (constant factor rule)
            if (node->left && node->left->type == NODE_CONST) {
                return create_node(NODE_MUL, 0, 0,
                                   copy_tree(node->left),
                                   integrate(node->right, var));
            }
            if (node->right && node->right->type == NODE_CONST) {
                return create_node(NODE_MUL, 0, 0,
                                   integrate(node->left, var),
                                   copy_tree(node->right));
            }
            fprintf(stderr, "Warning: Symbolic product integration (Integration by Parts) not fully supported.\n");
            return NULL;

        case NODE_POW:
            // Power Rule: ∫ (x^n) dx = (x^(n+1)) / (n+1)   [for n != -1]
            if (node->left && node->left->type == NODE_VAR && node->left->var_name == var &&
                node->right && node->right->type == NODE_CONST) {
                
                double n = node->right->val;
                if (n == -1.0) {
                    fprintf(stderr, "Error: ∫ x^-1 dx requires ln(x) support.\n");
                    return NULL;
                }
                
                double new_power = n + 1.0;
                Node* pow_term = create_node(NODE_POW, 0, 0,
                                             copy_tree(node->left),
                                             create_node(NODE_CONST, new_power, 0, NULL, NULL));
                
                return create_node(NODE_DIV, 0, 0,
                                   pow_term,
                                   create_node(NODE_CONST, new_power, 0, NULL, NULL));
            }
            return NULL;

        case NODE_EXP:
            // Rule: ∫ exp(x) dx = exp(x)
            if (node->left && node->left->type == NODE_VAR && node->left->var_name == var) {
                return copy_tree(node);
            }
            // Rule: ∫ exp(c * x) dx = (1/c) * exp(c * x)
            if (node->left && node->left->type == NODE_MUL) {
                Node* mul = node->left;
                if (mul->left && mul->left->type == NODE_CONST &&
                    mul->right && mul->right->type == NODE_VAR && mul->right->var_name == var) {
                    
                    double c = mul->left->val;
                    return create_node(NODE_MUL, 0, 0,
                                       create_node(NODE_CONST, 1.0 / c, 0, NULL, NULL),
                                       copy_tree(node));
                }
            }
            return NULL;

        default:
            return NULL;
    }
}


// ============================================================================
// 1. MATCHING HELPERS
// ============================================================================

static inline bool is_var(const Node* n, char var) {
    return n && n->type == NODE_VAR && (var == 0 || n->var_name == var);
}

static inline bool is_const(const Node* n, double val, bool check_val) {
    if (!n || n->type != NODE_CONST) return false;
    return check_val ? (n->val == val) : true;
}

static inline bool is_pow(const Node* n, char var, double exponent, bool check_exp) {
    if (!n || n->type != NODE_POW) return false;
    if (!is_var(n->left, var)) return false;
    return is_const(n->right, exponent, check_exp);
}

// ============================================================================
// 2. BUILDER HELPERS
// ============================================================================

static inline Node* make_pow(char var, double exponent) {
    return create_node(NODE_POW, 0, 0,
                       create_node(NODE_VAR, 0, var, NULL, NULL),
                       create_node(NODE_CONST, exponent, 0, NULL, NULL));
}

static inline Node* make_mul_const(double coeff, Node* term) {
    return create_node(NODE_MUL, 0, 0,
                       create_node(NODE_CONST, coeff, 0, NULL, NULL),
                       term);
}

// ============================================================================
// 3. MULTIPLICATION SIMPLIFIER
// ============================================================================

static Node* simplify_mul(Node* node) {
    Node* L = node->left;
    Node* R = node->right;

    // x^a * x^b -> x^(a + b)
    if (is_pow(L, 0, 0, false) && is_pow(R, L->left->var_name, 0, false)) {
        double new_pow = L->right->val + R->right->val;
        char var = L->left->var_name;
        free_tree(node);
        return make_pow(var, new_pow);
    }

    // x^a * x -> x^(a + 1)
    if (is_pow(L, 0, 0, false) && is_var(R, L->left->var_name)) {
        double new_pow = L->right->val + 1.0;
        char var = R->var_name;
        free_tree(node);
        return make_pow(var, new_pow);
    }

    // x * x^a -> x^(a + 1)
    if (is_var(L, 0) && is_pow(R, L->var_name, 0, false)) {
        double new_pow = R->right->val + 1.0;
        char var = L->var_name;
        free_tree(node);
        return make_pow(var, new_pow);
    }

    // x * x -> x^2
    if (is_var(L, 0) && is_var(R, L->var_name)) {
        char var = L->var_name;
        free_tree(node);
        return make_pow(var, 2.0);
    }

    return node;
}

// ============================================================================
// 4. ADDITION SIMPLIFIER
// ============================================================================

static Node* simplify_add(Node* node) {
    Node* L = node->left;
    Node* R = node->right;

    // x^n + x^n -> 2 * x^n
    if (is_pow(L, 0, 0, false) && is_pow(R, L->left->var_name, L->right->val, true)) {
        Node* term = L;
        node->left = NULL; // detach so free_tree doesn't destroy L
        free_tree(node);
        return make_mul_const(2.0, term);
    }

    // c1*x^n + x^n -> (c1 + 1)*x^n
    if (L && L->type == NODE_MUL && is_const(L->left, 0, false) &&
        is_pow(L->right, 0, 0, false) &&
        is_pow(R, L->right->left->var_name, L->right->right->val, true)) {
        
        double new_coeff = L->left->val + 1.0;
        Node* term = R;
        node->right = NULL;
        free_tree(node);
        return make_mul_const(new_coeff, term);
    }

    // Symmetric case: x^n + c1*x^n -> (c1 + 1)*x^n
    if (R && R->type == NODE_MUL && is_const(R->left, 0, false) &&
        is_pow(R->right, 0, 0, false) &&
        is_pow(L, R->right->left->var_name, R->right->right->val, true)) {
        
        double new_coeff = R->left->val + 1.0;
        Node* term = L;
        node->left = NULL;
        free_tree(node);
        return make_mul_const(new_coeff, term);
    }

    return node;
}

// ============================================================================
// 5. MAIN ENTRY POINT
// ============================================================================

Node* simplify_tree(Node* node) {
    if (!node) return NULL;

    node->left = simplify_tree(node->left);
    node->right = simplify_tree(node->right);

    if (node->type == NODE_MUL) return simplify_mul(node);
    if (node->type == NODE_ADD) return simplify_add(node);

    return node;
}

/* Helper to count nodes in an AST */
static int count_nodes(Node* node) {
    if (!node) return 0;
    return 1 + count_nodes(node->left) + count_nodes(node->right);
}

/* Repeatedly simplifies tree until no further structural reduction occurs */
Node* simplify_full(Node* node) {
    if (!node) return NULL;
    int prev_count = 0;
    int current_count = count_nodes(node);

    while (prev_count != current_count) {
        prev_count = current_count;
        node = simplify_tree(node);
        current_count = count_nodes(node);
    }
    return node;
}

void print_inorder(Node* node) {
    if (!node) return;

    if (node->type == NODE_EXP) {
        printf("exp(");
        print_inorder(node->left);
        printf(")");
        return;
    }

    int need_parens = (node->type != NODE_CONST && node->type != NODE_VAR);
    if (need_parens) printf("(");

    print_inorder(node->left);

    switch (node->type) {
        case NODE_CONST: printf("%.2f", node->val); break;
        case NODE_VAR:   printf("%c", node->var_name); break;
        case NODE_ADD:   printf(" + "); break;
        case NODE_SUB:   printf(" - "); break;
        case NODE_MUL:   printf(" * "); break;
        case NODE_DIV:   printf(" / "); break;
        case NODE_POW:   printf("^"); break;
        default: break;
    }

    print_inorder(node->right);

    if (need_parens) printf(")");
}