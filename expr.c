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
// 1. TERM COLLECTION STRUCTURES
// ============================================================================

typedef struct {
    char var;        // Variable name ('x', 'y', or 0 for constant)
    double exponent; // Exponent (1.0 for x, 0.0 for const, n for x^n)
    double coeff;    // Aggregated numeric multiplier
} Term;

typedef struct {
    Term* items;
    size_t count;
    size_t capacity;
} TermArray;

static void append_term(TermArray* arr, char var, double exp, double coeff) {
    // Combine with an existing term if variable and power match
    for (size_t i = 0; i < arr->count; i++) {
        if (arr->items[i].var == var && fabs(arr->items[i].exponent - exp) < 1e-9) {
            arr->items[i].coeff += coeff;
            return;
        }
    }
    // Expand buffer if needed
    if (arr->count >= arr->capacity) {
        arr->capacity = (arr->capacity == 0) ? 8 : arr->capacity * 2;
        arr->items = (Term*)realloc(arr->items, arr->capacity * sizeof(Term));
    }
    arr->items[arr->count++] = (Term){ .var = var, .exponent = exp, .coeff = coeff };
}

// ============================================================================
// 2. TREE FLATTENING (N-ARY ADDITION RECURSION)
// ============================================================================

static void collect_addition_terms(Node* node, TermArray* arr, double scale) {
    if (!node) return;

    // Case 1: Sub-addition -> Recurse down both sides
    if (node->type == NODE_ADD) {
        collect_addition_terms(node->left, arr, scale);
        collect_addition_terms(node->right, arr, scale);
        return;
    }

    // Case 2: Scaled term -> c * x^n or c * x
    if (node->type == NODE_MUL && node->left && node->left->type == NODE_CONST) {
        double current_coeff = node->left->val * scale;
        Node* term = node->right;

        if (term->type == NODE_VAR) {
            append_term(arr, term->var_name, 1.0, current_coeff);
            return;
        }
        if (term->type == NODE_POW && term->left && term->left->type == NODE_VAR &&
            term->right && term->right->type == NODE_CONST) {
            append_term(arr, term->left->var_name, term->right->val, current_coeff);
            return;
        }
    }

    // Case 3: Unscaled Power -> x^n
    if (node->type == NODE_POW && node->left && node->left->type == NODE_VAR &&
        node->right && node->right->type == NODE_CONST) {
        append_term(arr, node->left->var_name, node->right->val, scale);
        return;
    }

    // Case 4: Single Variable -> x
    if (node->type == NODE_VAR) {
        append_term(arr, node->var_name, 1.0, scale);
        return;
    }

    // Case 5: Bare Constant
    if (node->type == NODE_CONST) {
        append_term(arr, 0, 0.0, node->val * scale);
        return;
    }
}

// ============================================================================
// 3. CANONICAL TREE RECONSTRUCTION
// ============================================================================

static Node* build_term_node(Term t) {
    // Constant term
    if (t.var == 0 || t.exponent == 0.0) {
        return create_node(NODE_CONST, t.coeff, 0, NULL, NULL);
    }

    // Base variable or power node: x vs x^n
    Node* base = NULL;
    if (fabs(t.exponent - 1.0) < 1e-9) {
        base = create_node(NODE_VAR, 0, t.var, NULL, NULL);
    } else {
        base = create_node(NODE_POW, 0, 0,
                           create_node(NODE_VAR, 0, t.var, NULL, NULL),
                           create_node(NODE_CONST, t.exponent, 0, NULL, NULL));
    }

    // Return term directly if coefficient is 1.0
    if (fabs(t.coeff - 1.0) < 1e-9) {
        return base;
    }

    // Wrap with multiplier node: c * base
    return create_node(NODE_MUL, 0, 0,
                       create_node(NODE_CONST, t.coeff, 0, NULL, NULL),
                       base);
}

static Node* reconstruct_addition_tree(TermArray* arr) {
    Node* result = NULL;

    for (size_t i = 0; i < arr->count; i++) {
        if (fabs(arr->items[i].coeff) < 1e-9) continue; // Omit 0.0 terms

        Node* term_node = build_term_node(arr->items[i]);
        if (!result) {
            result = term_node;
        } else {
            result = create_node(NODE_ADD, 0, 0, result, term_node);
        }
    }

    return result ? result : create_node(NODE_CONST, 0.0, 0, NULL, NULL);
}

// ============================================================================
// MULTIPLICATION FLATTENING & RECONSTRUCTION
// ============================================================================

typedef struct {
    char var;
    double exponent;
} MulVariable;

typedef struct {
    MulVariable* items;
    size_t count;
    size_t capacity;
} MulVarArray;

static void add_variable_exponent(MulVarArray* arr, char var, double exp) {
    for (size_t i = 0; i < arr->count; i++) {
        if (arr->items[i].var == var) {
            arr->items[i].exponent += exp;
            return;
        }
    }
    if (arr->count >= arr->capacity) {
        arr->capacity = (arr->capacity == 0) ? 4 : arr->capacity * 2;
        arr->items = (MulVarArray*)realloc(arr->items, arr->capacity * sizeof(MulVariable));
    }
    arr->items[arr->count++] = (MulVariable){ .var = var, .exponent = exp };
}

static void collect_multiplication_terms(Node* node, MulVarArray* vars, double* total_coeff) {
    if (!node) return;

    // Case 1: Sub-multiplication -> Recurse left and right
    if (node->type == NODE_MUL) {
        collect_multiplication_terms(node->left, vars, total_coeff);
        collect_multiplication_terms(node->right, vars, total_coeff);
        return;
    }

    // Case 2: Constant factor
    if (node->type == NODE_CONST) {
        *total_coeff *= node->val;
        return;
    }

    // Case 3: Power term (x^n)
    if (node->type == NODE_POW && node->left && node->left->type == NODE_VAR &&
        node->right && node->right->type == NODE_CONST) {
        add_variable_exponent(vars, node->left->var_name, node->right->val);
        return;
    }

    // Case 4: Plain Variable (x)
    if (node->type == NODE_VAR) {
        add_variable_exponent(vars, node->var_name, 1.0);
        return;
    }
}

static Node* reconstruct_multiplication_tree(MulVarArray* vars, double total_coeff) {
    Node* result = NULL;

    // 1. Build variable terms (x^n or x)
    for (size_t i = 0; i < vars->count; i++) {
        if (fabs(vars->items[i].exponent) < 1e-9) continue; // x^0 = 1

        Node* var_node = NULL;
        if (fabs(vars->items[i].exponent - 1.0) < 1e-9) {
            var_node = create_node(NODE_VAR, 0, vars->items[i].var, NULL, NULL);
        } else {
            var_node = create_node(NODE_POW, 0, 0,
                                   create_node(NODE_VAR, 0, vars->items[i].var, NULL, NULL),
                                   create_node(NODE_CONST, vars->items[i].exponent, 0, NULL, NULL));
        }

        if (!result) {
            result = var_node;
        } else {
            result = create_node(NODE_MUL, 0, 0, result, var_node);
        }
    }

    // 2. Attach total numeric coefficient if not 1.0 (or if expression is purely constant)
    if (!result) {
        return create_node(NODE_CONST, total_coeff, 0, NULL, NULL);
    }

    if (fabs(total_coeff - 1.0) > 1e-9) {
        result = create_node(NODE_MUL, 0, 0,
                             create_node(NODE_CONST, total_coeff, 0, NULL, NULL),
                             result);
    }

    return result;
}

// ============================================================================
// 4. MAIN SIMPLIFIER ENTRY POINT
// ============================================================================

Node* simplify_tree(Node* node) {
    if (!node) return NULL;

    // 1. Post-order bottom-up recursion
    node->left = simplify_tree(node->left);
    node->right = simplify_tree(node->right);

    // 2. N-ary Multiplication Collector
    if (node->type == NODE_MUL) {
        MulVarArray vars = {0};
        double total_coeff = 1.0;

        collect_multiplication_terms(node, &vars, &total_coeff);
        Node* simplified = reconstruct_multiplication_tree(&vars, total_coeff);

        free(vars.items);
        free_tree(node); // Safely clean up uncollected binary subtrees
        return simplified;
    }

    // 3. N-ary Addition Collector
    if (node->type == NODE_ADD) {
        TermArray arr = {0};
        collect_addition_terms(node, &arr, 1.0);

        Node* simplified = reconstruct_addition_tree(&arr);

        free(arr.items);
        free_tree(node);
        return simplified;
    }

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