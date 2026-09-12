// Differentiation, integration, evaluation, printing
#include <stdbool.h>
#include <stdlib.h>
#include "expr.h"

Node *create_node(NodeType type, double val, char var_name, Node *left, Node *right)
{
    Node *n = (Node *)malloc(sizeof(Node));
    if (!n)
    {
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

void free_tree(Node *node)
{
    if (!node)
        return;
    free_tree(node->left);
    free_tree(node->right);
    free(node);
}

Node *copy_tree(Node *node)
{
    if (!node)
        return NULL;
    return create_node(node->type, node->val, node->var_name,
                       copy_tree(node->left), copy_tree(node->right));
}

double evaluate(Node *node, double x_val)
{
    if (!node)
        return 0.0;
    switch (node->type)
    {
    case NODE_CONST:
        return node->val;
    case NODE_VAR:
        return x_val;
    case NODE_ADD:
        return evaluate(node->left, x_val) + evaluate(node->right, x_val);
    case NODE_SUB:
        return evaluate(node->left, x_val) - evaluate(node->right, x_val);
    case NODE_MUL:
        return evaluate(node->left, x_val) * evaluate(node->right, x_val);
    case NODE_DIV:
        return evaluate(node->left, x_val) / evaluate(node->right, x_val);
    case NODE_POW:
        return pow(evaluate(node->left, x_val), evaluate(node->right, x_val));
    case NODE_EXP:
        return exp(evaluate(node->left, x_val));
    case NODE_SIN:
        return sin(evaluate(node->left, x_val));
    case NODE_COS:
        return cos(evaluate(node->left, x_val));
    default:
        return 0.0;
    }
}

/* Definite Numerical Integration using Simpson's 1/3 Rule */
double integrate_numerical(Node *root, double a, double b, int n)
{
    if (n % 2 != 0)
        n++; // Simpson's rule requires an even number of intervals

    double h = (b - a) / n;
    double sum = evaluate(root, a) + evaluate(root, b);

    for (int i = 1; i < n; i++)
    {
        double x = a + i * h;
        if (i % 2 == 0)
        {
            sum += 2.0 * evaluate(root, x);
        }
        else
        {
            sum += 4.0 * evaluate(root, x);
        }
    }

    return (h / 3.0) * sum;
}

Node *differentiate(Node *node, char var)
{
    if (!node)
        return NULL;

    switch (node->type)
    {
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
            Node *num = create_node(NODE_SUB, 0, 0,
                                    create_node(NODE_MUL, 0, 0, differentiate(node->left, var), copy_tree(node->right)),
                                    create_node(NODE_MUL, 0, 0, copy_tree(node->left), differentiate(node->right, var)));
            Node *den = create_node(NODE_POW, 0, 0, copy_tree(node->right), create_node(NODE_CONST, 2.0, 0, NULL, NULL));
            return create_node(NODE_DIV, 0, 0, num, den);
        }

    case NODE_EXP:
        /* Chain Rule: (exp(u))' = exp(u) * u' */
        return create_node(NODE_MUL, 0, 0,
                           create_node(NODE_EXP, 0, 0, copy_tree(node->left), NULL),
                           differentiate(node->left, var));

    case NODE_POW:
        /* Power Rule for polynomial terms (u^c): (u^n)' = n * u^(n-1) * u' */
        if (node->right && node->right->type == NODE_CONST)
        {
            double n = node->right->val;
            Node *pow_term = create_node(NODE_POW, 0, 0,
                                         copy_tree(node->left),
                                         create_node(NODE_CONST, n - 1.0, 0, NULL, NULL));
            Node *coeff_term = create_node(NODE_MUL, 0, 0,
                                           create_node(NODE_CONST, n, 0, NULL, NULL),
                                           pow_term);
            return create_node(NODE_MUL, 0, 0, coeff_term, differentiate(node->left, var));
        }

        /* Exponential Rule for constant base (a^u): (a^u)' = ln(a) * a^u * u' */
        if (node->left && node->left->type == NODE_CONST)
        {
            double a = node->left->val;
            Node *ln_base = create_node(NODE_CONST, log(a), 0, NULL, NULL);
            Node *pow_term = create_node(NODE_MUL, 0, 0, ln_base, copy_tree(node));
            return create_node(NODE_MUL, 0, 0, pow_term, differentiate(node->right, var));
        }

        return create_node(NODE_CONST, 0.0, 0, NULL, NULL);

    case NODE_SIN:
        /* d/dx(sin(u)) = cos(u) * u' */
        return create_node(NODE_MUL, 0, 0,
                           create_node(NODE_COS, 0, 0, copy_tree(node->left), NULL),
                           differentiate(node->left, var));

    case NODE_COS:
        /* d/dx(cos(u)) = -1 * sin(u) * u' */
        return create_node(NODE_MUL, 0, 0,
                           create_node(NODE_MUL, 0, 0,
                                       create_node(NODE_CONST, -1.0, 0, NULL, NULL),
                                       create_node(NODE_SIN, 0, 0, copy_tree(node->left), NULL)),
                           differentiate(node->left, var));

    default:
        return NULL;
    }
}

Node *integrate(Node *node, char var)
{
    if (!node)
        return NULL;

    switch (node->type)
    {
    case NODE_CONST:
        // Rule: ∫ c dx = c * x
        if (node->val == 0.0)
        {
            return create_node(NODE_CONST, 0.0, 0, NULL, NULL);
        }
        return create_node(NODE_MUL, 0, 0,
                           copy_tree(node),
                           create_node(NODE_VAR, 0, var, NULL, NULL));

    case NODE_VAR:
        // Rule: ∫ x dx = (1/2) * x^2
        if (node->var_name == var)
        {
            Node *pow_term = create_node(NODE_POW, 0, 0,
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
        if (node->left && node->left->type == NODE_CONST)
        {
            return create_node(NODE_MUL, 0, 0,
                               copy_tree(node->left),
                               integrate(node->right, var));
        }
        if (node->right && node->right->type == NODE_CONST)
        {
            return create_node(NODE_MUL, 0, 0,
                               integrate(node->left, var),
                               copy_tree(node->right));
        }
        fprintf(stderr, "Warning: Symbolic product integration (Integration by Parts) not fully supported.\n");
        return NULL;

    case NODE_POW:
        // Power Rule: ∫ (x^n) dx = (x^(n+1)) / (n+1)   [for n != -1]
        if (node->left && node->left->type == NODE_VAR && node->left->var_name == var &&
            node->right && node->right->type == NODE_CONST)
        {

            double n = node->right->val;
            if (n == -1.0)
            {
                fprintf(stderr, "Error: ∫ x^-1 dx requires ln(x) support.\n");
                return NULL;
            }

            double new_power = n + 1.0;
            Node *pow_term = create_node(NODE_POW, 0, 0,
                                         copy_tree(node->left),
                                         create_node(NODE_CONST, new_power, 0, NULL, NULL));

            return create_node(NODE_DIV, 0, 0,
                               pow_term,
                               create_node(NODE_CONST, new_power, 0, NULL, NULL));
        }
        return NULL;

    case NODE_EXP:
        // Rule: ∫ exp(x) dx = exp(x)
        if (node->left && node->left->type == NODE_VAR && node->left->var_name == var)
        {
            return copy_tree(node);
        }
        // Rule: ∫ exp(c * x) dx = (1/c) * exp(c * x)
        if (node->left && node->left->type == NODE_MUL)
        {
            Node *mul = node->left;
            if (mul->left && mul->left->type == NODE_CONST &&
                mul->right && mul->right->type == NODE_VAR && mul->right->var_name == var)
            {

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
// UNIFIED TERM COLLECTION & SIMPLIFICATION (PRESERVES ALL SUBTREES)
// ============================================================================

typedef enum
{
    TERM_CONST,    /* Bare constant c */
    TERM_POLY,     /* Polynomial term c * x^n */
    TERM_EXP_BASE, /* Exponential term c * a^(k*x) */
    TERM_EXP_NAT,  /* Natural exponential c * exp(x) */
    TERM_CUSTOM    /* Unrecognized subtree (fallback - never dropped) */
} TermType;

typedef struct
{
    TermType type;
    double coeff;      /* Scalar multiplier c */
    char var;          /* Variable identifier 'x' */
    double exponent;   /* Power n in x^n */
    double base;       /* Constant base a in a^(k*x) */
    double var_coeff;  /* Exponent scale k in a^(k*x) */
    Node *custom_node; /* AST subtree for fallback */
} Term;

typedef struct
{
    Term *items;
    size_t count;
    size_t capacity;
} TermArray;

static void append_term_item(TermArray *arr, Term t)
{
    /* Aggregate coefficients for matching terms */
    for (size_t i = 0; i < arr->count; i++)
    {
        Term *existing = &arr->items[i];
        if (existing->type == t.type)
        {
            if (t.type == TERM_CONST)
            {
                existing->coeff += t.coeff;
                return;
            }
            if (t.type == TERM_POLY && existing->var == t.var && fabs(existing->exponent - t.exponent) < 1e-9)
            {
                existing->coeff += t.coeff;
                return;
            }
            if (t.type == TERM_EXP_BASE && existing->var == t.var &&
                fabs(existing->base - t.base) < 1e-9 && fabs(existing->var_coeff - t.var_coeff) < 1e-9)
            {
                existing->coeff += t.coeff;
                return;
            }
            if (t.type == TERM_EXP_NAT && existing->var == t.var)
            {
                existing->coeff += t.coeff;
                return;
            }
        }
    }
    if (arr->count >= arr->capacity)
    {
        arr->capacity = (arr->capacity == 0) ? 8 : arr->capacity * 2;
        arr->items = (Term *)realloc(arr->items, arr->capacity * sizeof(Term));
    }
    arr->items[arr->count++] = t;
}

static void collect_addition_terms(Node *node, TermArray *arr, double scale)
{
    if (!node)
        return;

    if (node->type == NODE_ADD)
    {
        collect_addition_terms(node->left, arr, scale);
        collect_addition_terms(node->right, arr, scale);
        return;
    }

    if (node->type == NODE_SUB)
    {
        collect_addition_terms(node->left, arr, scale);
        collect_addition_terms(node->right, arr, -scale);
        return;
    }

    double c = scale;
    Node *term = node;
    if (node->type == NODE_MUL && node->left && node->left->type == NODE_CONST)
    {
        c = node->left->val * scale;
        term = node->right;
    }

    if (term->type == NODE_CONST)
    {
        append_term_item(arr, (Term){.type = TERM_CONST, .coeff = c * term->val});
        return;
    }

    if (term->type == NODE_VAR)
    {
        append_term_item(arr, (Term){.type = TERM_POLY, .coeff = c, .var = term->var_name, .exponent = 1.0});
        return;
    }

    if (term->type == NODE_POW)
    {
        /* Polynomial x^n */
        if (term->left && term->left->type == NODE_VAR && term->right && term->right->type == NODE_CONST)
        {
            append_term_item(arr, (Term){.type = TERM_POLY, .coeff = c, .var = term->left->var_name, .exponent = term->right->val});
            return;
        }
        /* Constant base a^x */
        if (term->left && term->left->type == NODE_CONST && term->right && term->right->type == NODE_VAR)
        {
            append_term_item(arr, (Term){.type = TERM_EXP_BASE, .coeff = c, .base = term->left->val, .var = term->right->var_name, .var_coeff = 1.0});
            return;
        }
        /* Constant base with scaled exponent a^(k*x) */
        if (term->left && term->left->type == NODE_CONST && term->right && term->right->type == NODE_MUL &&
            term->right->left && term->right->left->type == NODE_CONST &&
            term->right->right && term->right->right->type == NODE_VAR)
        {
            append_term_item(arr, (Term){.type = TERM_EXP_BASE, .coeff = c, .base = term->left->val, .var = term->right->right->var_name, .var_coeff = term->right->left->val});
            return;
        }
    }

    if (term->type == NODE_EXP && term->left && term->left->type == NODE_VAR)
    {
        append_term_item(arr, (Term){.type = TERM_EXP_NAT, .coeff = c, .var = term->left->var_name});
        return;
    }

    /* Preserve unclassified subtrees without dropping them */
    append_term_item(arr, (Term){.type = TERM_CUSTOM, .coeff = c, .custom_node = copy_tree(term)});
}

static Node *build_term_node(Term t)
{
    Node *base_node = NULL;

    switch (t.type)
    {
    case TERM_CONST:
        return create_node(NODE_CONST, t.coeff, 0, NULL, NULL);

    case TERM_POLY:
        if (fabs(t.exponent - 1.0) < 1e-9)
        {
            base_node = create_node(NODE_VAR, 0, t.var, NULL, NULL);
        }
        else
        {
            base_node = create_node(NODE_POW, 0, 0,
                                    create_node(NODE_VAR, 0, t.var, NULL, NULL),
                                    create_node(NODE_CONST, t.exponent, 0, NULL, NULL));
        }
        break;

    case TERM_EXP_BASE:
    {
        Node *exp_node = NULL;
        if (fabs(t.var_coeff - 1.0) < 1e-9)
        {
            exp_node = create_node(NODE_VAR, 0, t.var, NULL, NULL);
        }
        else
        {
            exp_node = create_node(NODE_MUL, 0, 0,
                                   create_node(NODE_CONST, t.var_coeff, 0, NULL, NULL),
                                   create_node(NODE_VAR, 0, t.var, NULL, NULL));
        }
        base_node = create_node(NODE_POW, 0, 0,
                                create_node(NODE_CONST, t.base, 0, NULL, NULL),
                                exp_node);
        break;
    }

    case TERM_EXP_NAT:
        base_node = create_node(NODE_EXP, 0, 0,
                                create_node(NODE_VAR, 0, t.var, NULL, NULL),
                                NULL);
        break;

    case TERM_CUSTOM:
        base_node = copy_tree(t.custom_node);
        break;
    }

    if (fabs(t.coeff - 1.0) < 1e-9)
        return base_node;

    return create_node(NODE_MUL, 0, 0,
                       create_node(NODE_CONST, t.coeff, 0, NULL, NULL),
                       base_node);
}

static Node *reconstruct_addition_tree(TermArray *arr)
{
    Node *result = NULL;

    for (size_t i = 0; i < arr->count; i++)
    {
        if (fabs(arr->items[i].coeff) < 1e-9)
            continue;

        Node *term_node = build_term_node(arr->items[i]);
        if (!result)
        {
            result = term_node;
        }
        else
        {
            result = create_node(NODE_ADD, 0, 0, result, term_node);
        }
    }

    return result ? result : create_node(NODE_CONST, 0.0, 0, NULL, NULL);
}

// ============================================================================
// MULTIPLICATION FLATTENING & RECONSTRUCTION
// ============================================================================

typedef enum
{
    MUL_POLY,     /* x^n */
    MUL_EXP_BASE, /* a^(k*x) */
    MUL_CUSTOM    /* Subtree factor */
} MulType;

typedef struct
{
    MulType type;
    char var;
    double exponent;
    double base;
    double var_coeff;
    Node *custom_node;
} MulItem;

typedef struct
{
    MulItem *items;
    size_t count;
    size_t capacity;
} MulVarArray;

static void add_poly_exponent(MulVarArray *arr, char var, double exp)
{
    for (size_t i = 0; i < arr->count; i++)
    {
        if (arr->items[i].type == MUL_POLY && arr->items[i].var == var)
        {
            arr->items[i].exponent += exp;
            return;
        }
    }
    if (arr->count >= arr->capacity)
    {
        arr->capacity = (arr->capacity == 0) ? 4 : arr->capacity * 2;
        arr->items = (MulItem *)realloc(arr->items, arr->capacity * sizeof(MulItem));
    }
    arr->items[arr->count++] = (MulItem){.type = MUL_POLY, .var = var, .exponent = exp};
}

static void add_base_exponent(MulVarArray *arr, double base, char var, double exp_coeff)
{
    for (size_t i = 0; i < arr->count; i++)
    {
        if (arr->items[i].type == MUL_EXP_BASE && fabs(arr->items[i].base - base) < 1e-9 && arr->items[i].var == var)
        {
            arr->items[i].var_coeff += exp_coeff;
            return;
        }
    }
    if (arr->count >= arr->capacity)
    {
        arr->capacity = (arr->capacity == 0) ? 4 : arr->capacity * 2;
        arr->items = (MulItem *)realloc(arr->items, arr->capacity * sizeof(MulItem));
    }
    arr->items[arr->count++] = (MulItem){.type = MUL_EXP_BASE, .base = base, .var = var, .var_coeff = exp_coeff};
}

static void add_custom_factor(MulVarArray *arr, Node *node)
{
    if (arr->count >= arr->capacity)
    {
        arr->capacity = (arr->capacity == 0) ? 4 : arr->capacity * 2;
        arr->items = (MulItem *)realloc(arr->items, arr->capacity * sizeof(MulItem));
    }
    arr->items[arr->count++] = (MulItem){.type = MUL_CUSTOM, .custom_node = copy_tree(node)};
}

static void collect_multiplication_terms(Node *node, MulVarArray *vars, double *total_coeff)
{
    if (!node)
        return;

    if (node->type == NODE_MUL)
    {
        collect_multiplication_terms(node->left, vars, total_coeff);
        collect_multiplication_terms(node->right, vars, total_coeff);
        return;
    }

    if (node->type == NODE_CONST)
    {
        *total_coeff *= node->val;
        return;
    }

    if (node->type == NODE_VAR)
    {
        add_poly_exponent(vars, node->var_name, 1.0);
        return;
    }

    if (node->type == NODE_POW && node->left && node->left->type == NODE_VAR &&
        node->right && node->right->type == NODE_CONST)
    {
        add_poly_exponent(vars, node->left->var_name, node->right->val);
        return;
    }

    if (node->type == NODE_POW && node->left && node->left->type == NODE_CONST &&
        node->right && node->right->type == NODE_VAR)
    {
        add_base_exponent(vars, node->left->val, node->right->var_name, 1.0);
        return;
    }

    if (node->type == NODE_POW && node->left && node->left->type == NODE_CONST &&
        node->right && node->right->type == NODE_MUL &&
        node->right->left && node->right->left->type == NODE_CONST &&
        node->right->right && node->right->right->type == NODE_VAR)
    {
        add_base_exponent(vars, node->left->val, node->right->right->var_name, node->right->left->val);
        return;
    }

    add_custom_factor(vars, node);
}

static Node *reconstruct_multiplication_tree(MulVarArray *vars, double total_coeff)
{
    Node *result = NULL;

    for (size_t i = 0; i < vars->count; i++)
    {
        MulItem *item = &vars->items[i];
        Node *factor_node = NULL;

        if (item->type == MUL_CUSTOM)
        {
            factor_node = copy_tree(item->custom_node);
        }
        else if (item->type == MUL_POLY)
        {
            if (fabs(item->exponent) < 1e-9)
                continue;

            if (fabs(item->exponent - 1.0) < 1e-9)
            {
                factor_node = create_node(NODE_VAR, 0, item->var, NULL, NULL);
            }
            else
            {
                factor_node = create_node(NODE_POW, 0, 0,
                                          create_node(NODE_VAR, 0, item->var, NULL, NULL),
                                          create_node(NODE_CONST, item->exponent, 0, NULL, NULL));
            }
        }
        else if (item->type == MUL_EXP_BASE)
        {
            if (fabs(item->var_coeff) < 1e-9)
                continue;

            Node *exponent_node = NULL;
            if (fabs(item->var_coeff - 1.0) < 1e-9)
            {
                exponent_node = create_node(NODE_VAR, 0, item->var, NULL, NULL);
            }
            else
            {
                exponent_node = create_node(NODE_MUL, 0, 0,
                                            create_node(NODE_CONST, item->var_coeff, 0, NULL, NULL),
                                            create_node(NODE_VAR, 0, item->var, NULL, NULL));
            }

            factor_node = create_node(NODE_POW, 0, 0,
                                      create_node(NODE_CONST, item->base, 0, NULL, NULL),
                                      exponent_node);
        }

        if (!result)
        {
            result = factor_node;
        }
        else
        {
            result = create_node(NODE_MUL, 0, 0, result, factor_node);
        }
    }

    if (!result)
        return create_node(NODE_CONST, total_coeff, 0, NULL, NULL);

    if (fabs(total_coeff - 1.0) > 1e-9)
    {
        result = create_node(NODE_MUL, 0, 0,
                             create_node(NODE_CONST, total_coeff, 0, NULL, NULL),
                             result);
    }

    return result;
}

// ============================================================================
// MAIN SIMPLIFIER ENTRY POINT
// ============================================================================

Node *simplify_tree(Node *node)
{
    if (!node)
        return NULL;

    node->left = simplify_tree(node->left);
    node->right = simplify_tree(node->right);

    if (node->type == NODE_ADD || node->type == NODE_SUB)
    {
        TermArray arr = {0};
        collect_addition_terms(node, &arr, 1.0);
        Node *simplified = reconstruct_addition_tree(&arr);

        for (size_t i = 0; i < arr.count; i++)
        {
            if (arr.items[i].type == TERM_CUSTOM && arr.items[i].custom_node)
            {
                free_tree(arr.items[i].custom_node);
            }
        }
        free(arr.items);
        free_tree(node);
        return simplified;
    }

    if (node->type == NODE_MUL)
    {
        MulVarArray vars = {0};
        double total_coeff = 1.0;
        collect_multiplication_terms(node, &vars, &total_coeff);

        if (vars.count > 0 || fabs(total_coeff - 1.0) > 1e-9)
        {
            Node *simplified = reconstruct_multiplication_tree(&vars, total_coeff);
            for (size_t i = 0; i < vars.count; i++)
            {
                if (vars.items[i].type == MUL_CUSTOM && vars.items[i].custom_node)
                {
                    free_tree(vars.items[i].custom_node);
                }
            }
            free(vars.items);
            free_tree(node);
            return simplified;
        }
        for (size_t i = 0; i < vars.count; i++)
        {
            if (vars.items[i].type == MUL_CUSTOM && vars.items[i].custom_node)
            {
                free_tree(vars.items[i].custom_node);
            }
        }
        free(vars.items);
    }

    return node;
}

/* Helper to count nodes in an AST */
static int count_nodes(Node *node)
{
    if (!node)
        return 0;
    return 1 + count_nodes(node->left) + count_nodes(node->right);
}

/* Repeatedly simplifies tree until no further structural reduction occurs */
Node *simplify_full(Node *node)
{
    if (!node)
        return NULL;
    int prev_count = 0;
    int current_count = count_nodes(node);

    while (prev_count != current_count)
    {
        prev_count = current_count;
        node = simplify_tree(node);
        current_count = count_nodes(node);
    }
    return node;
}

void print_inorder(Node *node)
{
    if (!node)
        return;

    if (node->type == NODE_EXP)
    {
        printf("exp(");
        print_inorder(node->left);
        printf(")");
        return;
    }

    int need_parens = (node->type != NODE_CONST && node->type != NODE_VAR);
    if (need_parens)
        printf("(");

    print_inorder(node->left);

    switch (node->type)
    {
    case NODE_CONST:
        printf("%.2f", node->val);
        break;
    case NODE_VAR:
        printf("%c", node->var_name);
        break;
    case NODE_ADD:
        printf(" + ");
        break;
    case NODE_SUB:
        printf(" - ");
        break;
    case NODE_MUL:
        printf(" * ");
        break;
    case NODE_DIV:
        printf(" / ");
        break;
    case NODE_POW:
        printf("^");
        break;
    case NODE_SIN:
        printf("sin(");
        print_inorder(node->left);
        printf(")");
        break;
    case NODE_COS:
        printf("cos(");
        print_inorder(node->left);
        printf(")");
        break;
    default:
        break;
    }

    print_inorder(node->right);

    if (need_parens)
        printf(")");
}