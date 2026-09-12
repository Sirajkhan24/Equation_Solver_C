#include <stdio.h>
#include <math.h>
#include "expr.h"
#include "parser.h"

int main(void) {
    system("cls");
    printf("=========================================================\n");
    printf("           SYMBOLIC CALCULUS VERIFICATION TEST           \n");
    printf("=========================================================\n\n");



    
    // 1. Define initial equation string f(x)
    const char* original_str = "2^x + 2^x";
    double test_x = 2.0;

    // 2. Parse string into AST f(x) and simplify initial structure
    Node* f_x = parse_expression(original_str);
    f_x = simplify_full(f_x);
    
    printf(" f(x): ");
    print_inorder(f_x);
    printf("\n");
    double val_fx = evaluate(f_x, test_x);
    printf(" f(%.1f) = %.4f\n\n", test_x, val_fx);


    // Clean up memory
    free_tree(f_x);
    return 0;
}