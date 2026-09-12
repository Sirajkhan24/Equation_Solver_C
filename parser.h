// Parser interface (string -> AST)

#ifndef PARSER_H
#define PARSER_H

#include "expr.h"

#define STACK_MAX 100

/* Parse infix math expression string into an AST */
Node* parse_expression(const char* expr);

#endif /* PARSER_H */