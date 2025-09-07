#include "harness.h"

#define CHECK_PARSE(var, expr, allowVars) \
    do { \
      var = Expr::From(expr, allowVars, false); \
      CHECK_TRUE(var != NULL); \
    } while(0)

#define CHECK_PARSE_ERR(expr, msg) \
    do { \
      std::string err; \
      Expr *e = Expr::Parse(expr, false, nullptr, &err); \
      CHECK_TRUE(e == NULL); \
      CHECK_TRUE(err.find(msg) != std::string::npos); \
    } while(0)

TEST_CASE(constant) {
  Expr *e;
  CHECK_PARSE(e, "pi", false);
  CHECK_EQ_EPS(e->Eval({}), 3.1415926);
}

TEST_CASE(literal) {
  Expr *e;
  CHECK_PARSE(e, "42", false);
  CHECK_TRUE(e->Eval({}) == 42);
  CHECK_PARSE(e, "42.5", false);
  CHECK_TRUE(e->Eval({}) == 42.5);
  CHECK_PARSE(e, "1_000_000", false);
  CHECK_TRUE(e->Eval({}) == 1000000);
}

TEST_CASE(unary_ops) {
  Expr *e;
  CHECK_PARSE(e, "-10", false);
  CHECK_TRUE(e->Eval({}) == -10);
}

TEST_CASE(binary_ops) {
  Expr *e;
  CHECK_PARSE(e, "1 + 2", false);
  CHECK_TRUE(e->Eval({}) == 3);
  CHECK_PARSE(e, "1 - 2", false);
  CHECK_TRUE(e->Eval({}) == -1);
  CHECK_PARSE(e, "3 * 4", false);
  CHECK_TRUE(e->Eval({}) == 12);
  CHECK_PARSE(e, "3 / 4", false);
  CHECK_TRUE(e->Eval({}) == 0.75);
}

TEST_CASE(parentheses) {
  Expr *e;
  CHECK_PARSE(e, "(1 + 2) * 3", false);
  CHECK_TRUE(e->Eval({}) == 9);
  CHECK_PARSE(e, "1 + (2 * 3)", false);
  CHECK_TRUE(e->Eval({}) == 7);
}

TEST_CASE(functions) {
  Expr *e;
  CHECK_PARSE(e, "sqrt(2)", false);
  CHECK_EQ_EPS(e->Eval({}), 1.414213);
  CHECK_PARSE(e, "square(3)", false);
  CHECK_EQ_EPS(e->Eval({}), 9);
  CHECK_PARSE(e, "sin(180)", false);
  CHECK_EQ_EPS(e->Eval({}), 0);
  CHECK_PARSE(e, "sin(90)", false);
  CHECK_EQ_EPS(e->Eval({}), 1);
  CHECK_PARSE(e, "cos(180)", false);
  CHECK_EQ_EPS(e->Eval({}), -1);
  CHECK_PARSE(e, "asin(1)", false);
  CHECK_EQ_EPS(e->Eval({}), 90);
  CHECK_PARSE(e, "acos(0)", false);
  CHECK_EQ_EPS(e->Eval({}), 90);
}

TEST_CASE(variable) {
  Expr *e;
  CHECK_PARSE(e, "var", true);
  CHECK_TRUE(e->op == Expr::Op::VARIABLE);
  CHECK_TRUE(std::string(e->s) == "var");
  CHECK_PARSE(e, "Var", true);
  CHECK_TRUE(e->op == Expr::Op::VARIABLE);
  CHECK_TRUE(std::string(e->s) == "Var");
  CHECK_PARSE(e, "a123", true);
  CHECK_TRUE(e->op == Expr::Op::VARIABLE);
  CHECK_TRUE(std::string(e->s) == "a123");
  CHECK_PARSE(e, "a_123_", true);
  CHECK_TRUE(e->op == Expr::Op::VARIABLE);
  CHECK_TRUE(std::string(e->s) == "a_123_");
}

TEST_CASE(precedence) {
  Expr *e;
  CHECK_PARSE(e, "2 + 3 * 4", false);
  CHECK_TRUE(e->Eval({}) == 14);
  CHECK_PARSE(e, "2 - 3 / 4", false);
  CHECK_TRUE(e->Eval({}) == 1.25);
  CHECK_PARSE(e, "-3 + 2", false);
  CHECK_TRUE(e->Eval({}) == -1);
  CHECK_PARSE(e, "2 + 3 - 4", false);
  CHECK_TRUE(e->Eval({}) == 1);
}

TEST_CASE(errors) {
  CHECK_PARSE_ERR("\x01",
                  "Unexpected character");
  CHECK_PARSE_ERR("notavar",
                  "Referencing variable 'notavar' is not allowed in this expression");
  CHECK_PARSE_ERR("1e2e3",
                  "'1e2e3' is not a valid number");
  CHECK_PARSE_ERR("_",
                  "'_' is not a valid operator");
  CHECK_PARSE_ERR("2 2",
                  "Expected an operator");
  CHECK_PARSE_ERR("2 + +",
                  "Expected an operand");
  CHECK_PARSE_ERR("( 2 + 2",
                  "Expected ')'");
  CHECK_PARSE_ERR("(",
                  "Expected ')'");
}
