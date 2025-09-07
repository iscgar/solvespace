//-----------------------------------------------------------------------------
// The symbolic algebra system used to write our constraint equations;
// routines to build expressions in software or from a user-provided string,
// and to compute the partial derivatives that we'll use when write our
// Jacobian matrix.
//
// Copyright 2008-2013 Jonathan Westhues.
//-----------------------------------------------------------------------------
#include "solvespace.h"

namespace SolveSpace {

static inline Expr *AllocExpr() {
    return (Expr *)Platform::AllocTemporary(sizeof(Expr));
}

ExprVector ExprVector::From(Expr *x, Expr *y, Expr *z) {
    ExprVector r = { x, y, z};
    return r;
}

ExprVector ExprVector::From(Vector vn) {
    ExprVector ve;
    ve.x = Expr::From(vn.x);
    ve.y = Expr::From(vn.y);
    ve.z = Expr::From(vn.z);
    return ve;
}

ExprVector ExprVector::From(hParam x, hParam y, hParam z) {
    ExprVector ve;
    ve.x = Expr::From(x);
    ve.y = Expr::From(y);
    ve.z = Expr::From(z);
    return ve;
}

ExprVector ExprVector::From(double x, double y, double z) {
    ExprVector ve;
    ve.x = Expr::From(x);
    ve.y = Expr::From(y);
    ve.z = Expr::From(z);
    return ve;
}

ExprVector ExprVector::Minus(ExprVector b) const {
    ExprVector r;
    r.x = x->Minus(b.x);
    r.y = y->Minus(b.y);
    r.z = z->Minus(b.z);
    return r;
}

ExprVector ExprVector::Plus(ExprVector b) const {
    ExprVector r;
    r.x = x->Plus(b.x);
    r.y = y->Plus(b.y);
    r.z = z->Plus(b.z);
    return r;
}

Expr *ExprVector::Dot(ExprVector b) const {
    Expr *r;
    r =         x->Times(b.x);
    r = r->Plus(y->Times(b.y));
    r = r->Plus(z->Times(b.z));
    return r;
}

ExprVector ExprVector::Cross(ExprVector b) const {
    ExprVector r;
    r.x = (y->Times(b.z))->Minus(z->Times(b.y));
    r.y = (z->Times(b.x))->Minus(x->Times(b.z));
    r.z = (x->Times(b.y))->Minus(y->Times(b.x));
    return r;
}

ExprVector ExprVector::ScaledBy(Expr *s) const {
    ExprVector r;
    r.x = x->Times(s);
    r.y = y->Times(s);
    r.z = z->Times(s);
    return r;
}

ExprVector ExprVector::WithMagnitude(Expr *s) const {
    Expr *m = Magnitude();
    return ScaledBy(s->Div(m));
}

Expr *ExprVector::Magnitude() const {
    Expr *r;
    r =         x->Square();
    r = r->Plus(y->Square());
    r = r->Plus(z->Square());
    return r->Sqrt();
}

Vector ExprVector::Eval(const ResolutionMap &resolutions) const {
    Vector r;
    r.x = x->Eval(resolutions);
    r.y = y->Eval(resolutions);
    r.z = z->Eval(resolutions);
    return r;
}

bool ExprVector::Resolve(const ResolutionMap &resolutions) {
    return x->Resolve(resolutions) && y->Resolve(resolutions) && z->Resolve(resolutions);
}

ExprQuaternion ExprQuaternion::From(hParam w, hParam vx, hParam vy, hParam vz) {
    ExprQuaternion q;
    q.w  = Expr::From(w);
    q.vx = Expr::From(vx);
    q.vy = Expr::From(vy);
    q.vz = Expr::From(vz);
    return q;
}

ExprQuaternion ExprQuaternion::From(Expr *w, Expr *vx, Expr *vy, Expr *vz)
{
    ExprQuaternion q;
    q.w = w;
    q.vx = vx;
    q.vy = vy;
    q.vz = vz;
    return q;
}

ExprQuaternion ExprQuaternion::From(Quaternion qn) {
    ExprQuaternion qe;
    qe.w = Expr::From(qn.w);
    qe.vx = Expr::From(qn.vx);
    qe.vy = Expr::From(qn.vy);
    qe.vz = Expr::From(qn.vz);
    return qe;
}

ExprVector ExprQuaternion::RotationU() const {
    ExprVector u;
    Expr *two = Expr::From(2);

    u.x = w->Square();
    u.x = (u.x)->Plus(vx->Square());
    u.x = (u.x)->Minus(vy->Square());
    u.x = (u.x)->Minus(vz->Square());

    u.y = two->Times(w->Times(vz));
    u.y = (u.y)->Plus(two->Times(vx->Times(vy)));

    u.z = two->Times(vx->Times(vz));
    u.z = (u.z)->Minus(two->Times(w->Times(vy)));

    return u;
}

ExprVector ExprQuaternion::RotationV() const {
    ExprVector v;
    Expr *two = Expr::From(2);

    v.x = two->Times(vx->Times(vy));
    v.x = (v.x)->Minus(two->Times(w->Times(vz)));

    v.y = w->Square();
    v.y = (v.y)->Minus(vx->Square());
    v.y = (v.y)->Plus(vy->Square());
    v.y = (v.y)->Minus(vz->Square());

    v.z = two->Times(w->Times(vx));
    v.z = (v.z)->Plus(two->Times(vy->Times(vz)));

    return v;
}

ExprVector ExprQuaternion::RotationN() const {
    ExprVector n;
    Expr *two = Expr::From(2);

    n.x =              two->Times( w->Times(vy));
    n.x = (n.x)->Plus (two->Times(vx->Times(vz)));

    n.y =              two->Times(vy->Times(vz));
    n.y = (n.y)->Minus(two->Times( w->Times(vx)));

    n.z =               w->Square();
    n.z = (n.z)->Minus(vx->Square());
    n.z = (n.z)->Minus(vy->Square());
    n.z = (n.z)->Plus (vz->Square());

    return n;
}

ExprVector ExprQuaternion::Rotate(ExprVector p) const {
    // Express the point in the new basis
    return (RotationU().ScaledBy(p.x)).Plus(
            RotationV().ScaledBy(p.y)).Plus(
            RotationN().ScaledBy(p.z));
}

ExprQuaternion ExprQuaternion::Times(ExprQuaternion b) const {
    Expr *sa = w, *sb = b.w;
    ExprVector va = { vx, vy, vz };
    ExprVector vb = { b.vx, b.vy, b.vz };

    ExprQuaternion r;
    r.w = (sa->Times(sb))->Minus(va.Dot(vb));
    ExprVector vr = vb.ScaledBy(sa).Plus(
                    va.ScaledBy(sb).Plus(
                    va.Cross(vb)));
    r.vx = vr.x;
    r.vy = vr.y;
    r.vz = vr.z;
    return r;
}

Expr *ExprQuaternion::Magnitude() const {
    return ((w ->Square())->Plus(
            (vx->Square())->Plus(
            (vy->Square())->Plus(
            (vz->Square())))))->Sqrt();
}

bool ExprQuaternion::Resolve(const ResolutionMap &resolutions) {
    return vx->Resolve(resolutions) && vy->Resolve(resolutions) &&
        vz->Resolve(resolutions) && w->Resolve(resolutions);
}

Expr *Expr::From(hParam p) {
    Expr *r = AllocExpr();
    r->op = Op::PARAM;
    r->parh = p;
    return r;
}

Expr *Expr::From(double v) {
    // Statically allocate common constants.
    // Note: this is only valid because AllocExpr() uses AllocTemporary(),
    // and Expr* is never explicitly freed.

    if(v == 0.0) {
        static Expr zero(0.0);
        return &zero;
    }

    if(v == 1.0) {
        static Expr one(1.0);
        return &one;
    }

    if(v == -1.0) {
        static Expr mone(-1.0);
        return &mone;
    }

    if(v == 0.5) {
        static Expr half(0.5);
        return &half;
    }

    if(v == -0.5) {
        static Expr mhalf(-0.5);
        return &mhalf;
    }

    Expr *r = AllocExpr();
    r->op = Op::CONSTANT;
    r->v = v;
    return r;
}

Expr *Expr::AnyOp(Op newOp, Expr *b) {
    Expr *r = AllocExpr();
    r->op = newOp;
    r->a = this;
    r->b = b;
    return r;
}

int Expr::Children() const {
    switch(op) {
        case Op::PARAM:
        case Op::PARAM_PTR:
        case Op::CONSTANT:
        case Op::VARIABLE:
            return 0;

        case Op::PLUS:
        case Op::MINUS:
        case Op::TIMES:
        case Op::DIV:
            return 2;

        case Op::NEGATE:
        case Op::SQRT:
        case Op::SQUARE:
        case Op::SIN:
        case Op::COS:
        case Op::ASIN:
        case Op::ACOS:
            return 1;
    }
    ssassert(false, "Unexpected operation");
}

int Expr::Nodes() const {
    switch(Children()) {
        case 0: return 1;
        case 1: return 1 + a->Nodes();
        case 2: return 1 + a->Nodes() + b->Nodes();
        default: ssassert(false, "Unexpected children count");
    }
}

Expr *Expr::DeepCopyWithParamsAsPointers(ParamList *firstTry, ParamList *thenTry,
                                         const ResolutionMap &resolutions,
                                         bool foldConstants) const {
    hParam eparh = op == Op::PARAM ? parh : NO_PARAMS;
    Op eop = this->op;
    if(eop == Op::VARIABLE) {
        auto it = resolutions.find(s);
        ssassert(it != resolutions.end(), "Failed to resolve named parameter");
        eparh = it->second;
        eop  = Op::PARAM;
    }
    Expr *n = AllocExpr();
    if(eop == Op::PARAM) {
        // A param that is referenced by its hParam gets rewritten to go
        // straight in to the parameter table with a pointer, or simply
        // into a constant if it's already known.
        Param *p = firstTry->FindByIdNoOops(eparh);
        if(!p) p = thenTry->FindById(eparh);
        if(p->known) {
            n->op = Op::CONSTANT;
            n->v = p->val;
        } else {
            n->op = Op::PARAM_PTR;
            n->parp = p;
        }
        return n;
    }

    *n = *this;
    int c = n->Children();
    if(c > 0) {
        n->a = a->DeepCopyWithParamsAsPointers(firstTry, thenTry, resolutions, foldConstants);
        bool hasConstants = n->a->op == Op::CONSTANT;
        if(c > 1) {
            n->b = b->DeepCopyWithParamsAsPointers(firstTry, thenTry, resolutions, foldConstants);
            hasConstants |= n->b->op == Op::CONSTANT;
        }
        if(hasConstants && foldConstants) {
            n = n->FoldConstants(false, 0);
        }
    }
    return n;
}

double Expr::Eval(const ResolutionMap &resolutions) const {
    switch(op) {
        case Op::PARAM:         return SK.GetParam(parh)->val;
        case Op::PARAM_PTR:     return parp->val;

        case Op::CONSTANT:      return v;
        case Op::VARIABLE:      {
            auto it = resolutions.find(s);
            ssassert(it != resolutions.end(), "Failed to resolve named parameter during Eval()");
            return SK.GetParam(it->second)->val;
        }

        case Op::PLUS:          return a->Eval(resolutions) + b->Eval(resolutions);
        case Op::MINUS:         return a->Eval(resolutions) - b->Eval(resolutions);
        case Op::TIMES:         return a->Eval(resolutions) * b->Eval(resolutions);
        case Op::DIV:           return a->Eval(resolutions) / b->Eval(resolutions);

        case Op::NEGATE:        return -(a->Eval(resolutions));
        case Op::SQRT:          return sqrt(a->Eval(resolutions));
        case Op::SQUARE:        { double r = a->Eval(resolutions); return r*r; }
        case Op::SIN:           return sin(a->Eval(resolutions));
        case Op::COS:           return cos(a->Eval(resolutions));
        case Op::ACOS:          return acos(a->Eval(resolutions));
        case Op::ASIN:          return asin(a->Eval(resolutions));
    }
    ssassert(false, "Unexpected operation");
}

Expr *Expr::PartialWrt(hParam p) const {
    Expr *da, *db;

    switch(op) {
        case Op::PARAM_PTR: return From(p == parp->h ? 1 : 0);
        case Op::PARAM:     return From(p == parh ? 1 : 0);

        case Op::CONSTANT:  return From(0.0);
        case Op::VARIABLE:  ssassert(false, "Variables should be resolved before calling PartialWrt");

        case Op::PLUS:      return (a->PartialWrt(p))->Plus(b->PartialWrt(p));
        case Op::MINUS:     return (a->PartialWrt(p))->Minus(b->PartialWrt(p));

        case Op::TIMES:
            da = a->PartialWrt(p);
            db = b->PartialWrt(p);
            return (a->Times(db))->Plus(b->Times(da));

        case Op::DIV:
            da = a->PartialWrt(p);
            db = b->PartialWrt(p);
            return ((da->Times(b))->Minus(a->Times(db)))->Div(b->Square());

        case Op::SQRT:
            return (From(0.5)->Div(a->Sqrt()))->Times(a->PartialWrt(p));

        case Op::SQUARE:
            return (From(2.0)->Times(a))->Times(a->PartialWrt(p));

        case Op::NEGATE:    return (a->PartialWrt(p))->Negate();
        case Op::SIN:       return (a->Cos())->Times(a->PartialWrt(p));
        case Op::COS:       return ((a->Sin())->Times(a->PartialWrt(p)))->Negate();

        case Op::ASIN:
            return (From(1)->Div((From(1)->Minus(a->Square()))->Sqrt()))
                        ->Times(a->PartialWrt(p));
        case Op::ACOS:
            return (From(-1)->Div((From(1)->Minus(a->Square()))->Sqrt()))
                        ->Times(a->PartialWrt(p));
    }
    ssassert(false, "Unexpected operation");
}

void Expr::ParamsUsedList(ParamSet *list, const ResolutionMap &resolutions) const {
    if(op == Op::VARIABLE) {
        auto it = resolutions.find(s);
        ssassert(it != resolutions.end(), "Unresolved variable reference");
        list->insert(it->second);
        return;
    }
    if(op == Op::PARAM || op == Op::PARAM_PTR) {
        // leaf: just add ourselves if we aren't already there
        hParam param = (op == Op::PARAM) ? parh : parp->h;
        list->insert(param);
        return;
    }

    int c = Children();
    if(c >= 1) {
        a->ParamsUsedList(list);
        if(c >= 2) b->ParamsUsedList(list);
    }
}

bool Expr::Tol(double a, double b) {
    return fabs(a - b) < 0.001;
}

bool Expr::IsZeroConst() const {
    return op == Op::CONSTANT && EXACT(v == 0.0);
}

Expr *Expr::FoldConstants(bool allocCopy, size_t depth) {
    Expr *n;
    
    if(!allocCopy) {
        n = this;
    } else {
        n = AllocExpr();
        *n = *this;
    }

    switch(op) {
        case Op::PARAM_PTR:
        case Op::PARAM:
        case Op::CONSTANT:
        case Op::VARIABLE:
            break;

        case Op::MINUS:
        case Op::TIMES:
        case Op::DIV:
        case Op::PLUS:
            if(depth > 0) {
                n->a = a->FoldConstants(allocCopy, depth - 1);
                n->b = b->FoldConstants(allocCopy, depth - 1);
            }
    
            // If both ops are known, then we can evaluate immediately
            if(n->a->op == Op::CONSTANT && n->b->op == Op::CONSTANT) {
                double nv = n->Eval();
                n->op = Op::CONSTANT;
                n->v = nv;
                break;
            }
            // x + 0 = 0 + x = x
            if(op == Op::PLUS && n->b->op == Op::CONSTANT && Tol(n->b->v, 0)) {
                *n = *(n->a); break;
            }
            if(op == Op::PLUS && n->a->op == Op::CONSTANT && Tol(n->a->v, 0)) {
                *n = *(n->b); break;
            }
            // 1*x = x*1 = x
            if(op == Op::TIMES && n->b->op == Op::CONSTANT && Tol(n->b->v, 1)) {
                *n = *(n->a); break;
            }
            if(op == Op::TIMES && n->a->op == Op::CONSTANT && Tol(n->a->v, 1)) {
                *n = *(n->b); break;
            }
            // 0*x = x*0 = 0
            if(op == Op::TIMES && n->b->op == Op::CONSTANT && Tol(n->b->v, 0)) {
                n->op = Op::CONSTANT; n->v = 0; break;
            }
            if(op == Op::TIMES && n->a->op == Op::CONSTANT && Tol(n->a->v, 0)) {
                n->op = Op::CONSTANT; n->v = 0; break;
            }

            break;

        case Op::SQRT:
        case Op::SQUARE:
        case Op::NEGATE:
        case Op::SIN:
        case Op::COS:
        case Op::ASIN:
        case Op::ACOS:
            if(depth > 0) {
                n->a = a->FoldConstants(allocCopy, depth - 1);
            }

            if(n->a->op == Op::CONSTANT) {
                double nv = n->Eval();
                n->op = Op::CONSTANT;
                n->v = nv;
            }
            break;
    }
    return n;
}

void Expr::Substitute(const SubstitutionMap &subMap) {
    ssassert(op != Op::PARAM_PTR, "Expected an expression that refer to params via handles");

    if(op == Op::PARAM) {
        auto it = subMap.find(parh);
        if(it != subMap.end()) {
            parh = it->second->h;
        }
    } else {
        int c = Children();
        if(c >= 1) {
            a->Substitute(subMap);
            if(c >= 2) b->Substitute(subMap);
        }
    }
}

bool Expr::Resolve(const ResolutionMap &resolutions) {
    if(op == Op::VARIABLE) {
        auto it = resolutions.find(s);
        if(it == resolutions.end()) {
            return false;
        }
        parh = it->second;
        op   = Op::PARAM;
    } else {
        const int c = Children();
        if(c >= 1) {
            if(!a->Resolve(resolutions)) {
                return false;
            }
            if(c >= 2 && !b->Resolve(resolutions)) {
                return false;
            }
        }
    }

    return true;
}

//-----------------------------------------------------------------------------
// If the expression references only one parameter that appears in pl, then
// return that parameter. If no param is referenced, then return NO_PARAMS.
// If multiple params are referenced, then return MULTIPLE_PARAMS.
//-----------------------------------------------------------------------------
const hParam Expr::NO_PARAMS       = { 0 };
const hParam Expr::MULTIPLE_PARAMS = { std::numeric_limits<decltype(hParam::v)>::max() };
hParam Expr::ReferencedParams(ParamList *pl, const ResolutionMap &resolutions) const {
    if(op == Op::VARIABLE) {
        auto it = resolutions.find(s);
        if(it == resolutions.end()) {
            return MULTIPLE_PARAMS;
        }
        return it->second;
    }
    if(op == Op::PARAM) {
        if(pl->FindByIdNoOops(parh)) {
            return parh;
        } else {
            return NO_PARAMS;
        }
    }
    ssassert(op != Op::PARAM_PTR, "Expected an expression that refer to params via handles");

    int c = Children();
    if(c == 0) {
        return NO_PARAMS;
    } else if(c == 1) {
        return a->ReferencedParams(pl, resolutions);
    } else if(c == 2) {
        hParam pa, pb;
        pa = a->ReferencedParams(pl, resolutions);
        pb = b->ReferencedParams(pl, resolutions);
        if(pa == NO_PARAMS) {
            return pb;
        } else if(pb == NO_PARAMS) {
            return pa;
        } else if(pa == pb) {
            return pa; // either, doesn't matter
        } else {
            return MULTIPLE_PARAMS;
        }
    } else ssassert(false, "Unexpected children count");
}


//-----------------------------------------------------------------------------
// Routines to pretty-print an expression. Mostly for debugging.
//-----------------------------------------------------------------------------

std::string Expr::Print() const {
    char c;
    switch(op) {
        case Op::PARAM:     return ssprintf("param(%08x)", parh.v);
        case Op::PARAM_PTR: return ssprintf("param(p%08x)", parp->h.v);

        case Op::CONSTANT:  return ssprintf("%.3f", v);
        case Op::VARIABLE:  return s;

        case Op::PLUS:      c = '+'; goto p;
        case Op::MINUS:     c = '-'; goto p;
        case Op::TIMES:     c = '*'; goto p;
        case Op::DIV:       c = '/'; goto p;
p:
            return "(" + a->Print() + " " + c + " " + b->Print() + ")";
            break;

        case Op::NEGATE:    return "(- " + a->Print() + ")";
        case Op::SQRT:      return "(sqrt " + a->Print() + ")";
        case Op::SQUARE:    return "(square " + a->Print() + ")";
        case Op::SIN:       return "(sin " + a->Print() + ")";
        case Op::COS:       return "(cos " + a->Print() + ")";
        case Op::ASIN:      return "(asin " + a->Print() + ")";
        case Op::ACOS:      return "(acos " + a->Print() + ")";
    }
    ssassert(false, "Unexpected operation");
}


//-----------------------------------------------------------------------------
// A parser; convert a string to an expression. Infix notation, with the
// usual shift/reduce approach. I had great hopes for user-entered eq
// constraints, but those don't seem very useful, so right now this is just
// to provide calculator type functionality wherever numbers are entered.
//-----------------------------------------------------------------------------

class ExprParser {
public:
    enum class TokenType {
        ERROR = 0,

        PAREN_LEFT,
        PAREN_RIGHT,
        BINARY_OP,
        UNARY_OP,
        OPERAND,

        END,
    };

    class Token {
    public:
        TokenType  type;
        Expr      *expr;

        static Token From(double constant);
        static Token From(TokenType type = TokenType::ERROR, Expr *expr = NULL);
        static Token From(TokenType type, Expr::Op op);
        bool IsError() const { return type == TokenType::ERROR; }
    };

    std::string::const_iterator it, end;
    std::vector<Token> stack;

    char ReadChar();
    char PeekChar();

    std::string ReadWord();
    void SkipSpace();

    Token PopOperator(std::string *error);
    Token PopOperand(std::string *error);

    int Precedence(Token token);
    Token LexNumber(std::string *error);
    Token Lex(bool allowVariables, std::string *error);
    bool Reduce(std::string *error);
    bool Parse(bool allowVariables, std::unordered_set<std::string> *variableRefs,
               std::string *error, size_t reduceUntil = 0);

    static Expr *Parse(const std::string &input, bool allowVariables,
                       std::unordered_set<std::string> *variableRefs,
                       std::string *error);
};

ExprParser::Token ExprParser::Token::From(double constant) {
    Token t = Token::From(TokenType::OPERAND, Expr::Op::CONSTANT);
    t.expr->v = constant;
    return t;
}

ExprParser::Token ExprParser::Token::From(TokenType type, Expr *expr) {
    Token t;
    t.type = type;
    t.expr = expr;
    return t;
}

ExprParser::Token ExprParser::Token::From(TokenType type, Expr::Op op) {
    Token t;
    t.type = type;
    t.expr = AllocExpr();
    t.expr->op = op;
    return t;
}

char ExprParser::ReadChar() {
    return *it++;
}

char ExprParser::PeekChar() {
    if(it == end) {
        return '\0';
    } else {
        return *it;
    }
}

std::string ExprParser::ReadWord() {
    std::string s;

    while(char c = PeekChar()) {
        if(!isalnum(c) && c != '_')
            break;
        s.push_back(ReadChar());
    }

    return s;
}

void ExprParser::SkipSpace() {
    while(char c = PeekChar()) {
        if(!isspace(c)) break;
        ReadChar();
    }
}

ExprParser::Token ExprParser::LexNumber(std::string *error) {
    std::string s;

    while(char c = PeekChar()) {
        if(!((c >= '0' && c <= '9') || c == 'e' || c == 'E' || c == '.' || c == '_')) break;
        if(c == '_') {
            ReadChar();
            continue;
        }
        s.push_back(ReadChar());
    }

    char *endptr;
    double d = strtod(s.c_str(), &endptr);

    Token t = Token::From();
    if(endptr == s.c_str() + s.size()) {
        t = Token::From(TokenType::OPERAND, Expr::Op::CONSTANT);
        t.expr->v = d;
    } else {
        *error = "'" + s + "' is not a valid number";
    }
    return t;
}

static const std::unordered_map<std::string, ExprParser::Token (*)()> KNOWN_IDENTIFIERS{
    {"sqrt", [] { return ExprParser::Token::From(ExprParser::TokenType::UNARY_OP, Expr::Op::SQRT); }},
    {"square", [] { return ExprParser::Token::From(ExprParser::TokenType::UNARY_OP, Expr::Op::SQUARE); }},
    {"sin", [] { return ExprParser::Token::From(ExprParser::TokenType::UNARY_OP, Expr::Op::SIN); }},
    {"cos", [] { return ExprParser::Token::From(ExprParser::TokenType::UNARY_OP, Expr::Op::COS); }},
    {"asin", [] { return ExprParser::Token::From(ExprParser::TokenType::UNARY_OP, Expr::Op::ASIN); }},
    {"acos", [] { return ExprParser::Token::From(ExprParser::TokenType::UNARY_OP, Expr::Op::ACOS); }},
    {"pi", [] { return ExprParser::Token::From(PI); }},
};

static const std::unordered_set<std::string> RESERVED_IDENTIFIERS{
    // Useful functions that we might want to implement at some point
    "neg",
    "abs",
    "sgn",
    "mod",
    "ceil",
    "floor",
    "log",
    "log2",
    "log10",
    "ln",
    "lb",
    "lg",
    "atan",
    "atan2",
    // Useful constants that we might want to introduce later
    "e", "tau", "phi",
};

ExprParser::Token ExprParser::Lex(bool allowVariables, std::string *error) {
    SkipSpace();

    Token t = Token::From();
    char c = PeekChar();
    if(isalpha(c)) {
        std::string s = ReadWord();
        auto it       = KNOWN_IDENTIFIERS.find(s);
        if(it == KNOWN_IDENTIFIERS.end()) {
            if(RESERVED_IDENTIFIERS.find(s) != RESERVED_IDENTIFIERS.end()) {
                *error = "'" + s + "' is reserved and cannot be used to name a variable";
            } else if(!allowVariables) {
                *error = "Referencing variable '" + s + "' is not allowed in this expression";
            } else {
                t = Token::From(TokenType::OPERAND, Expr::Op::VARIABLE);
                char *cstr = (char *)Platform::AllocTemporary(s.size() + 1);
                memcpy(cstr, s.c_str(), s.size() + 1);
                t.expr->s = cstr;
            }
        } else {
            t = (*it->second)();
        }
    } else if(isdigit(c) || c == '.') {
        return LexNumber(error);
    } else if(ispunct(c)) {
        ReadChar();
        if(c == '+') {
            t = Token::From(TokenType::BINARY_OP, Expr::Op::PLUS);
        } else if(c == '-') {
            t = Token::From(TokenType::BINARY_OP, Expr::Op::MINUS);
        } else if(c == '*') {
            t = Token::From(TokenType::BINARY_OP, Expr::Op::TIMES);
        } else if(c == '/') {
            t = Token::From(TokenType::BINARY_OP, Expr::Op::DIV);
        } else if(c == '(') {
            t = Token::From(TokenType::PAREN_LEFT);
        } else if(c == ')') {
            t = Token::From(TokenType::PAREN_RIGHT);
        } else {
            *error = "'" + std::string(1, c) + "' is not a valid operator";
        }
    } else if(c == '\0') {
        t = Token::From(TokenType::END);
    } else {
        *error = "Unexpected character '" + std::string(1, c) + "'";
    }

    return t;
}

ExprParser::Token ExprParser::PopOperand(std::string *error) {
    Token t = Token::From();
    if(stack.empty() || stack.back().type != TokenType::OPERAND) {
        *error = "Expected an operand";
    } else {
        t = stack.back();
        stack.pop_back();
    }
    return t;
}

ExprParser::Token ExprParser::PopOperator(std::string *error) {
    Token t = Token::From();
    if(stack.empty() || (stack.back().type != TokenType::UNARY_OP &&
                         stack.back().type != TokenType::BINARY_OP)) {
        *error = "Expected an operator";
    } else {
        t = stack.back();
        stack.pop_back();
    }
    return t;
}

int ExprParser::Precedence(Token t) {
    ssassert(t.type == TokenType::BINARY_OP ||
             t.type == TokenType::UNARY_OP ||
             t.type == TokenType::OPERAND,
             "Unexpected token type");

    if(t.type == TokenType::UNARY_OP) {
        return 30;
    } else if(t.expr->op == Expr::Op::TIMES ||
              t.expr->op == Expr::Op::DIV) {
        return 20;
    } else if(t.expr->op == Expr::Op::PLUS ||
              t.expr->op == Expr::Op::MINUS) {
        return 10;
    } else if(t.type == TokenType::OPERAND) {
        return 0;
    } else ssassert(false, "Unexpected operator");
}

bool ExprParser::Reduce(std::string *error) {
    Token a = PopOperand(error);
    if(a.IsError()) return false;

    Token op = PopOperator(error);
    if(op.IsError()) return false;

    Token r = Token::From(TokenType::OPERAND);
    switch(op.type) {
        case TokenType::BINARY_OP: {
            Token b = PopOperand(error);
            if(b.IsError()) return false;
            r.expr = b.expr->AnyOp(op.expr->op, a.expr);
            break;
        }

        case TokenType::UNARY_OP: {
            Expr *e = a.expr;
            switch(op.expr->op) {
                case Expr::Op::NEGATE: e = e->Negate(); break;
                case Expr::Op::SQRT:   e = e->Sqrt(); break;
                case Expr::Op::SQUARE: e = e->Times(e); break;
                case Expr::Op::SIN:    e = e->Times(Expr::From(PI/180))->Sin(); break;
                case Expr::Op::COS:    e = e->Times(Expr::From(PI/180))->Cos(); break;
                case Expr::Op::ASIN:   e = e->ASin()->Times(Expr::From(180/PI)); break;
                case Expr::Op::ACOS:   e = e->ACos()->Times(Expr::From(180/PI)); break;
                default: ssassert(false, "Unexpected unary operator");
            }
            r.expr = e;
            break;
        }

        default: ssassert(false, "Unexpected operator");
    }
    stack.push_back(r);

    return true;
}

bool ExprParser::Parse(bool allowVariables, std::unordered_set<std::string> *variableRefs,
                       std::string *error, size_t reduceUntil) {
    while(true) {
        Token t = Lex(allowVariables, error);
        switch(t.type) {
            case TokenType::ERROR:
                return false;

            case TokenType::END:
            case TokenType::PAREN_RIGHT:
                while(stack.size() > 1 + reduceUntil) {
                    if(!Reduce(error)) return false;
                }

                if(t.type == TokenType::PAREN_RIGHT) {
                    stack.push_back(t);
                }
                return true;

            case TokenType::PAREN_LEFT: {
                // sub-expression
                if(!Parse(allowVariables, variableRefs, error, /*reduceUntil=*/stack.size())) return false;

                if(stack.empty() || stack.back().type != TokenType::PAREN_RIGHT) {
                    *error = "Expected ')'";
                    return false;
                }
                stack.pop_back();
                break;
            }

            case TokenType::BINARY_OP:
                if((stack.size() > reduceUntil && stack.back().type != TokenType::OPERAND) ||
                   stack.size() == reduceUntil) {
                    if(t.expr->op == Expr::Op::MINUS) {
                        t.type = TokenType::UNARY_OP;
                        t.expr->op = Expr::Op::NEGATE;
                        stack.push_back(t);
                        break;
                    }
                }

                while(stack.size() > 1 + reduceUntil &&
                      Precedence(t) <= Precedence(stack[stack.size() - 2])) {
                    if(!Reduce(error)) return false;
                }

                stack.push_back(t);
                break;

            case TokenType::UNARY_OP:
            case TokenType::OPERAND:
                if(t.expr->op == Expr::Op::VARIABLE) {
                    variableRefs->insert(t.expr->s);
                }
                stack.push_back(t);
                break;
        }
    }

    return true;
}

Expr *ExprParser::Parse(const std::string &input, bool allowVariables,
                        std::unordered_set<std::string> *variableRefs, std::string *error) {
    ExprParser parser;
    parser.it  = input.cbegin();
    parser.end = input.cend();
    std::unordered_set<std::string> varRefs;
    if(!parser.Parse(allowVariables, &varRefs, error))
        return NULL;

    Token r = parser.PopOperand(error);
    if(r.IsError()) return NULL;
    if(variableRefs) {
        *variableRefs = std::move(varRefs);
    }
    return r.expr;
}

Expr *Expr::Parse(const std::string &input, bool allowVariables,
                  std::unordered_set<std::string> *variableRefs, std::string *error) {
    return ExprParser::Parse(input, allowVariables, variableRefs, error);
}

Expr *Expr::From(const std::string &input, bool allowVariables, bool popUpError,
                 std::unordered_set<std::string> *variableRefs) {
    std::string error;
    Expr *e = ExprParser::Parse(input, allowVariables, variableRefs, &error);
    if(!e) {
        dbp("Parse/lex error: %s", error.c_str());
        if(popUpError) {
            Error("Not a valid number or expression: '%s'.\n%s.",
                  input.c_str(), error.c_str());
        }
    }
    return e;
}

} // namespace SolveSpace
