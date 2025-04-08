//-----------------------------------------------------------------------------
// Data structures used frequently in the program, various kinds of vectors
// (of real numbers, not symbolic algebra stuff) and our templated lists.
//
// Copyright 2008-2013 Jonathan Westhues.
//-----------------------------------------------------------------------------
#ifndef SOLVESPACE_DSC_H
#define SOLVESPACE_DSC_H

#include <algorithm>
#include <limits>
#include <new>
#include <type_traits>
#include <vector>

/// Trait indicating which types are handle types and should get the associated operators.
/// Specialize for each handle type and inherit from std::true_type.
template<typename T>
struct IsHandleOracle : std::false_type {};

// Equality-compare any two instances of a handle type.
template<typename T>
static inline typename std::enable_if<IsHandleOracle<T>::value, bool>::type
operator==(T const &lhs, T const &rhs) {
    return lhs.v == rhs.v;
}

// Inequality-compare any two instances of a handle type.
template<typename T>
static inline typename std::enable_if<IsHandleOracle<T>::value, bool>::type
operator!=(T const &lhs, T const &rhs) {
    return !(lhs == rhs);
}

// Less-than-compare any two instances of a handle type.
template<typename T>
static inline typename std::enable_if<IsHandleOracle<T>::value, bool>::type
operator<(T const &lhs, T const &rhs) {
    return lhs.v < rhs.v;
}

class Vector;
class Vector4;
class Point2d;
class hEntity;
class hParam;

class Quaternion {
public:
    // a + (vx)*i + (vy)*j + (vz)*k
    double w, vx, vy, vz;

    static const Quaternion IDENTITY;

    static Quaternion From(double w, double vx, double vy, double vz);
    static Quaternion From(hParam w, hParam vx, hParam vy, hParam vz);
    static Quaternion From(Vector u, Vector v);
    static Quaternion From(Vector axis, double dtheta);

    Quaternion Plus(Quaternion b) const;
    Quaternion Minus(Quaternion b) const;
    Quaternion ScaledBy(double s) const;
    double Magnitude() const;
    Quaternion WithMagnitude(double s) const;

    // Call a rotation matrix [ u' v' n' ]'; this returns the first and
    // second rows, where that matrix is generated by this quaternion
    Vector RotationU() const;
    Vector RotationV() const;
    Vector RotationN() const;
    Vector Rotate(Vector p) const;

    Quaternion ToThe(double p) const;
    Quaternion Inverse() const;
    Quaternion Times(Quaternion b) const;
    Quaternion Mirror() const;
};

class Vector {
public:
    double x, y, z;

    static Vector From(double x, double y, double z);
    static Vector From(hParam x, hParam y, hParam z);
    static Vector AtIntersectionOfPlanes(Vector n1, double d1,
                                         Vector n2, double d2);
    static Vector AtIntersectionOfLines(Vector a0, Vector a1,
                                        Vector b0, Vector b1,
                                        bool *skew,
                                        double *pa=NULL, double *pb=NULL);
    static Vector AtIntersectionOfPlaneAndLine(Vector n, double d,
                                               Vector p0, Vector p1,
                                               bool *parallel);
    static Vector AtIntersectionOfPlanes(Vector na, double da,
                                         Vector nb, double db,
                                         Vector nc, double dc, bool *parallel);
    static void ClosestPointBetweenLines(Vector pa, Vector da,
                                         Vector pb, Vector db,
                                         double *ta, double *tb);

    double Element(int i) const;
    bool Equals(Vector v, double tol=LENGTH_EPS) const;
    bool EqualsExactly(Vector v) const;
    Vector Plus(Vector b) const;
    Vector Minus(Vector b) const;
    Vector Negated() const;
    Vector Cross(Vector b) const;
    double DirectionCosineWith(Vector b) const;
    double Dot(Vector b) const;
    Vector Normal(int which) const;
    Vector RotatedAbout(Vector orig, Vector axis, double theta) const;
    Vector RotatedAbout(Vector axis, double theta) const;
    Vector DotInToCsys(Vector u, Vector v, Vector n) const;
    Vector ScaleOutOfCsys(Vector u, Vector v, Vector n) const;
    double DistanceToLine(Vector p0, Vector dp) const;
    double DistanceToPlane(Vector normal, Vector origin) const;
    bool OnLineSegment(Vector a, Vector b, double tol=LENGTH_EPS) const;
    Vector ClosestPointOnLine(Vector p0, Vector deltal) const;
    double Magnitude() const;
    double MagSquared() const;
    Vector WithMagnitude(double s) const;
    Vector ScaledBy(double s) const;
    Vector ProjectInto(hEntity wrkpl) const;
    Vector ProjectVectorInto(hEntity wrkpl) const;
    double DivProjected(Vector delta) const;
    Vector ClosestOrtho() const;
    void MakeMaxMin(Vector *maxv, Vector *minv) const;
    Vector ClampWithin(double minv, double maxv) const;
    static bool BoundingBoxesDisjoint(Vector amax, Vector amin,
                                      Vector bmax, Vector bmin);
    static bool BoundingBoxIntersectsLine(Vector amax, Vector amin,
                                          Vector p0, Vector p1, bool asSegment);
    bool OutsideAndNotOn(Vector maxv, Vector minv) const;
    Vector InPerspective(Vector u, Vector v, Vector n,
                         Vector origin, double cameraTan) const;
    Point2d Project2d(Vector u, Vector v) const;
    Point2d ProjectXy() const;
    Vector4 Project4d() const;
};

inline double Vector::Element(int i) const {
    switch (i) {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default: ssassert(false, "Unexpected vector element index");
    }
}

inline bool Vector::Equals(Vector v, double tol) const {
    // Quick axis-aligned tests before going further
    const Vector dv = this->Minus(v);
    if (fabs(dv.x) > tol) return false;
    if (fabs(dv.y) > tol) return false;
    if (fabs(dv.z) > tol) return false;

    return dv.MagSquared() < tol*tol;
}

inline Vector Vector::From(double x, double y, double z) {
    return {x, y, z};
}

inline Vector Vector::Plus(Vector b) const {
    return {x + b.x, y + b.y, z + b.z};
}

inline Vector Vector::Minus(Vector b) const {
    return {x - b.x, y - b.y, z - b.z};
}

inline Vector Vector::Negated() const {
    return {-x, -y, -z};
}

inline Vector Vector::Cross(Vector b) const {
    return {-(z * b.y) + (y * b.z), (z * b.x) - (x * b.z), -(y * b.x) + (x * b.y)};
}

inline double Vector::Dot(Vector b) const {
    return (x * b.x + y * b.y + z * b.z);
}

inline double Vector::MagSquared() const {
    return x * x + y * y + z * z;
}

inline double Vector::Magnitude() const {
    return sqrt(x * x + y * y + z * z);
}

inline Vector Vector::ScaledBy(const double v) const {
    return {x * v, y * v, z * v};
}

inline void Vector::MakeMaxMin(Vector *maxv, Vector *minv) const {
    maxv->x = max(maxv->x, x);
    maxv->y = max(maxv->y, y);
    maxv->z = max(maxv->z, z);

    minv->x = min(minv->x, x);
    minv->y = min(minv->y, y);
    minv->z = min(minv->z, z);
}

struct VectorHash {
    size_t operator()(const Vector &v) const;
};

struct VectorPred {
    bool operator()(Vector a, Vector b) const;
};

class Vector4 {
public:
    double w, x, y, z;

    static Vector4 From(double w, double x, double y, double z);
    static Vector4 From(double w, Vector v3);
    static Vector4 Blend(Vector4 a, Vector4 b, double t);

    Vector4 Plus(Vector4 b) const;
    Vector4 Minus(Vector4 b) const;
    Vector4 ScaledBy(double s) const;
    Vector PerspectiveProject() const;
};

class Point2d {
public:
    double x, y;

    static Point2d From(double x, double y);
    static Point2d FromPolar(double r, double a);

    Point2d Plus(const Point2d &b) const;
    Point2d Minus(const Point2d &b) const;
    Point2d ScaledBy(double s) const;
    double DivProjected(Point2d delta) const;
    double Dot(Point2d p) const;
    double DistanceTo(const Point2d &p) const;
    double DistanceToLine(const Point2d &p0, const Point2d &dp, bool asSegment) const;
    double DistanceToLineSigned(const Point2d &p0, const Point2d &dp, bool asSegment) const;
    double Angle() const;
    double AngleTo(const Point2d &p) const;
    double Magnitude() const;
    double MagSquared() const;
    Point2d WithMagnitude(double v) const;
    Point2d Normal() const;
    bool Equals(Point2d v, double tol=LENGTH_EPS) const;
    bool LessThan(const Point2d &v, double tol=LENGTH_EPS) const;
};

// A simple list
template<class T>
class List {
    std::vector<T> elem;

public:
    bool IsEmpty() const {
        return elem.empty();
    }

    int Size() const {
        return (int)elem.size();
    }

    void ReserveMore(int howMuch) {
        if(howMuch > 0) {
            elem.reserve(elem.capacity() + howMuch);
        }
    }

    void Add(const T &t) {
        elem.push_back(t);
    }

    void Add(T &&t) {
        elem.emplace_back(std::forward<T>(t));
    }

    void AddToBeginning(const T *t) {
        elem.insert(elem.begin(), *t);
    }

    T *First() {
        return IsEmpty() ? nullptr : &elem.front();
    }
    const T *First() const {
        return IsEmpty() ? nullptr : &elem.front();
    }

    T *Last() {
        return IsEmpty() ? nullptr : &elem.back();
    }
    const T *Last() const {
        return IsEmpty() ? nullptr : &elem.back();
    }

    T &Get(size_t i) { return elem[i]; }
    T const &Get(size_t i) const { return elem[i]; }
    T &operator[](size_t i) { return Get(i); }
    T const &operator[](size_t i) const { return Get(i); }

    T *begin() { return elem.data(); }
    T *end() { return elem.data() + Size(); }
    const T *begin() const { return elem.data(); }
    const T *end() const { return elem.data() + Size(); }
    const T *cbegin() const { return begin(); }
    const T *cend() const { return end(); }

    void ClearTags() {
        for(auto &elt : *this) {
            elt.tag = 0;
        }
    }

    void Clear() {
        elem.clear();
    }

    void RemoveTagged() {
        auto newEnd = std::remove_if(this->begin(), this->end(), [](const T &t) {
            return t.tag != 0;
        });
        elem.resize(newEnd - begin());
    }

    void RemoveLast(int cnt) {
        ssassert(Size() >= cnt, "Removing more elements than the list contains");
        elem.resize(Size() - cnt);
    }

    void Reverse() {
        std::reverse(elem.begin(), elem.end());
    }
};

// A list, where each element has an integer identifier. The list is kept
// sorted by that identifier, and items can be looked up in log n time by
// id.
template<class T>
class IdList {
    static_assert(!std::is_reference<T>::value, "Cannot store reference types");

    using Handle      = decltype(std::remove_pointer<T>::type::h);
    using HandleValue = decltype(Handle::v);

    // XXX: Can't use this assertion because `IdList<Canvas::Stroke>` and `IdList<Canvas::Fill>`
    //      are instantiated inside `Canvas`, before `IsHandleOracle` can be specialised for
    //      `Canvas::hStroke` and `Canvas::hFill`.
    // static_assert(IsHandleOracle<Handle>::value, "Invalid handle type");
    static_assert(std::is_integral<HandleValue>::value && sizeof(HandleValue) <= sizeof(size_t),
                  "Invalid handle value type");

    struct Storage {
        Storage() : used_(false) {}

        explicit Storage(const T &t) : used_(true) {
            new(&data) T(t);
        }

        explicit Storage(T &&t) : used_(true) {
            new(&data) T(std::forward<T>(t));
        }

        Storage(const Storage &other) : used_(other.used()) {
            if(other.used()) {
                new(&data) T(*other.get());
            }
        }

        Storage(Storage &&other) noexcept(noexcept(T(std::move(std::declval<T>())))) : used_(other.used()) {
            if(other.used()) {
                new(&data) T(std::move(*other.get()));
                other.reset();
            }
        }

        ~Storage() {
            reset();
        }

        Storage &operator=(const Storage &other) {
            this->~Storage();
            return *new(this) Storage(other);
        }

        Storage &operator=(Storage &&other) noexcept(noexcept(Storage(std::forward<Storage>(other)))) {
            this->~Storage();
            return *new(this) Storage(std::forward<Storage>(other));
        }

        bool used() const noexcept {
            return used_;
        }

        void reset() {
            if(used()) {
                get()->~T();
                used_ = false;
            }
        }

        void reset(const T &t) {
            if(used()) {
                get()->~T();
            }
            new(&data) T(t);
            used_ = true;
        }

        void reset(T &&t) {
            if(used()) {
                get()->~T();
            }
            new(&data) T(std::forward<T>(t));
            used_ = true;
        }

        T *get() noexcept {
            return reinterpret_cast<T *>(&data);
        }
        const T *get() const noexcept {
            return reinterpret_cast<const T *>(&data);
        }

    private:
        typename std::aligned_storage<sizeof(T), alignof(T)>::type data;
        bool used_;
    };

    struct Target { uint32_t value; };

    std::vector<Storage> elemstore;
    std::vector<Target> targets;
    size_t used = 0;

public:
    // This iterator is only invalidated if items are removed or added before or at
    // its current position. References, however, can also be invalidated if the
    // underlying storage capacity changes.
    template<class ValueType>
    struct IdListIterator {
        using iterator_category = std::random_access_iterator_tag;
        using value_type        = ValueType;
        using difference_type   = std::ptrdiff_t;
        using pointer           = ValueType *;
        using reference         = ValueType &;

        using ListType = typename std::conditional<
            std::is_same<ValueType, T>::value, IdList<T>, const IdList<T>>::type;

    public:
        ValueType &operator*() noexcept {
            return *l->GetTarget(l->targets[idx]);
        }
        const ValueType &operator*() const noexcept {
            return *l->GetTarget(l->targets[idx]);
        }
        ValueType *operator->() noexcept {
            return l->GetTarget(l->targets[idx]);
        }
        const ValueType *operator->() const noexcept {
            return l->GetTarget(l->targets[idx]);
        }

        bool operator==(const IdListIterator &p) const noexcept {
            if(p.l != l) {
                ssassert(false, "Invalid comparison of iterators for different lists");
            }
            return p.idx == idx;
        }
        bool operator!=(const IdListIterator &p) const noexcept { return !operator==(p); }

        IdListIterator &operator++() noexcept {
            ++idx;
            return *this;
        }

        IdListIterator operator+(size_t i) const noexcept {
            return IdListIterator(l, idx + i);
        }

        // Needed for std::find_if of gcc used in entity.cpp GenerateEquations
        difference_type operator-(const IdListIterator &rhs) const noexcept {
            return idx - rhs.idx;
        }

        IdListIterator(ListType *l_, size_t idx_ = 0) noexcept : l(l_), idx(idx_) {}

    private:
        ListType *l;
        size_t idx;
    };

    using iterator = IdListIterator<T>;
    using const_iterator = IdListIterator<const T>;

    IdList() noexcept = default;

    IdList(IdList &&other) noexcept(noexcept(std::vector<Storage>(std::move(other.elemstore)))) :
            elemstore(std::move(other.elemstore)), targets(std::move(other.targets)), used(other.used) {
        other.used = 0;
    }

    IdList(const IdList &other) : elemstore(), targets(), used(other.Size()) {
        elemstore.reserve(used);
        targets.reserve(used);
        for(size_t i = 0; i < used; ++i) {
            elemstore.emplace_back(other.Get(i));
            targets.push_back({ uint32_t(i) });
        }
    }

    IdList &operator=(const IdList &other) {
        this->~IdList();
        return *new(this) IdList(other);
    }

    IdList &operator=(IdList &&other) noexcept(noexcept(IdList(std::forward<IdList>(other)))) {
        this->~IdList();
        return *new(this) IdList(std::forward<IdList>(other));
    }

    bool IsEmpty() const noexcept {
        return Size() == 0;
    }

    size_t Size() const noexcept {
        return used;
    }

    HandleValue MaximumId() const noexcept {
        if(IsEmpty()) {
            return 0;
        } else {
            return GetTarget(targets[Size()-1])->h.v;
        }
    }

    Handle AddAndAssignId(T &&t) {
        t.h.v = MaximumId() + 1;
        InsertAt(Size(), std::forward<T>(t));
        return t.h;
    }

    void ReserveMore(size_t howMuch) {
        const size_t total_reserve = Size() + howMuch;
        elemstore.reserve(total_reserve);
        targets.reserve(total_reserve);
    }

    void Add(const T &t) {
        const size_t idx = FindInsertionPoint(t.h);
        if(idx < Size()) {
            ssassert(GetTarget(targets[idx])->h.v != t.h.v, "Handle isn't unique");
        }
        InsertAt(idx, t);
    }

    void Add(T &&t) {
        const size_t idx = FindInsertionPoint(t.h);
        if(idx < Size()) {
            ssassert(GetTarget(targets[idx])->h.v != t.h.v, "Handle isn't unique");
        }
        InsertAt(idx, std::forward<T>(t));
    }

    T *FindById(Handle h) {
        T *t = FindByIdNoOops(h);
        ssassert(t != nullptr, "Cannot find handle");
        return t;
    }

    T *FindByIdNoOops(Handle h) noexcept {
        const size_t idx = FindInsertionPoint(h);
        if(idx >= Size()) {
            return nullptr;
        }
        T *t = GetTarget(targets[idx]);
        if(t->h.v != h.v) {
            return nullptr;
        }
        return t;
    }

    T &Get(size_t i) {
        ssassert(i < Size(), "Out of bounds access");
        return *GetTarget(targets[i]);
    }
    const T &Get(size_t i) const {
        ssassert(i < Size(), "Out of bounds access");
        return *GetTarget(targets[i]);
    }

    iterator begin() noexcept { return iterator(this); }
    iterator end() noexcept { return iterator(this, Size()); }
    const_iterator begin() const noexcept { return const_iterator(this); }
    const_iterator end() const noexcept { return const_iterator(this, Size()); }
    const_iterator cbegin() const noexcept { return const_iterator(this); }
    const_iterator cend() const noexcept { return const_iterator(this, Size()); }

    void ClearTags() noexcept {
        for(T &elt : *this) { elt.tag = 0; }
    }

    void Tag(Handle h, int tag) {
        auto it = FindByIdNoOops(h);
        if (it != nullptr) {
            it->tag = tag;
        }
    }

    void RemoveTagged() {
        size_t transfer_idx = 0;
        for(size_t i = 0; i < Size(); ++i) {
            // Manually shifting the items by bubbling them forwards into the free list
            // area means taht we're only moving each target index once, giving us a nice
            // O(n) complexity, instead of the naive O(n^2) solution, which is to repeatedly
            // perform removal by erasing each target index, like `RemoveById()` does.
            //
            // NOTE: this is only relevant if multiple elements are tagged. Otherwise, if
            // only a single element is tagged, due to the bubbling behaviour this way is
            // slower than simply calling `RemoveById()` on it.
            const auto target = targets[i];
            Storage &at_target = elemstore[target.value];
            if(at_target.get()->tag != 0) {
                at_target.reset();
            } else {
                if(transfer_idx < i) {
                    const auto transfer_target = targets[transfer_idx];
                    targets[transfer_idx] = target;
                    targets[i] = transfer_target;
                }
                ++transfer_idx;
            }
        }
        used = transfer_idx;
    }

    void RemoveById(Handle h) {
        const size_t idx = FindInsertionPoint(h);
        if(idx < Size()) {
            const Target target = targets[idx];
            Storage &at_target = elemstore[target.value];
            if(at_target.get()->h.v == h.v) {
                at_target.reset();
                --used;
                if(idx < used) {
                    targets.erase(targets.begin() + idx);
                    targets.push_back(target);
                }
            }
        }
    }

    void Clear() {
        elemstore.clear();
        targets.clear();
        used = 0;
    }

private:
    T *GetTarget(Target target) noexcept {
        return elemstore[target.value].get();
    }

    const T *GetTarget(Target target) const noexcept {
        return elemstore[target.value].get();
    }

    size_t FindInsertionPoint(Handle h) const {
        auto it = std::lower_bound(targets.begin(), targets.begin() + Size(), h,
                                   [this](Target target, Handle h) {
                                       return GetTarget(target)->h.v < h.v;
                                   });
        return it - targets.begin();
    }

    void InsertAt(size_t idx, const T &t) {
        const Target target = AllocForInsert(idx);
        elemstore[target.value].reset(t);
    }

    void InsertAt(size_t idx, T &&t) {
        const Target target = AllocForInsert(idx);
        elemstore[target.value].reset(std::forward<T>(t));
    }

    Target AllocForInsert(size_t idx) {
        if(used >= elemstore.size()) {
            targets.push_back({uint32_t(elemstore.size())});
            elemstore.emplace_back();
        }
        const Target target = targets.back();
        if(idx < used) {
            targets.pop_back();
            targets.insert(targets.begin() + idx, target);
        } else if(used < targets.size() - 1) {
            std::swap(targets[used], targets.back());
        }
        ++used;
        return target;
    }
};

class BandedMatrix {
public:
    enum {
        MAX_UNKNOWNS   = 16,
        RIGHT_OF_DIAG  = 1,
        LEFT_OF_DIAG   = 2
    };

    double A[MAX_UNKNOWNS][MAX_UNKNOWNS];
    double B[MAX_UNKNOWNS];
    double X[MAX_UNKNOWNS];
    int n;

    void Solve();
};

#define RGBi(r, g, b) RgbaColor::From((r), (g), (b))
#define RGBf(r, g, b) RgbaColor::FromFloat((float)(r), (float)(g), (float)(b))

// Note: sizeof(class RgbaColor) should be exactly 4
//
class RgbaColor {
public:
    uint8_t red, green, blue, alpha;

    float redF()   const { return (float)red   / 255.0f; }
    float greenF() const { return (float)green / 255.0f; }
    float blueF()  const { return (float)blue  / 255.0f; }
    float alphaF() const { return (float)alpha / 255.0f; }

    bool IsEmpty() const { return alpha == 0; }

    bool Equals(RgbaColor c) const {
        return
            c.red   == red   &&
            c.green == green &&
            c.blue  == blue  &&
            c.alpha == alpha;
    }

    RgbaColor WithAlpha(uint8_t newAlpha) const {
        RgbaColor color = *this;
        color.alpha = newAlpha;
        return color;
    }

    uint32_t ToPackedIntBGRA() const {
        return
            blue |
            (uint32_t)(green << 8) |
            (uint32_t)(red << 16) |
            (uint32_t)((255 - alpha) << 24);
    }

    uint32_t ToPackedInt() const {
        return
            red |
            (uint32_t)(green << 8) |
            (uint32_t)(blue << 16) |
            (uint32_t)((255 - alpha) << 24);
    }

    uint32_t ToARGB32() const {
        return
            blue |
            (uint32_t)(green << 8) |
            (uint32_t)(red << 16) |
            (uint32_t)(alpha << 24);
    }

    static RgbaColor From(int r, int g, int b, int a = 255) {
        RgbaColor c;
        c.red   = (uint8_t)r;
        c.green = (uint8_t)g;
        c.blue  = (uint8_t)b;
        c.alpha = (uint8_t)a;
        return c;
    }

    static RgbaColor FromFloat(float r, float g, float b, float a = 1.0) {
        return From(
            (int)(255.1f * r),
            (int)(255.1f * g),
            (int)(255.1f * b),
            (int)(255.1f * a));
    }

    static RgbaColor FromPackedInt(uint32_t rgba) {
        return From(
            (int)((rgba)       & 0xff),
            (int)((rgba >> 8)  & 0xff),
            (int)((rgba >> 16) & 0xff),
            (int)(255 - ((rgba >> 24) & 0xff)));
    }

    static RgbaColor FromPackedIntBGRA(uint32_t bgra) {
        return From(
            (int)((bgra >> 16) & 0xff),
            (int)((bgra >> 8)  & 0xff),
            (int)((bgra)       & 0xff),
            (int)(255 - ((bgra >> 24) & 0xff)));
    }
};

struct RgbaColorCompare {
    bool operator()(RgbaColor a, RgbaColor b) const {
        return a.ToARGB32() < b.ToARGB32();
    }
};

class BBox {
public:
    Vector minp;
    Vector maxp;

    static BBox From(const Vector &p0, const Vector &p1);

    Vector GetOrigin() const;
    Vector GetExtents() const;

    void Include(const Vector &v, double r = 0.0);
    bool Overlaps(const BBox &b1) const;
    bool Contains(const Point2d &p, double r = 0.0) const;
};

#endif
