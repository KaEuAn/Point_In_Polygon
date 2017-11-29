#include <iostream>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>

using std::vector;
using std::cout;
using std::cin;
using std::string;
using u32 = uint32_t;
using u64 = uint64_t;
using ld = long double;
using ll = long long;

enum ActionType{
    in,
    query,
    out
};

enum Answer{
    inside,
    border,
    outside
};


struct Point {
    int x;
    int y;
public:
    Point() : x(0), y(0) {}
    Point(int first, int second) : x(first), y(second){}
    Point(const Point& a) : x(a.x), y(a.y) {}

    ld operator%(const Point& second) const {
        return x * second.y - y * second.x;
    }
    ld operator*(const Point& second) const {
        return x * second.x + y * second.y;
    }


    bool operator ==(const Point& a) const {
        return (x == a.x && y == a.y);
    }
    bool operator <(const Point& a) const {
        if (x == a.x) {
            return y < a.y;
        }
        return x < a.x;
    }
    bool operator >(const Point& a) const {
        return a < *this;
    }

    bool operator != (const Point& a) {
        return ! (*this == a);
    }
    Point operator - (const Point& a) const {
        Point c(*this);
        c.x -= a.x;
        c.y -= a.y;
        return c;
    }
    Point& operator += (const Point& a) {
        x += a.x;
        y += a.y;
        return *this;
    }
    Point operator + (const Point& a) const {
        Point c(*this);
        return c += a;
    }
    Point operator / (double a) const {
        Point c(*this);
        c.x /= a;
        c.y /= a;
        return c;
    }
    Point operator * (double a) const {
        Point c(*this);
        c.x *= a;
        c.y *= a;
        return c;
    }
};

struct Segment {
    Point start;
    Point finish;
    Segment(const Point& a, const Point& b) : start(a), finish(b) {
        if (start > finish)
            std::swap(start, finish);
    }

    bool operator ==(const Segment& a) const {
        return start == a.start && finish == a.finish;
    }

};
struct SegmentComparator {
    bool operator() (const Segment& a, const Segment& b) {
        if (a.start.x == b.start.x) {
            if (a.start.y == b.start.y) {
                Point first = a.finish - a.start;
                Point second = b.finish - b.start;
                return first % second > 0;
            }
            return a.start.y < b.start.y;
        }
        const Segment &second = (a.start.x < b.start.x ? b : a);
        const Segment &first = (a.start.x < b.start.x ? a : b);
        ll y1 = (second.start.x - first.start.x) * (first.finish.y - first.start.y)
                + first.start.y * (first.finish.x - first.start.x);
        bool is = y1 < second.start.y * (first.finish.x - first.start.x);
        return (a.start.x == first.start.x) == is;
    }
};

struct Event{
    Point point;
    ActionType action;
    Segment* segment;
    u32 number;

    Event(Point a, ActionType b, Segment* c, u32 n = 0) : point(a), action(b), segment(c), number(n) {}
    bool operator <(const Event& a) const {
        if (point == a.point) {
            return action < a.action;
        }
        if (point.x == a.point.x)
            return point.y < a.point.y;
        return point.x < a.point.x;
    }
};

bool isHigher(Point a, Segment b) {
    if (b.start.x == b.finish.x)
        return a.y > b.finish.y;
    return (a.x - b.start.x) * (b.finish.y - b.start.y) +
                               b.start.y * (b.finish.x - b.start.x) < a.y * (b.finish.x - b.start.x);
}
bool contains(Point a, Segment b) {
    if (b.start.x == b.finish.x)
        return (a.x == b.start.x) && (a.y <= b.finish.y && b.start.y <= a.y);
    ll x = (a.x - b.start.x) * (b.finish.y - b.start.y);
    return x + b.start.y * (b.finish.x - b.start.x) == a.y * (b.finish.x - b.start.x);
}

template <typename StoredType, typename Comp = std::less<StoredType>>
class Node {
protected:
    using myNode = Node<StoredType, Comp>;
    //left < right
    Comp comp;
    myNode* left;
    myNode* right;
    StoredType key;
    long long priority;
    u32 countBehind;

    u32 countBehindFunc() const {
        return this ? countBehind : 0;
    }

    void updateCount() {
        countBehind = left->countBehindFunc() + right->countBehindFunc() + 1;
    }

    void updateAll() {
        updateCount();
    }


public:

    std::pair<myNode*, u32> find(Point findKey) {
        if (!this)
            return std::make_pair<myNode*, u32>(nullptr, 0);
        if (contains(findKey, key)) {
            return std::make_pair(this, 0);
        }

        if (! isHigher(findKey, key)) {
            auto x = left->find(findKey);
            if(! x.first)
                x.first = this;
            return x;
        }
        auto x = right->find(findKey);
        x.second += left->countBehindFunc() + 1;
        if(! x.first) {
            x.first = this;
        }
        return x;
    }

    explicit Node(Segment k) : left(nullptr), right(nullptr), key(k), priority(std::rand()), countBehind(1) {}
    ~Node() {
        if (left)
            delete left;
        if (right)
            delete right;
    }

    myNode* insert(myNode* a) {
        if (!this)
            return a;
        if (a->priority > this->priority) {
            split(this, a->key, a->left, a->right);
            a->updateAll();
            return a;
        }
        (comp(a->key, key) ? left: right) = (comp(a->key, key) ? left: right)->insert(a);
        updateAll();
        return this;
    }
    myNode* erase(Segment newKey) {
        if (newKey == key) {
            myNode* x = this;
            myNode* answer = merge(left, right);
            x->left = nullptr;
            x->right = nullptr;
            delete x;
            return answer;
        }
        if (comp(newKey, key)) {
            left = left->erase(newKey);
        } else
            right = right->erase(newKey);
        updateAll();
        return this;
    }

    template <typename CurrentType, typename Compar> friend void split(Node<CurrentType, Compar>*, CurrentType, Node<CurrentType, Compar>*&, Node<CurrentType, Compar>*&);
    template <typename CurrentType, typename Compar> friend Node<CurrentType, Compar>* merge(Node<CurrentType, Compar>* left, Node<CurrentType, Compar>* right);
    friend class SolveTaskPointsInPolygon;
};


template <typename CurrentType, typename Compar>
Node<CurrentType, Compar>* merge(Node<CurrentType, Compar>* left, Node<CurrentType, Compar>* right) {
    Node<CurrentType, Compar>* newRoot;
    if (!left || !right) {
        newRoot = (left ? left : right);
        return newRoot;
    }
    if (left->priority < right->priority) {
        right->left = merge(left, right->left);
        newRoot = right;
    } else{
        left->right = merge(left->right, right);
        newRoot = left;
    }
    newRoot->updateAll();
    return newRoot;
}
template <typename CurrentType, typename Compar>
void split(Node<CurrentType, Compar>* thisVertex, CurrentType inputKey, Node<CurrentType, Compar>*& left, Node<CurrentType, Compar>*& right) {
    Compar comp;
    if (!thisVertex) {
        left = nullptr;
        right = nullptr;
        return;
    }
    if (!comp(inputKey, thisVertex->key)) {
        split<CurrentType, Compar>(thisVertex->right, inputKey, thisVertex->right, right);
        left = thisVertex;
        left->updateAll();
    } else {
        split<CurrentType, Compar>(thisVertex->left, inputKey, left, thisVertex->left);
        right = thisVertex;
        right->updateAll();
    }
    thisVertex->updateAll();
}

template <typename StoredType, typename Comp = std::less<StoredType>>
class BinarySearchTree{
    using myNode = Node<StoredType, Comp>;
    myNode* root;
public:
    explicit BinarySearchTree() : root(nullptr) {}
    BinarySearchTree(long long size, Segment* array) {
        root = new myNode(array[0]);
        for (long long i = 1; i < size; ++i) {
            Node<StoredType, Comp>* ver = new myNode(array[i]);
            root = root->insert(ver);
        }
    }
    ~BinarySearchTree() {
        if (root)
            delete root;
    }

    void insert(Segment key) {
        myNode* newNode = new myNode(key);
        if (! root) {
            root = newNode;
            return;
        }
        root = root->insert(newNode);
    }
    void erase(Segment key) {
        root = root->erase(key);
    }
    std::pair<myNode*, u32> find(Point key) const {
        return root->find(key);
    }
};


class SolveTaskPointsInPolygon{
    vector<Segment> Polygon;
    std::vector<Event> events;
    vector<Answer> answers;

    void sortEvents() {
        std::sort(events.begin(), events.end());
    }

    void addSegment(Point first, Point second) {
        Polygon.push_back(Segment(first, second));
        events.push_back(Event(Polygon.back().start, in, &Polygon.back()));
        events.push_back(Event(Polygon.back().finish, out, &Polygon.back()));
    }
    void makeAnswer() {
        BinarySearchTree<Segment, SegmentComparator> tree;
        for(const auto& it: events) {
            if (it.action == in) {
                tree.insert(*(it.segment));
            } else if (it.action == out) {
                tree.erase(*(it.segment));
            } else {
                const auto& current = tree.find(it.point);
                if (current.first && contains(it.point, current.first->key))
                    answers[it.number] = border;
                else if (current.second % 2 == 0) {
                    answers[it.number] = outside;
                } else if (current.second % 2 == 1) {
                    answers[it.number] = inside;
                }
            }
        }
    }

public:
    vector<Answer> identifyWhereArePoints(const vector<Point> &verticies, const vector<Point> &queries) {
        u32 n = verticies.size();
        Polygon.reserve(n);
        Point first(verticies[0]);
        Point firstfirst(first);
        for (u32 i = 1; i < n; ++i) {
            Point second(verticies[i]);
            std::swap(first, second);
            if(first == second)
                continue;
            addSegment(second, first);
        }
        addSegment(first, firstfirst);
        u32 k = queries.size();
        answers.assign(k, outside);
        events.reserve(n + k);
        for (u32 i = 0; i < k; ++i) {
            Point queryPoint(queries[i]);
            Event queryEvent(queryPoint, query, nullptr, i);
            events.push_back(queryEvent);
        }
        sortEvents();
        makeAnswer();
        return answers;
    }


};

void input(vector<Point>& verticies, vector<Point>& queries) {
    u32 n;
    cin >> n;
    verticies.reserve(n);
    for (u32 i = 0; i < n; ++i) {
        ll a, b;
        cin >> a >> b;
        verticies.push_back(Point(a, b));
    }
    cin >> n;
    queries.reserve(n);
    for (u32 i = 0; i < n; ++i) {
        ll a, b;
        cin >> a >> b;
        queries.push_back(Point(a, b));
    }
}


void printAnswer(const vector<Answer>& answers) {
    for (u32 i = 0; i < answers.size(); ++i) {
        if (answers[i] == inside)
            cout << "INSIDE" << '\n';
        else if (answers[i] == border)
            cout << "BORDER" << '\n';
        else
            cout << "OUTSIDE" << '\n';
    }
}

int main() {
    u32 t;
    cin >> t;
    for (u32 testNumber = 0; testNumber < t; ++testNumber) {
        vector<Point> verticies;
        vector<Point> queries;
        input(verticies, queries);
        SolveTaskPointsInPolygon task;
        printAnswer(task.identifyWhereArePoints(verticies, queries));
        cout << '\n';
    }
    return 0;
}