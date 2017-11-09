#include <iostream>
#include <string>
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
const ld EPS = 1e-10;

enum Type{
    in,
    query,
    out
};

bool isZero(ld x) {
    return -EPS <= x && x <= EPS;
}
bool areEqual(ld x, ld y) {
    return isZero(x - y);
}
bool isLessThan(ld x, ld y) {
    return x < y - EPS;
}



struct Point {
    ld x;
    ld y;
public:
    Point() : x(0), y(0) {}
    Point(ld first, ld second) : x(first), y(second){}

    ld length() const {
        return std::sqrt(x * x + y * y);
    }
    ld operator%(const Point& second) const {
        return x * second.y - y * second.x;
    }
    ld operator*(const Point& second) const {
        return x * second.x + y * second.y;
    }
    ld len2() const {
        return x*x + y*y;
    }
    ld len() const {
        return std::sqrt(len2());
    }


    bool operator ==(const Point& a) const {
        return (areEqual(x, a.x) && areEqual(y, a.y));
    }
    bool operator <(const Point& a) const {
        if (areEqual(x, a.x)) {
            return isLessThan(y, a.y);
        }
        return x < a.x;
    }
    bool operator >(const Point& a) const {
        return a < *this;
    }

    bool ycomp(const Point& a) const {
        if (areEqual(y, a.y)) {
            return isLessThan(x, a.x);
        }
        return y < a.y;
    }
    bool yequal(const Point& a) const {
        return areEqual(y, a.y);
    }

    bool xcomp(const Point& a) const {
        if (areEqual(x, a.x)) {
            return isLessThan(y, a.y);
        }
        return x < a.x;
    }
    bool xequal(const Point& a) const {
        return areEqual(x, a.x);
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
    //friend bool areSegmentrsIntersects(const Point& a, const Point& b, const Point& c, const Point& d);
};
struct Segment {
    Point start;
    Point finish;
    Segment(const Point& a, const Point& b) : start(a), finish(b) {
        if (start > finish)
            std::swap(start, finish);
    }
    bool operator <(const Segment& a) const {
        if (start.yequal(a.start)) {
            return finish.ycomp(a.finish);
        }
        return start.ycomp(a.start);
    }
    bool operator ==(const Segment& a) const {
        return start == a.start && finish == a.finish;
    }

};

struct Event{
    Point point;
    Type action;
    Segment* segment;
    u32 number;

    Event(Point a, Type b, Segment* c, u32 n = 0) : point(a), action(b), segment(c), number(n) {}
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
    return isLessThan((a.x - b.start.x) / (b.finish.x - b.start.x) * (b.finish.y - b.start.y) + b.start.y, a.y);
}
bool contains(Point a, Segment b) {
    return areEqual((a.x - b.start.x) / (b.finish.x - b.start.x) * (b.finish.y - b.start.y) + b.start.y, a.y);
}

class Node {
protected:
    //left < right
    Node* left;
    Node* right;
    Segment key;
    long long priority;
    u32 countBehind;

    u32 countBehindFunc() const {
        return this ? countBehind : 0;
    }

    void updateCount() {
        countBehind = left->countBehindFunc() + right->countBehindFunc() + 1;
    }

    void update_all() {
        updateCount();
    }


public:

    std::pair<Node*, u32> find(Point findKey) {
        if (!this)
            return std::make_pair<Node*, u32>(nullptr, 0);
        if (contains(findKey, key)) {
            return std::make_pair(this, 0);
        }

        if (isHigher(findKey, key)) {
            auto x = right->find(findKey);
            x.second += left->countBehindFunc() + 1;
            if(! x.first) {
                x.first = this;
            }
            return x;
        }
        auto x = left->find(findKey);
        if(! x.first)
            x.first = this;
        return x;
    }

    explicit Node(Segment k) : left(nullptr), right(nullptr), key(k), priority(std::rand()), countBehind(1) {}
    ~Node() {
        if (left)
            delete left;
        if (right)
            delete right;
    }

    Node* insert(Node* a) {
        if (!this)
            return a;
        if (a->priority > this->priority) {
            split(this, a->key, a->left, a->right);
            a->update_all();
            return a;
        }
        (a->key < key ? left: right) = (a->key < key ? left: right)->insert(a);
        update_all();
        return this;
    }
    Node* erase(Segment newKey) {
        long long count_left = left->countBehindFunc() + 1;
        if (newKey == key) {
            Node* x = this;
            Node* answer = merge(left, right);
            x->left = nullptr;
            x->right = nullptr;
            delete x;
            x = nullptr;
            return answer;
        }
        if (newKey < key) {
            left = left->erase(newKey);
        } else
            right = right->erase(newKey);
        update_all();
        return this;
    }

    friend void split(Node*, Segment, Node*&, Node*&);
    friend Node* merge(Node* l, Node* r);
    friend class solveTask;
};



Node* merge(Node* l, Node* r) {
    Node* my;
    if (!l || !r) {
        my = (l ? l : r);
        return my;
    }
    if (l->priority < r->priority) {
        r->left = merge(l, r->left);
        my = r;
    } else{
        l->right = merge(l->right, r);
        my = l;
    }
    my->update_all();
    return my;
}
void split(Node* myVertex, Segment inputKey, Node*& l, Node*& r) {
    if (!myVertex) {
        l = nullptr;
        r = nullptr;
        return;
    }
    if (!(inputKey < myVertex->key)) {
        split(myVertex->right, inputKey, myVertex->right, r);
        l = myVertex;
        l->update_all();
    } else {
        split(myVertex->left, inputKey, l, myVertex->left);
        r = myVertex;
        r->update_all();
    }
    myVertex->update_all();
}


class BinarySearchTree{
    Node* root;
public:
    explicit BinarySearchTree() : root(nullptr) {}
    BinarySearchTree(long long size, Segment* array) {
        root = new Node(array[0]);
        for (long long i = 1; i < size; ++i) {
            Node* ver = new Node(array[i]);
            root = root->insert(ver);
        }
    }
    ~BinarySearchTree() {
        if (root)
            delete root;
    }

    void insert(Segment key) {
        Node* newNode = new Node(key);
        if (! root) {
            root = newNode;
            return;
        }
        root = root->insert(newNode);
    }
    void erase(Segment key) {
        root = root->erase(key);
    }
    std::pair<Node*, u32> find(Point key) const {
        return root->find(key);
    }
};


class solveTask{
    u32 n, k;
    vector<Segment> Polygon;
    std::vector<Event> events;
    vector<uint8_t> answers;

    void preparation() {
        std::sort(events.begin(), events.end());
    }

    void addSegment(Point first, Point second) {
        Polygon.push_back(Segment(first, second));
        events.push_back(Event(Polygon.back().start, in, &Polygon.back()));
        events.push_back(Event(Polygon.back().finish, out, &Polygon.back()));
    }
    void printAnswer() const {
        for (int i = 0; i < k; ++i) {
            if (answers[i] == 0)
                cout << "INSIDE" << '\n';
            else if (answers[i] == 1)
                cout << "BORDER" << '\n';
            else
                cout << "OUTSIDE" << '\n';
        }
    }
    void makeAnswer() {
        BinarySearchTree tree;
        for(auto it: events) {
            if (it.action == in) {
                tree.insert(*(it.segment));
            } else if (it.action == out) {
                tree.erase(*(it.segment));
            } else {
                auto fight = tree.find(it.point);
                if (fight.first && contains(it.point, fight.first->key))
                    answers[it.number] = 1;
                else if (fight.second % 2 == 0) {
                    answers[it.number] = 2;
                } else if (fight.second % 2 == 1) {
                    answers[it.number] = 0;
                }
            }
        }
    }

public:
    void makeOneTest() {
        cin >> n;
        Polygon.reserve(n);
        ld a, b;
        cin >> a >> b;
        Point first(a, b);
        for (u32 i = 1; i < n; ++i) {
            cin >> a >> b;
            Point second(a, b);
            std::swap(first, second);
            if(first == second)
                continue;
            addSegment(second, first);
        }
        addSegment(first, Polygon[0].start);
        cin >> k;
        answers.assign(k, 0);
        events.reserve(n + k);
        for (u32 i = 0; i < k; ++i) {
            ld a, b;
            cin >> a >> b;
            Point queryPoint(a, b);
            Event queryEvent(queryPoint, query, nullptr, i);
            events.push_back(queryEvent);
        }
        preparation();
        makeAnswer();
        printAnswer();
    }

};


int main() {
    solveTask task;
    u32 t;
    cin >> t;
    for (int testNumber = 0; testNumber < t; ++testNumber) {
        task.makeOneTest();
    }
    return 0;
}