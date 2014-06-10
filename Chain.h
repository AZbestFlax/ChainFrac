#ifndef CHAIN_H_INCLUDED
#define CHAIN_H_INCLUDED

#include <math.h>
#include <iostream>
#include <vector>

#define EPS 1e-8
/*

        Frac like

        A = k*(n+sqrt(f))/d

*/

using namespace std;

struct frac {

private:

    int f;
    int n;
    int d;
    int k;

    int gcd(int a, int b);
    int abs(int a);
    bool full_square();

public:
    frac();
    frac(int _f);
    frac(int _n, int _d);
    frac(int _f, int _n, int _d);
    frac(int _f, int _n, int _d, int _k);
    frac(const frac &);                     //конструктор копіювання
    ~frac();

    void set_param(int _f);
    void set_param(int _n, int _d);
    void set_param(int _f, int _n, int _d);
    void set_param(int _f, int _n, int _d, int _k);

    int get_k() {
        return k;
    }

    int get_n() {
        return n;
    }

    int get_f() {
        return f;
    }

    int get_d() {
        return d;
    }

    friend bool operator==(const frac& left, const frac& right);
    friend bool operator == (const frac& left, const int & value );
    friend bool operator != (const frac& left, const int & value );

    //friend frac& operator+=(frac& left, const frac& right);
    friend const frac operator+(const frac& left, const frac& right);
    friend const frac operator-(const frac& left, const frac& right);
    friend const frac operator-(const frac& left, const int& right);
    friend const frac operator+(const frac& left, const int& right);

    int integer();
    void simplify();
    void print();
    frac converse();

    frac& operator=(const frac& right) {
        //перевірка на самоприсвоювання
        if (this == &right) {
            return *this;
        }
        f = right.f;    n = right.n;    d = right.d;    k = right.k;
        return *this;
    }
};

frac::frac() {
    f = 0;  n = 0;  d = 1;  k = 1;
}

frac::frac(int _f) {
    f = _f; n = 0;  d = 1;  k = 1;
}

frac::frac(int _n, int _d) {
    f = 0;  n = _n; d = _d; k = 1;
}

frac::frac(int _f, int _n, int _d) {
    f = _f; n = _n; d = _d; k = 1;
}

frac::frac(int _f, int _n, int _d, int _k) {
    f = _f; n = _n; d = _d; k = _k;
}

frac::frac(const frac &otherFrac) {
    f = otherFrac.f;    n = otherFrac.n;
    d = otherFrac.d;    k = otherFrac.k;
}

frac::~frac() {
}

void frac::set_param(int _f) {
    f = _f; n = 0;  d = 1;  k = 1;
}

void frac::set_param(int _n, int _d) {
    f = 0;  n = _n; d = _d; k = 1;
}

void frac::set_param(int _f, int _n, int _d) {
    f = _f; n = _n; d = _d; k = 1;
}

void frac::set_param(int _f, int _n, int _d, int _k) {
    f = _f; n = _n; d = _d; k = _k;
}

bool operator==(const frac& left, const frac& right) {
    return (left.f == right.f)&&
           (left.n == right.n)&&
           (left.d == right.d)&&
           (left.k == right.k);
}

bool operator == (const frac& left, const int & value ) {
    return (left.f==0&&left.k==1&&left.d==1&&left.n==value);
}

bool operator != (const frac& left, const int & value ) {
    return !(left==value);
}
/*
const frac operator+(const frac& left, const frac& right) {
    //повернути конструктором
}
*/
const frac operator-(const frac& left, const int& right) {
    int _f = left.f*left.k*left.k;
    int _n = left.k*left.n - right*left.d;
    if (left.k<0) _n*=-1;
    int kk;
    (left.k<0)?kk=-1:kk=1;
    frac t(_f,_n,left.d,kk);
    t.simplify();
    return t;
};

const frac operator+(const frac& left, const int& right) {
    int _f = left.f*left.k*left.k;
    int _n = left.k*left.n + right*left.d;
    if (left.k<0) _n*=-1;
    int kk;
    (left.k<0)?kk=-1:kk=1;
    frac t(_f,_n,left.d,kk);
    t.simplify();
    return t;
};

int frac::integer() {
    return k*(sqrt(f)+n)/d;
}

int frac::gcd(int a, int b) {
    if (b==0) {return a;}
    else {gcd(b, a%b);}
}

int frac::abs(int a) {
    return (a<0)?-a:a;
}

bool frac::full_square() {
    return fabs((int)(sqrt(f))-sqrt(f))<EPS;
}

frac frac::converse() {
    int new_d = k*(f-n*n);
    frac* t;
    t = new frac(f,-n,new_d,d);
    t->simplify();
    return *t;
};

void frac::simplify() {
    int u,_n;
    if (d<0) {d*=-1; k*=-1;}
//========================================== if simple frac
        if (f==0) {
            _n = n*k;
            u = gcd(abs(_n),abs(d));
            n=_n/u;
            d=d/u;
            k=1;
        } else {
//========================================== radiacal is full square
            if (this->full_square()) {
                _n=k*((int)(sqrt(f))+n);
                f=0;
                u = gcd(_n,d);
                n=_n/u;
                d=d/u;
                k=1;
            }
//========================================== extract from radical
            u = sqrt(f);
            int result;
            result = 1;
            for (int i=2; i<=u; ++i) {
                if (f%(i*i)==0) {
                    result *= i;
                    f/=(i*i);
                }
            };
    //========================================= extract from brackets
            u = gcd(result, abs(n));
            k*=u;
            result/=u;
            n/=u;
            f*=(result*result);
    //========================================= reduce
            u = gcd(abs(k),abs(d));
            k/=u;
            d/=u;
        }
}

void frac::print() {
    simplify();
    if (*this==0) {cout<<"0\n"; return;}
    if (f==0) {
        if (d==1) cout << n << endl; else
        cout << n << "/" << d << endl;
    } else {
        if (k==1) {
            if (d==1) {
                if (n) cout << n << "+V" << f << endl;
                else cout << "V" << f << endl;
            } else cout << "(" << n << "+V" << f << ")/" << d << endl;
        } else {
            if (k==-1) {
                if (d==1) {
                    if (n>0) {cout << "-" << n;} else {cout << abs(n);}
                    cout << "-V" << f << endl;
                } else {
                    cout << "(";
                    if (n>0) {cout << "-" << n;} else {cout << abs(n);}
                    cout << "-V" << f << ")/" << d << endl;
                }
            } else {
                cout << k << "*(" << n << "+V" << f << ")/" << d << endl;
            }
        }
    }
}
//======================================================================================================================================================
//================================================ P O L Y N O M =======================================================================================
//======================================================================================================================================================

int gcd(int a, int b) {
    return (b)?gcd(b,a%b):a;
}

int abs(int a) {
    return (a<0)?-a:a;
}

class poly {
private:
    void del_lead_zero();
public:
    vector <int> p;
    poly();
    poly(int a);                // {a=0}
    poly(int a, int b);         // {a*x+b=0}
    poly(int a, int b, int c);  // {a*x^2+b*x+c=0}
    poly(vector <int> &v);
//    poly(poly& x);

    void print(void);
    void simplify();
    int deg(void);
    int get_w_by_number(int n);

    poly& operator=(const poly& right) {
        //перевірка на самоприсвоювання
        if (this == &right) {
            return *this;
        }
        p = right.p;
        return *this;
    }

    friend bool operator==(const poly& left, const poly& right);
    friend bool operator != (const poly& left, const poly& value );
    friend const poly operator-(const poly& x);
    friend poly& operator+=(poly& left, const poly& right);
    friend poly& operator*=(poly& left, const poly& right);
    friend poly& operator*=(poly& left, const int& right);
    friend const poly operator+(const poly& left, const poly& right);
    friend const poly operator-(const poly& left, const poly& right);
    friend const poly operator*(const poly& left, const poly& right);
    friend const poly operator*(const poly& left, const int& right);

};

poly::poly() {
    p.push_back(0);
}

poly::poly(int a) {
    p.push_back(a);
}

poly::poly(int a, int b) {
    p.push_back(b);
    p.push_back(a);
}

poly::poly(int a, int b, int c) {
    p.push_back(c);
    p.push_back(b);
    p.push_back(a);
}

poly::poly(vector <int> &v) {
    int i;
    for (i=0;i<v.size();++i) {
        p.push_back(v[i]);
    }
}

void poly::del_lead_zero() {
    int i=this->p.size()-1;
    while(!this->p[i] && i>0) --i;
    ++i;
    this->p.resize(i);
}

void poly::simplify() {
    if (this->p.size()<2) return;
    int i=this->p.size()-1, g=abs(this->p[i]);
    --i;
    for(;i>=0;--i) g=gcd(g,abs(p[i]));
    for(i=0;i<this->p.size();++i) this->p[i]/=g;
}

int poly::deg() {
    return p.size();
}

void poly::print() {
    int i=p.size()-1;
    if (i>1) {
        cout<<p[i]<<"w^"<<i;
        --i;
    }
    for (;i>1;--i) {
        if (p[i]>0) cout <<"+";
        if (p[i]) cout<<p[i]<<"w^"<<i;
    }
    if (p[1]>0 && p.size()>2) cout << "+";
    if (p[1]) cout<<p[1]<<"w";
    if (p[0]>0 && p.size()>1) cout << "+";
    if (p[0]) cout<<p[0];
}

int poly::get_w_by_number(int n) {
    return p[n];
}

bool operator==(const poly& left, const poly& right) {
    return (left.p == right.p);
}

bool operator != (const poly& left, const poly& right) {
    return !(left.p==right.p);
}

const poly operator-(const poly& x) {
    vector <int> v;
    int i=0;
    for (;i<x.p.size();++i) {
        v.push_back(-x.p[i]);
    }
    return poly(v);
}

const poly operator+(const poly& left, const poly& right) {
    poly t(left.p[0]+right.p[0]);
    int i=1;
    int mn=min(left.p.size(),right.p.size());
    for(;i<mn;++i) {
        t.p.push_back(left.p[i]+right.p[i]);
    }
    if (left.p.size()>right.p.size()) {
        for(;i<left.p.size();++i) t.p.push_back(left.p[i]);
    } else {
        for(;i<right.p.size();++i) t.p.push_back(right.p[i]);
    }
    t.del_lead_zero();
    return t;
};

const poly operator-(const poly& left, const poly& right) {
    poly t=-right;
    t=left+t;
    return t;
}

const poly operator*(const poly& left, const poly& right) {
    poly t(0);
    t.p.resize(left.p.size()+right.p.size());
    int i=0,j;
    for(;i<left.p.size();++i) {
        j=0;
        for(;j<right.p.size();++j) {
            t.p[i+j]+=left.p[i]*right.p[j];
        }
    }
    t.del_lead_zero();
    return t;
};

const poly operator*(const poly& left, const int& right) {
    poly t(0);
    int i=0;
    for(;i<left.p.size();++i) t.p[i]=left.p[i]*right;
    t.del_lead_zero();
    return t;
};

poly& operator+=(poly& left, const poly& right) {
    int i=0;
    vector <int> v;
    int mn=min(left.p.size(),right.p.size());
    for(;i<mn;++i) {
        v.push_back(left.p[i]+right.p[i]);
    }
    if (right.p.size()>left.p.size()) {
        for(;i<right.p.size();++i) v.push_back(right.p[i]);
    } else {
        for(;i<left.p.size();++i) v.push_back(left.p[i]);
    }
    left.p=v;
    left.del_lead_zero();
    return left;
}

poly& operator*=(poly& left, const poly& right) {
    vector <int> v;
    v.resize(left.p.size()+right.p.size());
    int i=0,j;
    for(;i<left.p.size();++i) {
        j=0;
        for(;j<right.p.size();++j) {
            v[i+j]+=left.p[i]*right.p[j];
        }
    }
    left.p=v;
    left.del_lead_zero();
    return left;
}

poly& operator*=(poly& left, const int& right) {
    int i=0;
    for(;i<left.p.size();++i) left.p[i]*=right;
    left.del_lead_zero();
    return left;
}

//======================================================================================================================================================
//========================================= P O L Y N O M - F R A C ====================================================================================
//======================================================================================================================================================
class polyfrac {
private:
    poly n,d;
public:
    polyfrac();
    polyfrac(poly &_n, poly &_d);
    friend polyfrac& operator+=(polyfrac& left, const int& right);
    polyfrac& operator=(const polyfrac& right) {
        //перевірка на самоприсвоювання
        if (this == &right) {
            return *this;
        }
        n = right.n;
        d = right.d;
        return *this;
    }

    void print(void);
    void converse(void);
    poly to_solve(void); //формування многочлена для розв*язку
    frac solve(bool sign);        //пошук розв*язку
};

polyfrac::polyfrac() {
    n=poly(0);
    d=poly(1);
}

polyfrac::polyfrac(poly &_n, poly &_d) {
    n=_n;
    d=_d;
}

void polyfrac::print(void) {
    cout<<"(";
    n.print();
    cout << ")/(";
    d.print();
    cout<<")";
}

polyfrac& operator+=(polyfrac& left, const int& right) {
    poly t=left.d;
    t*=right;
    left.n+=t;
    return left;
}

void polyfrac::converse(void) {
    poly t=this->n;
    this->n = this->d;
    this->d=t;
}

poly polyfrac::to_solve() {
    poly t=this->d;
    t*=poly(1,0);
    t=t-this->n;
    t.simplify();
    return t;
}

frac polyfrac::solve(bool sign) {
    poly t=this->to_solve();
    if (t.p.size()==3) {
        int f=t.p[1]*t.p[1]-4*t.p[0]*t.p[2];
        int n=-t.p[1];
        int d=2*t.p[2];
        frac x;
        sign=(sign)==(t.p[2]>0);
        (sign)?x.set_param(f,n,d):x.set_param(f,-n,-d);
        x.simplify();
        return x;
    }
    return frac(0);
}
//======================================================================================================================================================
//================================================ C H A I N ===========================================================================================
//======================================================================================================================================================
class chain {
private:
    frac r;
    vector <int> v;
    bool period;
    int startPeriod, endPeriod;

    int find(const frac& p, const vector <frac>& _w) {
        int i;
        for (i=0; i<_w.size(); ++i) {
            if (p==_w[i]) return i;
        }
        return i;
    }

public:
    chain();
    chain(const frac _f);
    chain(vector <int> _v);
    chain(vector <int> _v, const int& _b, const int& _e);
    ~chain();

    void convertToChain();
    bool convertToFrac();
    void printChain();
    void printFrac();
};

chain::chain() {
    frac _f(0,0,1,1);
    period = false;
    startPeriod = 0;
    endPeriod = 0;
    r = _f;
}

chain::chain(const frac _f) {
    r = _f;
    period = false;
    startPeriod = 0;
    endPeriod = 0;
}

chain::chain(vector <int> _v) {
    v = _v;
    period = false;
    startPeriod = 0;
    endPeriod = 0;
}

chain::chain(vector <int> _v, const int& _b, const int& _e) {
    this->v = _v;
    period = true;
    startPeriod = _b;
    endPeriod = _e;
}

chain::~chain() {
}

void chain::convertToChain() {
    vector <frac> w;
    frac x(r);
    long int q, b;
    do {
        w.push_back(x);
        q = x.integer();
        v.push_back(q);
        x = x-q;
        if (x!=0) x = x.converse();
        b = find(x, w);
    } while (b==w.size()&&(x!=0));
    period = b!=w.size();
    if (period) {
        startPeriod = b;
        endPeriod = w.size()-1;
    }
}

void chain::printChain() {
    cout << "[";
    for (int i=0; i<v.size(); ++i) {
        if (period&&(i==startPeriod)) cout << "(";
        cout << v[i];
        if (i==0) {cout << ";";} else {cout << ",";}
    }
    cout << "\b";
    if (period) cout << ")";
    cout << "]" << endl;
}

bool chain::convertToFrac() {
    int e=v.size();
    if (!period&&(e>0)) {
        frac x(v[e-1],1);
        for (int i=e-2; i>=0; --i) {
            x = x.converse();
            x = x + v[i];
        }
        r = x;
    } else { //є період, дріб містить радикал
        poly x(1);
        poly y(1,0);
        polyfrac a(x,y);
        e=this->endPeriod;
        for(;e>startPeriod;--e) {
            a+=v[e];
            a.converse();
        }
        a+=this->v[e];
        --e;
        frac fx=a.solve(v[e]>0);
        for(;e>=0;--e) {
            fx=fx.converse();
            fx=fx+v[e];
        }
        this->r=fx;
    }
    return !period;
};

void chain::printFrac() {
    r.print();
}

#endif // CHAIN_H_INCLUDED
