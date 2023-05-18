
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <string>
#include <unordered_map>
#include <random>
#include <ctime>
#include <vector>
#include <scip/scipdefplugins.h>
#include <scip/scip.h>
#define a 4          //每个vm的大小
#define s_mount 4    //物理机s的数量
#define r_mount 1000//请求r的数量
#define Cs 24       //每个物理机s的大小
#define t_ 10
#define T 3000        //系统容忍的最大延迟时间
#define slight_dynamic 0
#define magnitude_dynamic 1
using namespace std;
typedef struct node {
    int j;
    int k;
    float p;
};
typedef struct cost_ {
    float data;
    cost_ * next;
}cost_;
typedef struct time_ {
    float data;
    time_* next;
}time_;

//生产随机小数
bool cmp(node x, node y) {
    return(x.p > y.p);
}
float getrandom(float min, float max, default_random_engine::result_type seed = 0) {
    static default_random_engine e(seed);
    static default_random_engine::result_type last_seed = seed;
    if (seed != last_seed) {
        e = default_random_engine(seed);
        last_seed = seed;
    }
    uniform_real_distribution<float> u(min, max);
    return(u(e));
}
//生成随机整数
int getrandomint(int min, int max, default_random_engine::result_type seed = 0) {
    static default_random_engine e(seed);
    static default_random_engine::result_type last_seed = seed;
    if (seed != last_seed) {
        e = default_random_engine(seed);
        last_seed = seed;
    }
    uniform_int_distribution<int> u(min, max);
    return(u(e));
}

//物理机
class S {
friend class problemsovler;
private:
    int C;     //物理机容量
    int n;     //最大容纳VM数量
    vector<bool> x; //x数组
    vector<float> l;//负载
    vector<int> o_vm;
public:
    S() :C(Cs){
        n = C / a;
        x.resize(n);
        l.resize(n);
    };
};
//请求类
class RR {
    friend class problemsovler;
    friend bool compare1(RR r1, RR r2);
private:
    int end;
    cost_* cost_n;
    cost_* f_r;//cost 链表
    float t_r;//更新所需时间
    int bs;   //更新前上一个分配的物理机号
    int bn;   //更新前上一个分配的物理机上的vm号
public:
    RR() {
        f_r = new cost_;
        f_r->next = nullptr;
        t_r = t_;
        cost_n = f_r;
    }
    RR(int start, int end, float cost) {
        this->end = end;
        f_r = new cost_;
        f_r->data = cost;
        f_r->next = nullptr;
        t_r = t_;
        cost_n = f_r;
    }
    int update() {
        if (cost_n->next != nullptr) {//更新
            
            cost_n = cost_n->next;
            return(0);
        }
        else return(1);//结束
    }
};
bool compare1(RR r1, RR r2) {
    return (r1.cost_n->data < r2.cost_n->data);
}
class problemsovler {
private:
    ifstream inFile;
    vector<RR> r_all; //all请求集合
    vector<RR> r; //正在处理的请求集合
    unordered_map<string, int> map;
    int num; //迭代次数
    vector<S> s; //物理机集合
    int m_vm = Cs / a; //vm数量
    int m_s = s_mount;//s数量
    int m_r = r_mount;//r数量
    int r_m = 0;
public:
    problemsovler() {
        inFile = ifstream("data.csv", ios::in);
        if (!inFile) {
            cout << "打开文件失败！" << endl;
            exit(1);
        }
        string line;
        string word;
        istringstream sin;//将整行字符串line读入到字符串istringstream中
        while (getline(inFile, line)) {
            vector<int> temp(4);
            int i = 0;
            float fl;
            sin.clear();
            sin.str(line);
            for (; i < 4; i++) {
                getline(sin, word, ',');
                temp[i] = atoi(word.data());
            }
            for (; i < 11; i++)
                getline(sin, word, ',');
            fl = stof(word);

            if (map.find(to_string(temp[2]) + "+" + to_string(temp[3])) == map.end()) {//该task第一次出现
                r_all.emplace_back(temp[0] / 1000000 - 600, temp[1] / 1000000 - 600, fl);
                map[to_string(temp[2]) + "+" + to_string(temp[3])] = r_all.size() - 1;
            }

            else {
                int index = map[to_string(temp[2]) + "+" + to_string(temp[3])];

                cost_* j = r_all[index].f_r;
                while (j->next != nullptr) { j = j->next; }

                j->next = new cost_;
                j->next->data = fl;//更新cost列表
                j->next->next = nullptr;
            }
        }
        num = 0;
        s.resize(m_s);//初始化s数组
    }
    void show() {
        cout << "第" << num << "次迭代：" << endl;
        cout << "x[][]:" << endl;
        int amount = 0;
        for (int j = 0; j < m_s; j++) {
            for (int k = 0; k < m_vm; k++) {
                cout << s[j].x[k] << " , ";
                if (s[j].x[k] == 1) amount++;
            }
            cout << endl;
        }
        cout << "vm负载：";
        for (int j = 0; j < m_s; j++) {
            for (int k = 0; k < m_vm; k++) {
                cout << s[j].l[k] << " , ";
            }
            cout << endl;
        }
        cout << "badput rates: "<<endl;
        float badput = 0;
        for (int j = 0; j < m_s; j++) {
            for (int k = 0; k < m_vm; k++) {
                if ((s[j].l[k] >= a)&&(s[j].x[k]==1)) {
                    badput+=(s[j].l[k]-a)/a;
                } 
            }
        }
        cout << badput / amount << endl;
        cout << "请求总数为" << m_r << endl;
        cout << "x求和为" << amount << endl;
     
    }

    //读入初始数量的r并为所有r随机分配vm
    void random_assign(int mode) {
        if (mode == slight_dynamic) {
            for (; r_m < 0.8 * r_mount; r_m++) {
                r.push_back(r_all[r_m]);
            }
        }
        else {
            for (; r_m < 0.5 * r_mount; r_m++) {
                r.push_back(r_all[r_m]);
            }
        }
        sort(r.begin(), r.end(),compare1);
        //first fit算法
        for (auto i = r.begin(); i != r.end(); i++) {
            i->bs = getrandomint(0, s_mount - 1, time(0));
            for (int j = 0; j < Cs / a; j++) {
                if (s[i->bs].l[j] + i->cost_n->data < a) {
                    s[i->bs].o_vm.push_back(j);
                    s[i->bs].l[j] += i->cost_n->data;
                    i->bn = j;
                    break;
                }
            }
        }
        m_r = r.size();
    }
    //更新请求的状态
    void update_all(int mode) {
        //更新所有请求的状态
        

        for (auto i = r.begin(); i != r.end();) {
            auto tem = i->cost_n->data;
            if (i->update()) {//返回1说明该删除了
                s[i->bs].l[i->bn] -= i->cost_n->data;
                i = r.erase(i);
            }
            else {//已经update过了，直接i++
                s[i->bs].l[i->bn] = s[i->bs].l[i->bn] - tem + i->cost_n->data;
                i++;
            }
        }
        //新的请求加入
        if (mode == slight_dynamic) {
            for (int i = 0; i < 0.04 * r_mount; i++) {
                r_m++;
                r.push_back(r_all[r_m]);
                RR& x = r.back();
                x.bs = getrandomint(0, s_mount - 1, time(0));
                while (s[x.bs].o_vm.size() == 0) {
                    x.bs = (x.bs + 1) % s_mount;
                }
                x.bn = s[x.bs].o_vm[getrandomint(0,s[x.bs].o_vm.size()-1, time(0))];
                s[x.bs].l[x.bn] += x.cost_n->data;
            }
        }
        else {
            for (int i = 0; i < 0.1 * r_mount; i++) {
                r_m++;
                r.push_back(r_all[r_m]);
                RR& x = r.back();
                x.bs = getrandomint(0, s_mount - 1, time(0));
                while (s[x.bs].o_vm.size() == 0) {
                    x.bs = (x.bs + 1) % s_mount;
                }
                x.bn = s[x.bs].o_vm[getrandomint(0, s[x.bs].o_vm.size() - 1, time(0))];
                s[x.bs].l[x.bn] += x.cost_n->data;
            }
        }
        m_r = r.size();
        num++; //迭代次数++
    }
    void solve() { //主要的算法函数
        vector<vector<vector<bool>>> b;
        //初始化b数组
        b.resize(m_r);
        for (int i = 0; i < m_r; i++) {
            b[i].resize(m_s);
            for (int j = 0; j < m_s; j++) {
                b[i][j].resize(m_vm);
            }
        }
        //b数组赋值
        for (int i = 0; i < m_r; i++) {
            b[i][r[i].bs][r[i].bn] = 1;
        }
        /*利用线性规划求解器求解:
          最优化目标： x数组求和
          输入：b[][][],输出:y[][][]
          约束：1.对于每个i，有且仅有有一对j，k使得y[i][j][k]=1
                2.对于每个j，k，对r[i].f_r求和小于a
                3.y[i][j][k]（1-b[i][j][k]）*t_r对i，j，k求和小于T/t
                4.y[i][j][k]<=x[j][k]    */
        auto b0 = b;
        int change_nums = 0;
        problem_solver2(b);
        //problem_solver3(b,b0);
         //利用求解器的结果得到x与y  y的结果赋值到r数组中的ys与yn
        for (int j = 0; j < b[0].size(); j++) {
            for (int k = 0; k < b[0][0].size(); k++) {
                s[j].x[k] = 0;
                s[j].l[k] = 0;
                
            }
            s[j].o_vm.clear();
        }
        for (int i = 0; i < b.size(); i++) {
            for (int j = 0; j < b[0].size(); j++) {
                for (int k = 0; k < b[0][0].size(); k++) {
                    if (b[i][j][k] == 1) {
                        if (b0[i][j][k] != 1) change_nums++;
                        r[i].bs = j;
                        r[i].bn = k;
                        s[j].x[k] = 1;
                        s[j].l[k] += r[i].cost_n->data;
                    }
                }
            }
        }
        cout << "迁移时间:"<<change_nums*10<<"ms"<<endl;
        for (int j = 0; j < b[0].size(); j++) {
            for (int k = 0; k < b[0][0].size(); k++) {
                if (s[j].x[k] == 1) {
                    s[j].o_vm.push_back(k);
                }
            }
        }
    }
    //对于过载的vm，开启一个新的vm，转一半的请求过去
    int problem_solver2(vector<vector<vector<bool>>>& b) {
        int count = 0;
        for (int i = 0; i < s_mount; i++) {
            for (int j = 0; j < b[0][0].size(); j++) {
                if (s[i].l[j] > a) {
                    int newi = getrandomint(0, s_mount - 1, time(0));
                    int newj = getrandomint(0, b[0][0].size() - 1, time(0));
                    while (s[newi].x[newj] == 1) {
                        newi = (newi + 2) % b[0].size();
                        newj = (newj + 1) % b[0][0].size();
                    }
                    for (int m = 0; m < b.size(); m++) {
                        if (b[m][i][j] == 1 && getrandom(0,1,time(0))>0.5 &&count <T/t_) {
                            count++;
                            b[m][i][j] = 0;
                            b[m][newi][newj] = 1;
                        }
                    }
                }
            }
            
        }
        return 0;
    }

    int problem_solver3(vector<vector<vector<bool>>>& b, vector<vector<vector<bool>>>& b0) {
        vector<int> temp;
        for (int i = 0; i < b.size(); i++) {
            for (int j = 0; j < b[0].size(); j++) {
                for (int k = 0; k < b[0][0].size(); k++) {
                    b[i][j][k] = 0;
                }
            }
        }
        for (int j = 0; j < b[0].size(); j++) {
            for (int k = 0; k < b[0][0].size(); k++) {
                s[j].x[k] = 0;
                s[j].l[k] = 0;

            }
            s[j].o_vm.clear();
        }
        sort(r.begin(), r.end(), compare1);
        for (auto i = r.begin(); i != r.end(); i++) {
            i->bs = getrandomint(0, s_mount - 1, time(0));
            for (int j = 0; j < Cs / a; j++) {
                if (s[i->bs].l[j] + i->cost_n->data < a) {
                    s[i->bs].o_vm.push_back(j);
                    s[i->bs].l[j] += i->cost_n->data;
                    i->bn = j;
                    b[i - r.begin()][i->bs][i->bn] = 1;
                    break;
                }
            }
        }
        return 0;
    }
    int problem_solver( vector<vector<vector<bool>>>& b) {
        SCIP* scip = nullptr;

        /* initialize SCIP environment */
        SCIP_CALL(SCIPcreate(&scip));
        /* include default plugins */
        SCIP_CALL(SCIPincludeDefaultPlugins(scip));
        SCIP_CALL(SCIPcreateProbBasic(scip, "SCIP_example"));
        //定义一个最小化问题
        SCIP_CALL(SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE));
        string name;
        //定义y变量
        vector<vector<vector< SCIP_VAR*>>> y_vars(b.size(),vector<vector< SCIP_VAR*>>(b[0].size(),vector< SCIP_VAR*>(b[0][0].size())));
        //将y[i][j][k]全部加入scip求解器中
        for (int i = 0; i < b.size(); i++) {
            for (int j = 0; j < b[0].size(); j++) {
                for (int k = 0; k < b[0][0].size(); k++) {
                    SCIP_VAR* var = nullptr;
                    name = "y[" + to_string(i) + "][" + to_string(j) + "][" + to_string(k) + "]";
                    SCIP_CALL(SCIPcreateVarBasic(
                        scip,                  // SCIP environment
                        &var,                  // reference to the variable
                        name.data(), // name of the variable
                        0.0,                   // Lower bound of the variable
                        1.0,                   // upper bound of the variable
                        0.0,                   // Obj. coefficient. Doesn't really matter for this problem
                        SCIP_VARTYPE_CONTINUOUS    // CONTINUOUS variable
                        ));
                    SCIP_CALL(SCIPaddVar(scip, var));
                    y_vars[i][j][k] = var;
                } 
            }
        }
        //定义x变量
        vector<vector< SCIP_VAR*>> x_vars(b[0].size(), vector<SCIP_VAR*>(b[0][0].size()));
        //将y[i][j][k]全部加入scip求解器中
        for (int j = 0; j < b[0].size(); j++) {
            for (int k = 0; k < b[0][0].size(); k++) {
                SCIP_VAR* var = nullptr;
                name = "x[" + to_string(j) + "][" + to_string(k) + "]";
                SCIP_CALL(SCIPcreateVarBasic(
                    scip,                  // SCIP environment
                    &var,                  // reference to the variable
                    name.data(), // name of the variable
                    0.0,                   // Lower bound of the variable
                    1.0,                   // upper bound of the variable
                    1.0,                   // Obj. coefficient. Doesn't really matter for this problem
                    SCIP_VARTYPE_CONTINUOUS    // continuous variable
                ));
                SCIP_CALL(SCIPaddVar(scip, var));
                x_vars[j][k] = var;
            }
        }
        vector<SCIP_CONS*> con1;  //第一类约束:对于每个i，有且仅有有一对j，k使得y[i][j][k]=1
        for (int i = 0; i < b.size(); i++) {
            SCIP_CONS* cons = nullptr;
            name = "r[" + to_string(i) + "]";
            SCIP_CALL(SCIPcreateConsBasicLinear(scip,
                &cons,
                name.data(),
                0,
                nullptr,
                nullptr,
                1.0,//下界
                1.0));// 上界
            for (int j = 0; j < b[0].size(); j++) {
                for (int k = 0; k < b[0][0].size(); k++) {
                    SCIP_CALL(SCIPaddCoefLinear(scip, cons, y_vars[i][j][k], 1.0));
                }
            }
            SCIP_CALL(SCIPaddCons(scip, cons));
            con1.push_back(cons);
        }
        
        vector<SCIP_CONS*> con2;//对于每个j，k，对y[i]的f_r求和小于a
        for (int j = 0; j < b[0].size(); j++) {
            for (int k = 0; k < b[0][0].size(); k++) {
                SCIP_CONS* cons = nullptr;
                name = "y[" + to_string(j) + "][" + to_string(k) + "]";
                SCIP_CALL(SCIPcreateConsBasicLinear(scip,
                &cons,
                name.data(),
                0,
                nullptr,
                nullptr,
                -SCIPinfinity(scip),//下界
                0));// 上界
                for (int i = 0; i < b.size(); i++) {
                    SCIP_CALL(SCIPaddCoefLinear(scip, cons, y_vars[i][j][k], r[i].cost_n->data));
                }
                SCIP_CALL(SCIPaddCoefLinear(scip, cons, x_vars[j][k],-a));
                SCIP_CALL(SCIPaddCons(scip, cons));
                con2.push_back(cons);
            }
        }

        vector<SCIP_CONS*> con3; //3.y[i][j][k]（1 - b[i][j][k]）对i，j，k的t_r求和小于T
        SCIP_CONS* cons = nullptr;
        name = "T_t";
        SCIP_CALL(SCIPcreateConsBasicLinear(scip,
            &cons,
            name.data(),
            0,
            nullptr,
            nullptr,
            -SCIPinfinity(scip),//下界
            T));// 上界
        for (int i = 0; i < b.size(); i++) {
            for (int j = 0; j < b[0].size(); j++) {
                for (int k = 0; k < b[0][0].size(); k++) {
                    SCIP_CALL(SCIPaddCoefLinear(scip, cons, y_vars[i][j][k], r[i].t_r*(1-b[i][j][k])));
                }
            }
        }
        SCIP_CALL(SCIPaddCons(scip, cons));
        con3.push_back(cons);
        
        vector<SCIP_CONS*> con4; //y[i][j][k] <= x[j][k]
        for (int i = 0; i < b.size(); i++) {
            for (int j = 0; j < b[0].size(); j++) {
                for (int k = 0; k < b[0][0].size(); k++) {
                    SCIP_CONS* cons = nullptr;
                    name = "y[" + to_string(i) + "][" + to_string(j) + "][" + to_string(k) + "]";
                    SCIP_CALL(SCIPcreateConsBasicLinear(scip,
                        &cons,
                        name.data(),
                        0,
                        nullptr,
                        nullptr,
                        -SCIPinfinity(scip),//下界
                        0));// 上界
                    SCIP_CALL(SCIPaddCoefLinear(scip, cons, y_vars[i][j][k], 1));
                    SCIP_CALL(SCIPaddCoefLinear(scip, cons, x_vars[j][k], -1));
                    SCIP_CALL(SCIPaddCons(scip, cons));
                    con4.push_back(cons);
                }
            }
        }

        SCIP_CALL(SCIPsolve(scip));
        SCIP_SOL* sol;
        sol = SCIPgetBestSol(scip);

        //先把b全部设为0
        for (int i = 0; i < b.size(); i++) {
            for (int j = 0; j < b[0].size(); j++) {
                for (int k = 0; k < b[0][0].size(); k++) {
                    b[i][j][k] = 0;
                }
            }
        }
       // cout << "y:" << endl;
        //根据b的小数解得到b的整数解
        for (int i = 0; i < b.size(); i++) {
            vector<node> vec;//保存概率
            for (int j = 0; j < b[0].size(); j++) {
                for (int k = 0; k < b[0][0].size(); k++) {
                    //cout << SCIPgetSolVal(scip, sol, y_vars[i][j][k]) << ",";
                    if (SCIPgetSolVal(scip, sol, y_vars[i][j][k]) !=0) {
                        node n;
                        n.j = j;
                        n.k = k;
                        n.p = SCIPgetSolVal(scip, sol, y_vars[i][j][k]);
                        vec.push_back(n);
                    }
                }
                //cout << endl;
            }
            //cout << endl;
            //vec数组的求和应该为1
            float q = getrandom(0, 1, time(0));
            sort(vec.begin(), vec.end(), cmp);
            for (auto u : vec) {
                if (u.p > q) {
                    b[i][u.j][u.k] = 1;
                    break;
                }
                else {
                    q -= u.p;
                }
            }
        }
        cout << "x[][]: " << endl;
        for (int j = 0; j < b[0].size(); j++) {
            for (int k = 0; k < b[0][0].size(); k++) {
                cout << SCIPgetSolVal(scip, sol, x_vars[j][k]) << " , ";
            }
            cout << endl;
        }
        for (int i = 0; i < b.size(); i++) {
            for (int j = 0; j < b[0].size(); j++) {
                for (int k = 0; k < b[0][0].size(); k++) {
                    SCIP_CALL(SCIPreleaseVar(scip, &y_vars[i][j][k]));
                }
            }
        }
        for (int j = 0; j < b[0].size(); j++) {
            for (int k = 0; k < b[0][0].size(); k++) {
                SCIP_CALL(SCIPreleaseVar(scip, &x_vars[j][k]));
            }
        }
        for (auto& constr : con1)
        {
            SCIP_CALL(SCIPreleaseCons(scip, &constr));
        }
        con1.clear();
        for (auto& constr : con2)
        {
            SCIP_CALL(SCIPreleaseCons(scip, &constr));
        }
        con2.clear();
        for (auto& constr : con3)
        {
            SCIP_CALL(SCIPreleaseCons(scip, &constr));
        }
        con3.clear();
        for (auto& constr : con4)
        {
            SCIP_CALL(SCIPreleaseCons(scip, &constr));
        }
        con4.clear();
        SCIP_CALL(SCIPfree(&scip));
        return(0);
    }
};
int main() {
    problemsovler P;
    P.random_assign(magnitude_dynamic);
    //P.update_all(magnitude_dynamic);
    for (int i = 0; i < 6; i++) {
        cout << "********************************************************************************************" << endl;
        P.show();
        P.solve();
        P.show();
        P.update_all(magnitude_dynamic);
    }
    return(0);
}