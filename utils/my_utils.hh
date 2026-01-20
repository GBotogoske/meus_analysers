#ifndef MY_UTILS_HH
#define MY_UTILS_HH

#include <vector>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>
#include <queue>
#include <stdexcept>


class QPoint
{
    public:
        QPoint()
        {
        }
        QPoint(double x0,double y0,double z0, double q0)
        {
            x=x0;
            y=y0;
            z=z0;
            q=q0;
        }
        ~QPoint() {}
        double x,y,z,q;
};


class QCluster : public std::vector<QPoint>
{
    public:
        int objID = -1;
        int type = 0; // 0 --> track, 1-->shower, 2-->PFP, 3-->Slice

        QCluster() = default;

        QCluster(const QCluster& other)
            : std::vector<QPoint>(other), objID(other.objID), type(other.type)
        {}

        QCluster& operator=(const QCluster& other)
        {
            if (this != &other)
            {
                std::vector<QPoint>::operator=(other);
                objID = other.objID;
                type = other.type;
            }
            return *this;
        }
};

class QFlash
{
public:
    QFlash() = default;
    QFlash(const QFlash& qflash) = default;
    ~QFlash() = default;

    int flashID = -1;
    std::vector<double> PE_CH;

    double y = 0.0, z = 0.0, y_err = 0.0, z_err = 0.0;
    double time = 0.0, time_err = 0.0;

    void norm_this_flash()
    {
        if (PE_CH.empty()) return;
        double sum = 0.0;
        for (size_t i = 0; i < PE_CH.size(); ++i) sum += PE_CH[i];
        if (sum <= 0.0) return;
        for (double& v : PE_CH) v /= sum;
    }


};

static inline bool isZero(double x, double eps=1e-12)
{
    return std::abs(x) <= eps;
}

struct Match 
{
    std::vector<int> row2col; // row -> col match
    std::vector<int> col2row; // col -> row match
    int size = 0;
};

static bool tryKuhn(int r,
                    const std::vector<std::vector<int>>& adj,
                    std::vector<int>& col2row,
                    std::vector<char>& vis)
{
    if (vis[r])// se ja visitou linha
    {
        return false; // casamento nao feito
    } 
    vis[r] = true; // marca que ja visitou para nao repetir o loop

    for (int c : adj[r]) // varre todas os indices de coluna com zero
    {
        if (col2row[c] == -1 || tryKuhn(col2row[c], adj, col2row, vis)) //se aquele coluna c ainda nao teve match ou caso ja esteja ocupada tenta realocaressa linha ja ocupada
        {
            col2row[c] = r; // casamento feito
            return true; //retorna verdadeiro
        }
    }
    return false; //casamento nao feito
}

static Match maxMatchingZeros(const std::vector<std::vector<double>>& M, double eps=1e-12)
{
    int p = (int)M.size();
    std::vector<std::vector<int>> adj(p); // representa em cada indice ( linha) em quais indices de coluna tem zeros
    for (int r = 0; r < p; ++r)
    {
        for (int c = 0; c < p; ++c)
        {
            if (isZero(M[r][c], eps))
            {
                adj[r].push_back(c);
            } 
        }     
    }
        
    Match mt;
    mt.row2col.assign(p, -1); // elementos de match
    mt.col2row.assign(p, -1); // elementos de match

    for (int r = 0; r < p; ++r) 
    {
        std::vector<char> vis(p, 0); // checar se ja visitou linha
        if (tryKuhn(r, adj, mt.col2row, vis)) // checar se posso casar a linha r 
        {
            mt.size++; //quantos matchs feitos
        }
    }

    for (int c = 0; c < p; ++c) 
    {
        if (mt.col2row[c] != -1)
        {
            mt.row2col[ mt.col2row[c] ] = c; // calcula a "transposta"
        }
    }

    return mt;
}

struct Cover
{
    std::vector<char> rowCov; // 1 = coberta
    std::vector<char> colCov; // 1 = coberta
    int size = 0;
};

static Cover minCoverFromMatching(const std::vector<std::vector<double>>& M,
                                 const Match& mt,
                                 double eps=1e-12)
{
    int p = (int)M.size();

    // adj de zeros --> novamente calculo em cada linha quais os indices com colunas
    std::vector<std::vector<int>> adj(p);
    for (int r = 0; r < p; ++r)
    {
        for (int c = 0; c < p; ++c)
        {
            if (isZero(M[r][c], eps)) 
            {
                adj[r].push_back(c);
            }
        }
            
    }
    
    // BFS em caminhos alternantes a partir de linhas livres
    std::vector<char> visR(p, 0), visC(p, 0);
    std::queue<int> q;

    ///varre as linhas
    for (int r = 0; r < p; ++r) 
    {
        if (mt.row2col[r] == -1) // linha não casada
        { 
            visR[r] = 1;
            q.push(r);
        }
    }

    //varre todas as linhas nao casadas
    while (!q.empty()) 
    {
        int r = q.front(); 
        q.pop();

        // linha -> col por arestas de zero NÃO casadas
        for (int c : adj[r]) 
        {
            if (mt.row2col[r] == c) 
            {
                continue; // aresta casada, pula
            }
            if (!visC[c])  // se nao chequei ainda essa coluna
            {
                visC[c] = 1;
                // col -> linha pela aresta casada (se existir)
                int r2 = mt.col2row[c];
                if (r2 != -1 && !visR[r2]) // se coluna esta casada e linha nao casada
                {
                    visR[r2] = 1;
                    q.push(r2);
                }
            }
        }
    }

    // cover mínimo = (linhas NÃO visitadas) U (colunas visitadas)
    Cover cov;
    cov.rowCov.assign(p, 0);
    cov.colCov.assign(p, 0);
    cov.size = 0;

    for (int r = 0; r < p; ++r) 
    { 
        cov.rowCov[r] = !visR[r]; 
        cov.size += cov.rowCov[r]; 
    }
    for (int c = 0; c < p; ++c) 
    { 
        cov.colCov[c] =  visC[c]; 
        cov.size += cov.colCov[c]; 
    }

    return cov;
}

static void step1_reduceRows(std::vector<std::vector<double>>& M)
{
    for (auto &row : M)
    {
        double mn = *std::min_element(row.begin(), row.end());
        for (double &x : row) 
        {
            x -= mn;
        }
    }
}

static void step2_reduceCols(std::vector<std::vector<double>>& M)
{
    int p = (int)M.size();
    std::vector<double> colMin(p, std::numeric_limits<double>::infinity());

    for (int i = 0; i < p; ++i)
    {
        for (int j = 0; j < p; ++j)
        {
            colMin[j] = std::min(colMin[j], M[i][j]);
        }      
    }
        
    for (int i = 0; i < p; ++i)
    {
        for (int j = 0; j < p; ++j)
        {
            M[i][j] -= colMin[j];
        }       
    }
        
}

static void step4_adjust(std::vector<std::vector<double>>& M,
                         const Cover& cov)
{
    int p = (int)M.size();
    double k = std::numeric_limits<double>::infinity();

    // menor elemento NÃO coberto
    for (int i = 0; i < p; ++i)
    {
        if (cov.rowCov[i]) 
        {
            continue;
        }
        for (int j = 0; j < p; ++j) 
        {
            if (cov.colCov[j]) 
            {
                continue;
            }
            k = std::min(k, M[i][j]);
        }
    }

    // ajuste
    for (int i = 0; i < p; ++i)
    {
        for (int j = 0; j < p; ++j) 
        {
            if (!cov.rowCov[i] && !cov.colCov[j])
            {
                M[i][j] -= k;  // não coberto
            }     
            else if (cov.rowCov[i] && cov.colCov[j])
            {
                M[i][j] += k;   // coberto duas vezes
            }
        }
    }
}

struct HungarianResult 
{
    std::vector<int> assign; // assign[row] = col (ou -1)
    double cost = 0.0;
};

// Minimiza custo
inline HungarianResult hungarian_min(const std::vector<std::vector<double>>& C_in, double eps=1e-12)
{
    int m = (int)C_in.size();
    int n = m ? (int)C_in[0].size() : 0;
    if (m == 0 || n == 0) return {{}, 0.0};

    // verifica retangular
    for (const auto& row : C_in)
    {
        if ((int)row.size() != n) throw std::runtime_error("Matriz nao eh retangular.");
    }

    //pega valor maximo
    int p = std::max(m, n);

    double maxC = 0.0;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            maxC = std::max(maxC, std::abs(C_in[i][j]));

    const double BIG = (maxC + 1.0) * p * 10.0; // grande e FINITO

    // STEP 0: quadratiza de mxn para pxp onde p = max(m,n)
    std::vector<std::vector<double>> M(p, std::vector<double>(p, BIG));
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            M[i][j] = C_in[i][j];
        }       
    }
        
    // STEP 1-2: reduções
    step1_reduceRows(M);
    step2_reduceCols(M);

    // Loop STEP 3/4 até matching perfeito (size == p)
    while (true) 
    {
        Match mt = maxMatchingZeros(M, eps);
        if (mt.size == p) // caso exato 1 para 1
        {
            // STEP 5: assignment vem do matching
            HungarianResult res;
            res.assign.assign(m, -1);
            res.cost = 0.0;

            for (int i = 0; i < m; ++i) 
            {
                int j = mt.row2col[i];
                if (j >= 0 && j < n) 
                {
                    res.assign[i] = j;
                    res.cost += C_in[i][j];
                } 
                else
                {
                    res.assign[i] = -1; // casou com dummy
                }
            }
            return res;
        }

        //caso tenha menos matchs que o numero de linhas
        Cover cov = minCoverFromMatching(M, mt, eps);
        // aqui, cov.size == mt.size (Kőnig). Se cov.size < p -> STEP 4
        step4_adjust(M, cov);
    }
}


#endif