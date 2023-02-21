#include <chrono>
#include <iostream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <cassert>
#include <queue>
#include <algorithm>
#include <sstream>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>


#include <thread>

#include "utils.hpp"
#include "treedecomposition.hpp"



namespace sspp {


TreeDecomposition::TreeDecomposition(size_t n) : n(n), bs(1), width(n-1), tree(2), bags({{},{}}) {
    bags[1].resize(n);
    for(size_t i = 1; i <= n; i++) {
        bags[1].push_back(i);
    }
}

TreeDecomposition::TreeDecomposition(const Graph& graph, double time) {
    bags.push_back({});
    n = graph.n();
    if (n == 0) {
        bs = 0;
        width = -1;
        tree = Graph(1);
        return;
    }
    if (n == 1) {
        bs = 1;
        tree = Graph(1);
        bags.push_back({1});
        return;
    }
    assert(n >= 2);
    auto es = graph.Edges();
    int m = es.size();
    // pipes:
    // parent to child is 0,1
    // child to parent is 2,3
    // i.e. if pipefds_1 == [0,1,2,3] then
    //     parent         child1 
    //        1 -----------> 0
    //        2 <----------- 3
    int pipefds[4];
    int ret = pipe(pipefds);
    if(ret == -1) {
        throw std::runtime_error("pipe() failed.");
    }
    ret = pipe(pipefds + 2);
    if(ret == -1) {
        throw std::runtime_error("pipe() failed.");
    }
    pid_t pid = 0;
    pid = fork();
    switch(pid) {
        case -1:
            throw std::runtime_error("fork() failed.");
        case 0:
            // child
            // duplicate
            ret = dup2(pipefds[0], STDIN_FILENO);
            if(ret == -1) {
                throw std::runtime_error("dup2() failed.");
            }
            ret = dup2(pipefds[3], STDOUT_FILENO);
            if(ret == -1) {
                throw std::runtime_error("dup2() failed.");
            }
            // close pipes
            close(pipefds[0]);
            close(pipefds[1]);
            close(pipefds[2]);
            close(pipefds[3]);
            break;
        default:
            // parent
            // close pipes
            close(pipefds[0]);
            close(pipefds[3]);
    }
    if(pid == 0) { // we are in the recursive case
        string cmd = "./flow_cutter_pace17";
        execlp(cmd.c_str(), cmd.c_str(), NULL);
    } else { 
        // we need to give input to the child
        auto start = std::chrono::system_clock::now();
        FILE *file = fdopen(pipefds[1], "w");
        fprintf(file, "p tw %d %d\n", n, m);
        for (auto e : es) {
            fprintf(file, "%d %d\n", e.F+1, e.S+1);
        }
        fflush(file);
        fclose(file); 
        cout<<"c o Primal edges "<<es.size()<<endl;
        file = fdopen(pipefds[2], "r");
        size_t read = 0;
        size_t len = 0;
        char * line = NULL;
        bool found_first = false;
        while (!found_first) {
            read = getline(&line, &len, file);
            if(read == -1) {
                throw std::runtime_error("getline() failed.");
            }
            assert(read > 0);
            found_first = strncmp("c status", line, 8) == 0;
        }
        auto passed = std::chrono::system_clock::now() - start;
        auto remaining = std::chrono::duration<double>(time) - passed;
        std::this_thread::sleep_for(remaining);
        kill(pid,SIGTERM); 
        string tmp;
        while (true) {
            read = getline(&line, &len, file);
            if(read == -1) {
                if(errno != 0) {
                    throw std::runtime_error("getline() failed.");
                } else {
                    break;
                }
            }
            std::stringstream ss(line);
            ss>>tmp;
            if (tmp == "c") continue;
            if (tmp == "s") {
                ss>>tmp;
                assert(tmp == "td");
                int nn;
                ss>>bs>>width>>nn;
                assert(nn == n);
                width--;
                bags.resize(bs + 1);
                tree = Graph(bs + 1);
            } else if (tmp == "b") {
                int bid;
                ss>>bid;
                vector<int> bag;
                int v;
                while (ss>>v) {
                    bags[bid].push_back(v-1);
                }
                std::sort(bags[bid].begin(), bags[bid].end());
            } else {
                int a = stoi(tmp);
                int b;
                ss>>b;
                tree.AddEdge(a, b);
            }
        }
        fclose(file);
        free(line);
        int status = 0;
        ret = waitpid(pid, &status, 0);
        if(ret == -1) {
            throw std::runtime_error("waitpid() failed.");
        }
        assert(status >= 0);
        assert(WIFEXITED(status));
        cout << "c o width " << Width() << endl;
        assert(Verify(graph));
    }
}

TreeDecomposition& TreeDecomposition::operator=(const TreeDecomposition& other) {
    n = other.n;
    bs = other.bs;
    width = other.width;
    tree = other.tree;
    bags = other.bags;
    return *this;
}

int TreeDecomposition::Width() const {
	return width;
}

bool TreeDecomposition::InBag(int b, int v) const {
	assert(1 <= b && b <= bs && 0 <= v && v < n);
	return BS(bags[b], v);
}

bool TreeDecomposition::dfs(int x, int v, int p, vector<int>& u) const {
	assert(InBag(x, v));
	assert(u[x] != v);
	u[x] = v;
	bool ok = true;
	for (int nx : tree.Neighbors(x)) {
		if (InBag(nx, v) && nx != p) {
			if (u[nx] == v) {
				return false;
			}
			ok = dfs(nx, v, x, u);
			if (!ok) {
				return false;
			}
		}
	}
	return true;
}

bool TreeDecomposition::Verify(const Graph& graph) const {
  assert(tree.n() == bs + 1);
  if(n != graph.n()) {
    std::cerr << n << " " << graph.n() << std::endl;
  }
	assert(n == graph.n());
	vector<vector<char>> aps(n);
	for (int i = 0; i < n; i++) {
		aps[i].resize(n);
	}
	for (auto bag : bags) {
		for (int v : bag) {
			for (int u : bag) {
				aps[v][u] = 1;
			}
		}
	}
	for (int i = 0; i < n; i++) {
		if (aps[i][i] == 0) { std::cerr << "symfail "<< i << std::endl; return false; }
	}
	for (auto e : graph.Edges()) {
		if (aps[e.F][e.S] == 0) { std::cerr << "edgefail " <<  e.F << " " << e.S << std::endl;  return false; }
	}
	vector<int> u(bs+1);
	for (int i = 1; i <= bs; i++) {
		u[i] = -1;
	}
	for (int i = 0; i < n; i++) {
		bool f = false;
		for (int j = 1; j <= bs; j++) {
			if (InBag(j, i)) {
				if (!f) {
					bool ok = dfs(j, i, 0, u);
					if (!ok) { std::cerr << "dfs1fail" << std::endl; return false; }
					f = true;
				}
				if (u[j] != i) {
                    std::cerr << "dfs2fail " << u[j] << " " << i << std::endl; 
					return false;
				}
			}
		}
	}
	return true;
}

vector<vector<int>>& TreeDecomposition::Bags() {
	return bags;
}

Graph TreeDecomposition::Chordal() const {
	Graph ret(n);
	for (const auto& bag : bags) {
		for (int i = 0; i < (int)bag.size(); i++) {
			for (int j = i+1; j < (int)bag.size(); j++) {
				assert(bag[i] >= 0 && bag[i] < n && bag[j] >= 0 && bag[j] < n);
				ret.AddEdge(bag[i], bag[j]);
			}
		}
	}
	return ret;
}

int TreeDecomposition::nbags() const {
	return bs;
}

int TreeDecomposition::nverts() const {
	return n;
}

const std::unordered_set<int>& TreeDecomposition::Neighbors(int b) const {
	assert(b >= 1 && b <= bs);
	return tree.Neighbors(b);
}

int TreeDecomposition::CenDfs(int b, int p, int& cen) const {
	assert(b >= 1 && b <= bs);
	assert(p >= 0 && p <= bs);
	assert(cen == 0);
	int intro = 0;
	for (int nb : Neighbors(b)) {
		if (nb == p) continue;
		int cintro = CenDfs(nb, b, cen);
		intro += cintro;
		if (cintro >= n/2) {
			assert(cen);
			return intro;
		}
	}
	for (int v : bags[b]) {
		if (p == 0 || !InBag(p, v)) {
			intro++;
		}
	}
	if (intro >= n/2) {
		cen = b;
	}
	return intro;
}

int TreeDecomposition::Centroid() const {
	int cen = 0;
	CenDfs(1, 0, cen);
	assert(cen >= 1 && cen <= bs);
	return cen;
}

void TreeDecomposition::OdDes(int b, int p, int d, vector<int>& ret) const {
  assert(b >= 1 && b <= bs);
  assert(p >= 0 && p <= bs);
  assert(d >= 1);
  bool new_vs = false;
  for (int v : bags[b]) {
    if (ret[v] == 0) {
      new_vs = true;
    } else {
      assert(ret[v] <= d);
      assert(binary_search(bags[p].begin(), bags[p].end(), v));
    }
  }
  if (new_vs) {
    d++;
    for (int v : bags[b]) {
      if (ret[v] == 0) {
        ret[v] = d;
      }
    }
  }
  for (int nb : Neighbors(b)) {
    if (nb == p) continue;
    OdDes(nb, b, d, ret);
  }
}

vector<int> TreeDecomposition::GetOrd() const {
  int centroid = Centroid();
  assert(centroid >= 1 && centroid <= bs);
  vector<int> ret(n);
  OdDes(centroid, 0, 1, ret);
  for (int i = 0; i < n; i++) {
    assert(ret[i] > 0);
  }
  return ret;
}

}