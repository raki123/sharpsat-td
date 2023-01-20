#include "treewidth.hpp"

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

namespace sspp {

namespace decomp {

TreeDecomposition Treedecomp(const Graph& graph, double time, string tmp_dir) {
	int n = graph.n();
	if (n == 0) {
		TreeDecomposition dec(0, 0);
		return dec;
	}
	if (n == 1) {
		TreeDecomposition dec(1, 1);
		dec.SetBag(1, {0});
		return dec;
	}
	if (time < 0.099) {
		TreeDecomposition dec(1, n);
		vector<int> all;
		for (int i = 0; i < n; i++) {
			all.push_back(i);
		}
		dec.SetBag(1, all);
		return dec;
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
		int claim_width = 0;
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
		TreeDecomposition dec(0, 0);
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
				int bs,nn;
				ss>>bs>>claim_width>>nn;
				assert(nn == n);
				claim_width--;
				dec = TreeDecomposition(bs, nn);
			} else if (tmp == "b") {
				int bid;
				ss>>bid;
				vector<int> bag;
				int v;
				while (ss>>v) {
					bag.push_back(v-1);
				}
				dec.SetBag(bid, bag);
			} else {
				int a = stoi(tmp);
				int b;
				ss>>b;
				dec.AddEdge(a, b);
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
		assert(dec.Width() <= claim_width);
		cout << "c o width " << dec.Width() << endl;
		assert(dec.Verify(graph));
		return dec;
	}
}
} // namespace decomp
} // namespace sspp