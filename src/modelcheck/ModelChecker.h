#ifndef MODELCHECKER_H
#define MODELCHECKER_H

#include <iostream>
#include <string>
#include <vector>
#include "Node.h"
#include "buildTree.h"

using namespace std;

class ModelChecker {
public:
	ModelChecker(string formulaString){
		roots.push_back(buildTree(formulaString,0));
		nTime = 0;
		this->formulaString = formulaString;
		}
	~ModelChecker(){
		for ( vector<Node*>::iterator it = roots.begin(); it != roots.end(); ++it ) {
			delete * it;
			}
		roots.clear();
		}

	
	int update(double *x, int isLast){
		roots[nTime]->update(x, isLast);
		int tf = roots[0]->evalNode();
		
		if(tf==-1){
			roots.push_back(roots[nTime]->duplicate());
			roots[nTime]->link(roots[nTime+1]);
			}
		
		nTime = nTime + 1;
		return tf;
		}
	
private:
		vector<Node*> roots;
		string formulaString;
		int nTime;
	};

#endif