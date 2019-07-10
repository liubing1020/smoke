#ifndef BUILDTREE_H
#define BUILDTREE_H

#include "Node.h"
#include <string>
#include <iostream>

using namespace std;

int* getBal(string,char,char);
int isBalanced(string,char,char);


Node* buildTree(string fStr,int timeLimit){
	char firstChar = fStr[0];
	char lastChar = fStr[fStr.length()-1];
	Node* root = NULL;
	
	//std::cout << "Str: " << fStr << endl;
	//printf("Str: %s\n",fStr.c_str());
	//printf("Time lim passed: %d\n", timeLimit);
	//printf("First: %c, Last: %c, \n",firstChar,lastChar);
	//printf("Time limit: %d\n",timeLimit);
	
	while (fStr[0]==' '){
		fStr = fStr.substr(1);
		}
	while (fStr[fStr.size()-1]==' '){
		fStr = fStr.substr(0,fStr.length()-1);
		}

	firstChar = fStr[0];
	lastChar = fStr[fStr.length()-1];


	switch(firstChar){
		case '!':
			root = new NotNode(timeLimit);
			//printf("-> Making ! node...\n");
			root->setChild1(buildTree(fStr.substr(1),timeLimit));
			return root;
		case 'X':
			root = new XNode(timeLimit);
			//printf("-> Making X node...\n");
			root->setChild1(buildTree(fStr.substr(1),timeLimit));
			return root;
			
		case 'F':
			root = new FNode(timeLimit);
		case 'G':
			if(!root) root = new GNode(timeLimit);
		case 'I':
			if(!root) root = new INode(timeLimit);
			
			if(fStr[1]=='{'){
				if(fStr[2]=='}'){
					root->setChild1(buildTree(fStr.substr(3),timeLimit));
					}
				else {
					int pos = fStr.find("}");
					int timeLim = atoi(fStr.substr(2,pos-2).c_str());
					if(timeLim > 0){
						root->setTimeLimit(timeLim);
						root->setChild1(buildTree(fStr.substr(pos+1),timeLim));
						}
					}
				}
			else {
				root->setChild1(buildTree(fStr.substr(1),timeLimit));
				}

			return root;
		case '[':
			if(lastChar == ']'){
				int varID;
				double lb, ub;
				sscanf(fStr.c_str(),"[%d,%lf,%lf]",&varID,&lb,&ub);
				root = new AtomicNode(varID,lb,ub);
				//printf("-> Making atomic node...\n");
				return root;
				}
			break;
		case '(':
			if(lastChar == ')'){
				int* parBal = getBal(fStr,'(',')');
				int balOK = 1;
				for(int j=0;j<fStr.length();j++){
					//printf("%d,",parBal[j]);
					if(parBal[j]<1){
						balOK = 0;
						break;
						}
					}
				if(balOK){
					root = buildTree(fStr.substr(1,fStr.length()-2),timeLimit);
					return root;
					}
				}
		}
	
	string strRight = fStr;
	string strLeft;
	int tokenPos;
	char token;
	while(true) {
		tokenPos = strRight.find_first_of("&|~UR");
		token = strRight[tokenPos];
		if (tokenPos <= 0){
			break;
			}
		
		strRight = strRight.substr(tokenPos+1);
		strLeft = fStr.substr(0,fStr.size()-strRight.size()-1);
		//printf("\nStr left: %s\n",strLeft.c_str());
		//printf("Str right: %s\n",strRight.c_str());
		//printf("Bal: %d\n",isBalanced(strLeft,'(',')'));

		if(isBalanced(strLeft,'(',')')==0){		
			switch(token){
				case '&':
					root = new AndNode(timeLimit);
					//printf("-> Making & node...\n");
					break;
				case '|':
					root = new OrNode(timeLimit);
					//printf("-> Making | node...\n");
					break;
				case '~':
					root = new ImplNode(timeLimit);
					//printf("-> Making ~ node...\n");
					break;
				case 'U':
					root = new UNode(timeLimit);
					//printf("-> Making U node...\n");
					break;
				}

			root->setChild1(buildTree(strLeft,timeLimit));
			root->setChild2(buildTree(strRight,timeLimit));
			return root;
			}
		}
	}

int isBalanced(string str, char openBrac, char closeBrac){
	int isBal = 0;
	isBal += (str[0]==openBrac);
	for(int i=1;i<str.length();i++){
		isBal += (str[i]==openBrac) - (str[i-1]==closeBrac);
		}
	isBal -= (str[str.length()-1]==closeBrac);
	return isBal;
	}

int* getBal(string str,char openBrac, char closeBrac){
	int* balVec = new int[str.length()];
	
	balVec[0] = (str[0]==openBrac);
	for(int i=1;i<str.length();i++){
		balVec[i] = balVec[i-1] + (str[i]==openBrac) - (str[i-1]==closeBrac);
		}
	return balVec;
	}

#endif