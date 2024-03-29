#ifndef NODE_H
#define NODE_H

#include <iostream>

class Node {
	public:
		Node(){
			known = truth = 0;
			child1 = child2 = previous = next = NULL;
			time = 1;
			}
		Node(int timeL) : timeLimit(timeL) {
			known = truth = 0;
			child1 = child2 = previous = next = NULL;
			time = 1;
			}
		~Node(){
			delete child1;
			delete child2;
			}
		void setChild1(Node *node){
			this->child1 = node;
			}
		
		void setChild2(Node *node){
			this->child2 = node;
			}
		
		void setTimeLimit(int timeLim){
			this->timeLimit = timeLim;
			}
		
		virtual int evalNode(){
			//printf("-> Eval generic node\n");
			if(known) {
				return truth;
				}
			else {
				return -1;
				}
		}
		
		virtual Node* duplicate(){
			Node* copyNode = new Node(timeLimit);
			if(child1)
				copyNode->setChild1(child1->duplicate());
			if(child2)
				copyNode->setChild2(child2->duplicate());
			return copyNode;
			}
		
		void link(Node *node){
			if(child1)
				(this->child1)->link(node->child1);
			if(child2)
				(this->child2)->link(node->child2);
			this->next = node;
			(this->next)->setTime(time+1);
			node->previous = this;
		}
		
		void setTime(int t){
			time = t;
			}
		
		virtual void update(double* x, int isLast){
			//printf("-> Update generic node\n");
			this->isLast = isLast;
			if(child1)
				(this->child1)->update(x,isLast);
			if(child2)
				(this->child2)->update(x,isLast);
			}
		
	protected:
		int truth, known;
		int time, timeLimit;
		int isLast;
		Node *child1, *child2;
		Node *previous, *next;
	};
	
	
class INode : public Node {
	public:
		INode() : Node(){}
		INode(int timeL) : Node(timeL) {}
		
		int evalNode(){
			if(known) {
				return truth;
			}
			if(time==timeLimit){
				int tf = child1->evalNode();
				known  = 1;
				truth = tf;
				return tf;
			}
			else if (next){
				int tf = next->evalNode();
				if (tf == -1)
					return tf;
				else {
					known  = 1;
					truth = tf;
					return tf;
					}
				}
			else {
				return -1;
			}
		}
		
		
		Node* duplicate(){
			Node* copyNode = new INode(timeLimit);
			if(child1)
				copyNode->setChild1(child1->duplicate());
			if(child2)
				copyNode->setChild2(child2->duplicate());
			return copyNode;
		}
		
		
};	
	
	
class XNode : public Node {
	public:
		XNode() : Node(){}
		XNode(int timeL) : Node(timeL) {}
		
		int evalNode(){
			//printf("-> Eval X node time %d\n",time);
			if (isLast || (time==timeLimit)){
				known = 1;
				truth = 0;
				return 0;
				}
			else if(!next){
				known = 0;
				return -1;
				}
			else {
				//int tfx = next->evalNode();
				int tfx = ((XNode*)next)->evalNext();
				if (tfx==-1){
					known = 0;
					return -1;
					}
				else {
					known = 1;
					truth = tfx;
					return tfx;
					}
				}
		}
		
		int evalNext(){
			int tf = child1->evalNode();
			return tf;
			}
		
		Node* duplicate(){
			Node* copyNode = new XNode(timeLimit);
			if(child1)
				copyNode->setChild1(child1->duplicate());
			if(child2)
				copyNode->setChild2(child2->duplicate());
			return copyNode;
		}
		
		
};

class FNode : public Node {
	public:
	FNode() : Node(){}
	FNode(int timeL) : Node(timeL) {}
	int evalNode(){
		if(known) {
			return truth;
			}
		//printf("-> Eval F node time %d\n",time);
		//printf("-> Eval F node\n");
		int tf = (this->child1)->evalNode();
		// If we are at end or at time limit
		if((this->isLast) || ((this->timeLimit > 0) && (this->time==this->timeLimit))){
			this->known = 1;
			this->truth = tf;
			return tf;
			}
		else {
			// Current time is true
			if(tf==1){
				this->known = 1;
				this->truth = 1;
				return 1;
				}
			// Current time is not true
			else {
				// There's no next node to ask
				if(!(this->next)){
					this->known = 0;
					return -1;
					}
				// We ask the next node
				else {
					int tfx = (this->next)->evalNode();
					if(tfx==-1){
						this->known = 0;
						return -1;
						}
					else {
						this->known = 1;
						this->truth = tfx;
						return tfx;
						}
					}
				}
			}
		}
	Node* duplicate(){
		Node* copyNode = new FNode(timeLimit);
		if(child1)
			copyNode->setChild1(child1->duplicate());
		if(child2)
			copyNode->setChild2(child2->duplicate());
		return copyNode;
		}
};

class GNode : public Node {
	public:
	GNode() : Node(){}
	GNode(int timeL) : Node(timeL) {}
	int evalNode(){
		//printf("-> Eval G node\n");
		
		
		if(known) {
			return truth;
		}
		
		int tf = child1->evalNode();
		if(isLast || ((timeLimit > 0) && (time==timeLimit))){
			this->known = 1;
			this->truth = tf;
			return tf;
			}
		else {
			if(tf==0){
				this->known = 1;
				this->truth = 0;
				return 0;
				}
			else {
				// There's no next node to ask
				if(!(this->next)){
					this->known = 0;
					return -1;
					}
				// We ask the next node
				else {
					int tfx = (this->next)->evalNode();
					if(tfx==-1){
						this->known = 0;
						return -1;
						}
					else {
						this->known = 1;
						this->truth = tfx;
						return tfx;
						}
					}
				}
			}
	}
	
	Node* duplicate(){
		Node* copyNode = new GNode(timeLimit);
		if(child1)
			copyNode->setChild1(child1->duplicate());
		if(child2)
			copyNode->setChild2(child2->duplicate());
		return copyNode;
		}
};

class UNode : public Node {
	public:
	UNode() : Node(){}
	UNode(int timeL) : Node(timeL) {}
	int evalNode(){
		//printf("-> Eval U node\n");
		if(known) {
			return truth;
		}
	}
	Node* duplicate(){
		Node* copyNode = new UNode(timeLimit);
		if(child1)
			copyNode->setChild1(child1->duplicate());
		if(child2)
			copyNode->setChild2(child2->duplicate());
		return copyNode;
		}
};

class AndNode : public Node {
	public:
	AndNode() : Node(){}
	AndNode(int timeL) : Node(timeL) {}
	int evalNode(){
		//printf("-> Eval & node\n");
		if(known) {
			return truth;
		}
		int tf1 = (this->child1)->evalNode();
		if (tf1==0){
			this->known = 1;
			this->truth = 0;
			return 0;
			}
		else {
			int tf2 = (this->child2)->evalNode();
			if(tf2==0){
				this->known = 1;
				this->truth = 0;
				return 0;
				}
			else if((tf1==-1)||(tf2==-1)){
				this->known = 0;
				return -1;
				}
			else {
				this->known = 1;
				this->truth = 1;
				return 1;
				}
			}
		}
	Node* duplicate(){
		Node* copyNode = new AndNode(timeLimit);
		if(child1)
			copyNode->setChild1(child1->duplicate());
		if(child2)
			copyNode->setChild2(child2->duplicate());
		return copyNode;
	}
};

class OrNode : public Node {
public:
	OrNode() : Node(){}
	OrNode(int timeL) : Node(timeL) {}
	int evalNode(){
		//printf("-> Eval | node\n");
		if(known) {
			return truth;
		}
		int tf1 = (this->child1)->evalNode();
		if(tf1==1){
			this->known = 1;
			this->truth = 1;
			return 1;
			}
		else {
			int tf2 = (this->child2)->evalNode();
			if(tf2==1){
				this->known = 1;
				this->truth = 1;
				return 1;
				}
			else if ((tf1==-1)||(tf2==-1)){
				this->known = 0;
				return -1;
				}
			else {
				this->known = 1;
				this->truth = 0;
				return 0;
				}
			}
		}
	Node* duplicate(){
		Node* copyNode = new OrNode(timeLimit);
		if(child1)
			copyNode->setChild1(child1->duplicate());
		if(child2)
			copyNode->setChild2(child2->duplicate());
		return copyNode;
	}
};

class NotNode : public Node {
	public:
	NotNode() : Node(){}
	NotNode(int timeL) : Node(timeL) {}
	int evalNode(){
		//printf("-> Eval ! node\n");
		if(known) {
			return truth;
		}
		int tf = (this->child1)->evalNode();
		if(tf==-1){
			this->known = false;
			return -1;
			}
		else {
			this->known = true;
			this->truth = !tf;
			return !tf;
			}
		}
	Node* duplicate(){
		Node* copyNode = new NotNode(timeLimit);
		if(child1)
			copyNode->setChild1(child1->duplicate());
		if(child2)
			copyNode->setChild2(child2->duplicate());
		return copyNode;
	}
};

class ImplNode : public Node {
	public:
		ImplNode() : Node(){}
		ImplNode(int timeL) : Node(timeL) {}
		int evalNode(){
			//printf("-> Eval ~ node\n");
			if(known) {
				return truth;
			}
			int tf1 = (this->child1)->evalNode();
			if(tf1==0){
				this->known = 1;
				this->truth = 1;
				return 1;
				}
			else {
				int tf2 = (this->child2)->evalNode();
				if(tf2==1){
					this->known = 1;
					this->truth = 1;
					return 1;
					}
				else if ((tf1==-1)||(tf2==-1)){
					this->known = 0;
					return -1;
					}
				else {
					this->known = 1;
					this->truth = 0;
					return 0;
					}
				}
			}
		Node* duplicate(){
			Node* copyNode = new ImplNode(timeLimit);
			if(child1)
				copyNode->setChild1(child1->duplicate());
			if(child2)
				copyNode->setChild2(child2->duplicate());
			return copyNode;
		}
};



class AtomicNode : public Node {
	public:
		AtomicNode(int v, double l, double u) : varID(v), lb(l), ub(u) {}
		void update(double* x, int isLast){
			//printf("--> Update atomic node\n");
			known = 1;
			if((x[varID-1] >= lb) && ((ub==-1) || (x[varID-1] <= ub))){
				truth = 1;
				}
			else {
				truth = 0;
				}
			}
		int evalNode(){
			return truth;
			}

		Node* duplicate(){
			Node* copyNode = new AtomicNode(varID,lb,ub);
			return copyNode;
		}
	private:
		int varID;
		double lb, ub;
};


class Terminator : public Node {
	int evalNode(){}
};

#endif