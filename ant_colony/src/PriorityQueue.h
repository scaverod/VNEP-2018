#ifndef PRIORITY_QUEUE
#define PRIORITY_QUEUE

#include "heap.h"
#include <limits>

class PriorityQueue {
	const double MINIMO;
	Heap* heap_;
	double* cost_;
public: 
	PriorityQueue(int maxItems, double **cost) : 
    MINIMO(std::numeric_limits<double>::min()), cost_(new double[maxItems]())
	{
		heap_ = HeapInit(maxItems+1);
		HeapSetKeys(heap_, cost_);
    *cost = cost_;
	}
	int key_top(){
		return HeapMin(heap_);
	}
	//double value_top(){

	//}
	void pop(){
		if(!empty())
			HeapDelMin(heap_);
		else
			exit(999);
	}
	void insert(int key){
		HeapInsert(heap_,key);
	}
	bool empty(){
		return HeapSize(heap_) == 0;
	}
	void decrease(int key){
		HeapDecKey(heap_, key);
	}
	void increase(int key){
		HeapIncKey(heap_, key);
	}
	/*void remove(int key){
		double temp = cost_[key];
		cost_[key] = MINIMO;
		decrease(key);
		
		if( key != key_top() )
			exit(100);
		pop();
		cost_[key] = temp;
	}*/
	~PriorityQueue(){
    delete [] cost_;
		HeapFree(heap_);
	}
};
#endif 
