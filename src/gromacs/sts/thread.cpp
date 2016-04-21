#include "thread.h"

#include "sts.h"

#include <string>

int thread_local Thread::id_ = 0;

void Thread::doWork() {
    for (int i=0; ; i++) {
        int c = STS::waitOnStepCounter(i);
        if (c<0) break; //negative task counter signals to terminate the thread
        processQueue();
    }
}

void Thread::processQueue() {
    STS* sts = STS::getCurrentInstance();
    while(sts->runNextSubTask());
}
