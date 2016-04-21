#include <sched.h>

#include "thread.h"

#include "sts.h"

int thread_local Thread::id_ = 0;

void Thread::doWork() {
    for (int i=0; ; i++) {
        int c = STS::waitOnStepCounter(i);
        if (c<0) break; //negative task counter signals to terminate the thread
        processQueue();
    }
}

void Thread::processQueue() {
    unsigned int s = STS::getCurrentInstance()->getNumSubTasks(Thread::id_);
    while(nextSubtaskId_<s) {
        processTask();
    }
    nextSubtaskId_=0;
}

void Thread::processTask() {
    SubTask* subtask = STS::getCurrentInstance()->getSubTask(Thread::id_, nextSubtaskId_++);
    auto startWaitTime = sts_clock::now();
    ITaskFunctor *task = STS::getCurrentInstance()->getTaskFunctor(subtask->getTaskId());
    auto startTaskTime = sts_clock::now();
    subtask->waitTime_ = startTaskTime - startWaitTime;
    task->run(subtask->getRange());
    subtask->runTime_ = sts_clock::now() - startTaskTime;
    STS::getCurrentInstance()->markSubtaskComplete(subtask->getTaskId());
}

