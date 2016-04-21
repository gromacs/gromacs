#include "task.h"

bool SubTask::run() {
    isDone_ = task_->run(range_, timeData_);
    return isDone_;
}
const Task* SubTask::getTask() const {
    return task_;
}
bool SubTask::isDone() const {
    return isDone_;
}
void SubTask::setDone(bool isDone) {
    isDone_ = isDone;
}
bool SubTask::isReady() const {
    return task_->isReady();
}
