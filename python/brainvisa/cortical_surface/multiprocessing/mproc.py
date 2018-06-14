from __future__ import print_function

import multiprocessing, os
import six.moves.queue as queue

def os_system ( primal_command, sujet, keyword ) :
    os.system ( primal_command )
    return 'OK'

class Worker(multiprocessing.Process):

    def __init__(self, work_queue, result_queue, function):
        multiprocessing.Process.__init__(self)
        self.work_queue = work_queue
        self.result_queue = result_queue
        self.kill_received = False
        self.function = function

    def run (self) :
        while not self.kill_received:
            try:
                job = self.work_queue.get_nowait()
            except queue.Empty:
                break
            print('jobs:', len(job))
            result = self.function ( *job )
            self.result_queue.put ( (job[-1], job[-2], result) )

def MultiProcExecute ( function, jobs, num_processes = 8 ) :

    # load up work queue
    work_queue = multiprocessing.Queue()
    for job in jobs:
        work_queue.put(job)

    # create a queue to pass to workers to store the results
    result_queue = multiprocessing.Queue()

    # spawn workers
    for i in range(num_processes):
        worker = Worker ( work_queue, result_queue, function )
        worker.start()

    # collect the results off the queue
    results = []
    while len(results) < len(jobs):
        stillRunning = True
        while ( stillRunning ) :
            try :
                result = result_queue.get()
                results.append(result)
                stillRunning = False
            except IOError:
                pass

    return results