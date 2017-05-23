import threading
import time


class MySimulator (threading.Thread):

    def __init__(self):
        threading.Thread.__init__(self)
        self.data = 0

    def run(self):
        print('Starting simulation loop')
        T = 0
        while T < 10:
            mylock.acquire()
            try:
                time.sleep(1.0)
                self.data += 1
                T += 0.1
                print('Stepping done')
            finally:
                mylock.release()
        print('Stopped simulation loop')


class MyRender(threading.Thread):

    def __init__(self, simulator):
        threading.Thread.__init__(self)
        self.simulator = simulator

    def run(self):
        print('Starting render loop')
        while True:
            mylock.acquire()
            print('  Starting rendering frame')
            print("  Drawn simulation data = ", self.simulator.data)
            print('  Stopped rendering')
            mylock.release()
        print('Stopped render loop')


if __name__ == '__main__':
    mylock = threading.Lock()

    simulator = MySimulator()
    simulator.start()

    render = MyRender(simulator)
    render.start()

