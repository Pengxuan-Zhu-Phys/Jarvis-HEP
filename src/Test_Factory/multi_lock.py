#!/usr/bin/env python3 

import multiprocessing, time 

def worker(dict, lock):
    while True:
        lock.acquire(timeout=5)
        number = dict.get("ticket")
        if number > 0:
            time.sleep(1)
            number -= 1 
            print("{}, ticket = {}".format(multiprocessing.current_process().name, number))
            dict.update({"ticket": number})
        else:
            break
        lock.release()

def main():
    lock = multiprocessing.Lock()
    manager = multiprocessing.Manager()
    manager_dict = manager.dict(ticket=5)
    print(manager_dict)
    jobs = [multiprocessing.Process(target=worker, args=(manager_dict, lock,), name="saller: {}".format(item)) for item in range(10)]
    for proc in jobs:
        proc.start()
    for proc in jobs:
        proc.join()

    print("All done, manager_dict => {}".format(manager_dict))

if __name__ == "__main__":
    main()