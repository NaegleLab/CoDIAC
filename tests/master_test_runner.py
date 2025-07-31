# run_all_tests.py

import run_setup_task_tests as test_part1
import run_adjfile_tests as test_part2
import run_contact_maps_tests as test_part3

def run_all():
    print("Running test part 1...")
    test_part1.main()
    
    print("Running test part 2...")
    test_part2.main()
    
    print("Running test part 3...")
    test_part3.main()

if __name__ == "__main__":
    run_all()
