
import argparse
from inversion_tool.core.inversion_tools import Inverter

def main():

    parser = argparse.ArgumentParser(description="Run the simulation")

    # Add any arguments or options you need for your simulation
    parser.add_argument("run_dir", type=str, help="Path to the run directory")
    parser.add_argument("--mode", type=str, help="Mode of simulation. Either inversion or grid")

    args = parser.parse_args()

    # Call the function to run the simulation with the provided arguments

    inv = Inverter(args.run_dir)

    if args.mode == 'grid':
        print("running in grid mode")
        a = inv.run_grid()
    else:
        a = inv.run()



if __name__ == "__main__":
    main()