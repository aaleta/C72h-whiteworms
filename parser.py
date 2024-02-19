import argparse


def parse_args(args):
    parser = argparse.ArgumentParser()

    # Type of analysis
    parser.add_argument(
        "-a",
        "--analysis-type",
        type=str,
        help="one of [protected|botnet|threshold]"
    )

    # Simulation parameters
    parser.add_argument(
        "-n",
        "--network",
        type=str,
        help="name of the graph to be used"
    )

    parser.add_argument(
        '-bb',
        '--beta_b',
        type=float,
        help='black worm beta value'
    )

    parser.add_argument(
        '-bw',
        '--beta_w',
        type=float,
        help='white worm beta value'
    )

    parser.add_argument(
        "-e",
        "--epsilon",
        type=float,
        help='ethics rate'
    )

    parser.add_argument(
        "-g",
        "--gamma",
        type=float,
        help='prompted transition rate'
    )

    parser.add_argument(
        "-m",
        "--mu",
        type=float,
        help='forced fix rate'
    )

    parser.add_argument(
        "-i",
        "--iterations",
        type=int,
        default=100,
        help='number of iterations'
    )

    parser.add_argument(
        "-nb",
        "--n_black",
        type=int,
        default=1,
        help='number of black worms'
    )

    parser.add_argument(
        "-nw",
        "--n_white",
        type=int,
        default=1,
        help='number of white worms'
    )

    parser.add_argument(
        "-th",
        "--threshold",
        type=float,
        default=0,
        help='fraction of nodes with black worm'
    )

    return parser.parse_args(args)
