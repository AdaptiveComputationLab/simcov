import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--dim" , nargs='+', type=int, help="Dimensions: x y z", default=[15000, 15000, 1])
    parser.add_argument("--whole-lung-dim" , nargs='+', type=int, help="Whole lung dimensions: x y z", default=[48000, 40000, 20000])
    parser.add_argument("--timesteps" , type=int, help="Number of timesteps (assuming 1 min per timestep)", default=4000)
    parser.add_argument("--infection-coords" , type=str, help="Location (distribution) of initial infections", default="uniform:1")
    parser.add_argument("--initial-infection" , type=int, help="Number of virions at initial infection locations", default=1000)
    parser.add_argument("--incubation-period" , type=int, help="Average number of time steps to expressing virions after cell is infected", default=480)
    parser.add_argument("--apoptosis-period" , type=int, help="Average number of time steps to death after apoptosis is induced", default=180)
    parser.add_argument("--expressing-period" , type=int, help="Average number of time steps to death after a cell starts expresssing", default=900)
    parser.add_argument("--infectivity" , type=float, help="Factor multiplied by number of virions to determine probability of infection", default=0.001)
    parser.add_argument("--infectivity-multiplier" , type=float, help="Multiplier reducing infectivity where inflammatory signal is present", default=1.0)
    parser.add_argument("--virion-production" , type=float, help="Number of virions produced by expressing cell each time step", default=1.1)
    parser.add_argument("--virion-production-multiplier" , type=float, help="Multiplier reducing virion production rate where inflammatory signal is present", default=1.0)
    parser.add_argument("--virion-clearance" , type=float, help="Fraction by which virion count drops each time step", default=0.004)
    parser.add_argument("--virion-diffusion" , type=float, help="Fraction of virions that diffuse into all neighbors each time step", default=0.15)
    parser.add_argument("--chemokine-production" , type=float, help="Amount of chemokine produced by expressing cells each time step", default=1.0)
    parser.add_argument("--chemokine-decay" , type=float, help="Amount by which chemokine concentration drops each time step", default=0.01)
    parser.add_argument("--chemokine-diffusion" , type=float, help="Fraction of chemokine concentration that diffuses into all neighbors each time step", default=1.0)
    parser.add_argument("--air-diffusion" , type=float, help="Fraction of chemokine concentration that diffuses into all neighbors through air each time step", default=1.0)
    parser.add_argument("--min-chemokine" , type=float, help="Minimum chemokine concentration that triggers a T cell", default=1e-6)
    parser.add_argument("--antibody-factor" , type=int, help="Impact of antibodies; multiplier for virion decay (setting to 1 means this has no effect)", default=1)
    parser.add_argument("--antibody-period" , type=int, help="Number of time steps before antibodies start to be produced", default=5760)
    parser.add_argument("--tcell-generation-rate" , type=int, help="Number of tcells generated at each timestep for the whole lung", default=105000)
    parser.add_argument("--tcell-initial-delay" , type=int, help="Number of time steps before T cells start to be produced", default=10080)
    parser.add_argument("--tcell-vascular-period" , type=int, help="Average number of time steps to death for a T cell in the vasculature", default=5760)
    parser.add_argument("--tcell-tissue-period" , type=int, help="Average number of time steps to death after a T cell extravasates", default=1440)
    parser.add_argument("--tcell-binding-period" , type=int, help="Number of time steps a T cell is bound to an epithelial cell when inducing apoptosis", default=10)
    parser.add_argument("--max-binding-prob" , type=int, help="Max probability of a T cell binding to an infected cell in one time step", default=1)
    parser.add_argument("--tcells-follow-gradient" , type=str, help="T cells in tissue follow the chemokine gradient", default="false")
    parser.add_argument("--seed" , type=int, help="Random seed", default=29)
    parser.add_argument("--sample-period" , type=int, help="Number of timesteps between samples (set to 0 to disable sampling)", default=10080)
    parser.add_argument("--sample-resolution" , type=int, help="Resolution for sampling", default=1)
    parser.add_argument("--max-block-dim" , type=int, help="Max. block dimension - larger means more locality but worse load balance. Set to 0 for largest possible", default=10)
    parser.add_argument("--output" , type=str, help="Output directory (automatically generated)", default="results")

    args = vars(parser.parse_args())
    for key in args:
        print(key.replace("_", "-"), end='=')
        if isinstance(args[key], list):
            args[key] = list(map(str, args[key]))
            print(' '.join(args[key]))
        else:
            print(args[key])