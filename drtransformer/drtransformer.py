#!/usr/bin/env python
#
# DrTransformer -- cotranscriptional folding.
#
import logging
import os
import sys
import math
import argparse
import numpy as np
from packaging import version

# Debugging only
from datetime import datetime

import RNA
from . import __version__, _MIN_VRNA_VERSION 
from .utils import parse_vienna_stdin, get_tkn_simulation_files
from .landscape import TrafoLandscape
from .pathfinder import top_down_coarse_graining

def restricted_float(x):
    y = float(x)
    if y < 0.0 or y > 1.0:
        raise argparse.ArgumentTypeError(f"{x} not in range (0.0, 1.0)")
    return y

def sorted_trajectories(nlist, tfile, plot_cgm = None, mapping = None):
    """ Yields the time course information using a treekin output file.

    Args:
      nlist (list): a list of nodes sorted by their energy
      tfile (str): treekin-output file name.
      softmap (dict, optional): A mapping to transfer occupancy between
        states. Likely not the most efficient implementation.

    Yields:
      list: ID, time, occupancy, structure, energy
    """
    with open(tfile) as tkn:
        for line in tkn:
            if line[0] == '#':
                continue
            course = list(map(float, line.strip().split()))
            time = course[0]
            if plot_cgm:
                # Combine all results that belong together.
                for e, occu in enumerate(course[1:]):
                    if nlist[e] not in plot_cgm:
                        continue
                    node = nlist[e]
                    for n in plot_cgm[node]:
                        ci = nlist.index(n) + 1
                        occu += course[ci]/len(mapping[n])
                    yield time, nlist[e], occu 
            else:
                for e, occu in enumerate(course[1:]):
                    # is it above visibility threshold?
                    yield time, nlist[e], occu 
    return

def parse_model_details(parser):
    """ ViennaRNA Model Details Argument Parser.  """
    model = parser.add_argument_group('ViennaRNA model details')

    model.add_argument("-T", "--temp", type = float, default = 37.0, metavar = '<flt>',
        help = 'Rescale energy parameters to a temperature of temp C.')

    model.add_argument("-4", "--noTetra", action = "store_true",
        help = 'Do not include special tabulated stabilizing energies for tri-, tetra- and hexaloop hairpins.')

    model.add_argument("-d", "--dangles", type = int, default = 2, metavar = '<int>',
        help = 'How to treat "dangling end" energies for bases adjacent to helices in free ends and multi-loops.')

    model.add_argument("--noGU", action = "store_true",
        help = 'Do not allow GU/GT pairs.')

    model.add_argument("--noClosingGU", action = "store_true",
        help = 'Do not allow GU/GT pairs at the end of helices.')

    model.add_argument("-P", "--paramFile", action = "store", default = None, metavar = '<str>',
        help = 'Read energy parameters from paramfile, instead of using the default parameter set.')

def parse_drtrafo_args(parser):
    """ A collection of arguments that are used by DrTransformer """

    environ = parser.add_argument_group('DrTransformer dependencies')
    output = parser.add_argument_group('DrTransformer output')
    trans = parser.add_argument_group('Transcription parameters')
    algo  = parser.add_argument_group('DrTransformer algorithm')

    ###################
    # Default options #
    ###################
    parser.add_argument('--version', action = 'version', 
            version = '%(prog)s ' + __version__)
    parser.add_argument("-v", "--verbose", action = 'count', default = 0,
            help = """Track process by writing verbose output to STDOUT during
            calculations. Use --logfile if you want to see *just* verbose
            information via STDOUT.""")

    ########################
    # DrTransformer output #
    ########################
    output.add_argument("--name", default = '', metavar = '<str>',
            help = """Name your output files, name the header of your plots, etc.
            this option overwrites the fasta-header.""")

    output.add_argument("--stdout", default = None, action = 'store',
            choices = ('log', 'drf', 'OFF'),
            help = """Choose STDOUT formats to follow the cotranscriptional
            folding progress in real time: *log*: a human readable output
            format.  *drf*: DrForna visalization input format. *OFF*: actively
            suppress output. The default (None) switches between *OFF*, if --logfile is
            specified, or *log* otherwise.""")

    output.add_argument("--logfile", action = "store_true",
            help = """Write verbose information to a file:
            {--outdir}/{--name}.log""")

    output.add_argument("--tmpdir", default = '', action = 'store', metavar = '<str>',
            help = """Specify path for storing landscape files for debugging. These
            files will not be removed when the program terminates. """)

    output.add_argument("--outdir", default = '', action = 'store', metavar = '<str>',
            help = """Place regular output files, into this directory. Creates
            the directory if it does not exist. """)

    output.add_argument("--no-timecourse", action = "store_true",
            help = """Do not produce the time-course file (outdir/name.drf).""")

    output.add_argument("--plot-minh", type = float, default = None, metavar = '<flt>',
            help = """Coarsify the output based on a minimum barrier height. In contrast
            to t-fast, this does *not* speed up the simulation.""")

    output.add_argument("--draw-graphs", action = "store_true",
            #help = """Export every landscape as json file. Uses --tempdir. """)
            help = argparse.SUPPRESS)

    output.add_argument("--t-lin", type = int, default = 30, metavar = '<int>',
            help = """Evenly space output *--t-lin* times during transcription on a linear time scale.""")

    output.add_argument("--t-log", type = int, default = 300, metavar = '<int>',
            help = """Evenly space output *--t-log* times after transcription on a logarithmic time scale.""")

    ############################
    # Transcription parameters #
    ############################
    trans.add_argument("--t-ext", type = float, default = 0.02, metavar = '<flt>',
            help = """Transcription speed, i.e. time per nucleotide extension
            [seconds per nucleotide].""")

    trans.add_argument("--t-end", type = float, default = 60, metavar = '<flt>',
            help = "Post-transcriptional simulation time [seconds].")

    trans.add_argument("--pause-sites", nargs='+', metavar='<int>=<flt>',
            help="""Transcriptional pausing sites.  E.g. \"--pause-sites 82=2e3
            84=33\" alters the simulation time at nucleotides 82 and 84 to 2000
            and 33 seconds, respectively. """)

    trans.add_argument("--start", type = int, default = 1, metavar = '<int>',
            help = "Start transcription at this nucleotide.")

    trans.add_argument("--stop", type = int, default = None, metavar = '<int>',
            help = "Stop transcription before this nucleotide")

    ###########################
    # DrTransformer algorithm #
    ###########################
    algo.add_argument("--o-prune", type = restricted_float, default = 0.1, metavar = '<flt>',
            help = """Occupancy threshold to prune structures from the 
            network. The structures with lowest occupancy are removed until
            at most o-prune occupancy has beem removed from the total population. """)

    algo.add_argument("--t-fast", type = float, default = 0.001, metavar = '<flt>',
            help = """Folding times faster than --t-fast are considered
            instantaneous.  Structural transitions that are faster than
            --t-fast are considerd part of the same macrostate. Directly
            translates to an energy barrier separating conforations using:
            dG = -RT*ln((1/t-fast)/k0). None: t-fast = 1/k_0 """)

    algo.add_argument("--minh", type = float, default = None, metavar = '<flt>',
            # An alternative to specify --t-fast in terms of a barrier height.
            help = argparse.SUPPRESS)

    algo.add_argument("--force", action = "store_true", 
            # Enforce a setting against all warnings.
            help = argparse.SUPPRESS)

    # NOTE: findpath width is flexible by default, easy to implement though.
    #algo.add_argument("--findpath-search-width", type = int, default = 0, metavar = '<int>',
    #        help = """Search width for the *findpath* heuristic. Higher values
    #        increase the chances to find energetically better transition state
    #        energies.""")

    algo.add_argument("--min-fraying", type = int, default = 6, metavar = '<int>',
            help = """Minimum number of freed bases during helix fraying.
            Fraying helices can vary greatly in length, starting with at
            least two base-pairs. This parameter defines the minimum amount of
            bases freed by helix fraying. For example, 6 corresponds to a
            stack of two base-pairs and a loop region of 2 nucleotides. If less
            bases are freed and there exists a nested stacked helix, this helix
            is considered to fray as well.""")

    algo.add_argument("--k0", type = float, default = 2e5, metavar = '<flt>',
            help = """Arrhenius rate constant (pre-exponential factor). Adjust
            this constant of the Arrhenius equation to relate free energy
            changes to experimentally determined folding time [seconds].""")
    return

def write_output(data, stdout = False, fh = None):
    # Helper function to print data to filehandle, STDOUT, or both.
    if stdout:
        sys.stdout.write(data)
    if fh:
        fh.write(data)
    return

def set_handle_verbosity(h, v):
    if v == 0:
        h.setLevel(logging.WARNING)
    elif v == 1:
        h.setLevel(logging.INFO)
    elif v == 2:
        h.setLevel(logging.DEBUG)
    elif v >= 3:
        h.setLevel(logging.NOTSET)

def main():
    """ DrTransformer - cotranscriptional folding. 
    """
    parser = argparse.ArgumentParser(
        formatter_class = argparse.ArgumentDefaultsHelpFormatter,
        description = 'DrTransformer: RNA folding kinetics during transcription.')
    parse_drtrafo_args(parser)
    parse_model_details(parser)
    args = parser.parse_args()

    # ~~~~~~~~~~~~~
    # Logging Setup 
    # ~~~~~~~~~~~~~
    title = 'DrTransformer: RNA folding kinetics during transcription.'
    logger = logging.getLogger('drtransformer')
    logger.setLevel(logging.DEBUG)

    banner = "{} {}".format(title, __version__)
    ch = logging.StreamHandler()
    formatter = logging.Formatter('# %(levelname)s - %(message)s')
    set_handle_verbosity(ch, args.verbose)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info(banner)

    (name, fullseq) = parse_vienna_stdin(sys.stdin, chars='ACGUNacgun')

    if args.plot_minh:
        args.plot_minh = int(round(args.plot_minh*100))
    # Adjust arguments, prepare simulation
    if args.name == '':
        args.name = name
    else:
        name = args.name

    if os.path.split(args.name)[0]:
        raise SystemExit('ERROR: Argument "--name" must not contain filepath.')

    if args.outdir:
        if not os.path.exists(args.outdir):
            os.makedirs(args.outdir)
        filepath = args.outdir + '/' + args.name
    else:
        filepath = args.name

    dfh = open(filepath + '.drf', 'w') if not args.no_timecourse else None
    lfh = open(filepath + '.log', 'w') if args.logfile else None

    # Adjust filehandle-stuff
    if args.stdout is None and lfh is None:
        args.stdout = 'log'

    if args.tmpdir:
        _tmpdir = args.tmpdir
        if not os.path.exists(_tmpdir):
            os.makedirs(_tmpdir)

    # Adjust simulation parameters
    _RT = 0.61632077549999997
    if args.temp != 37.0:
        kelvin = 273.15 + args.temp
        _RT = (_RT / 310.15) * kelvin

    if args.stop is None:
        args.stop = len(fullseq) + 1

    if args.t_end < args.t_ext:
        raise SystemExit('ERROR: Conflicting Settings: ' + \
                'Arguments must be such that "--t-end" >= "--t-ext"')

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Adjust the simulation window for treekin:
    #
    #   k0
    #   2e5  atu/s    50 nuc/s               1e-inf /s
    #   4000 /ext     1  nuc/ext             1e-inf /ext
    #   |------$------|-------------------$----------->
    #          k-fast                     k-slow
    #          --minh                     --maxh
    #          |------|---------------|--->
    #          t0     t_ext 0.02 s    t_end = 86400 s
    #   <----->    simulation             <---------->
    #   instant                             rate = 0
    #
    # (1) The rate of a spontaneous folding event
    # (>= k-fast) has to be much faster than the
    # rate of chain elongation (kx).
    #
    # (2) The rate of an effectively 0 folding
    # event (< k-slow) has to be much slower than
    # the rate of chain elongation (kx), it
    # should also be much slower than the
    # post-transcriptional simulation time --t-end
    #
    # Parameters:
    # k0 = maximum folding rate /s
    # t-ext = time for chain elongation
    # t-end = post-transcriptional simulation time
    # t0 = first simulation output time (<< t8)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

    if args.minh:
        logger.warning('Overwriting t-fast parameter.')
        args.t_fast = 1/(args.k0 * math.exp(-args.minh/_RT))
    else:
        args.minh = max(0, -_RT * math.log(1 / args.t_fast / args.k0))
    logger.info(f'--t-fast: {args.t_fast} s => {args.minh} kcal/mol barrier height ' + 
                f'and {1/args.t_fast} /s rate at k0 = {args.k0}')

    if not args.force and args.t_fast and args.t_fast * 10 > args.t_ext:
        raise SystemExit('ERROR: Conflicting Settings: ' + 
                'Arguments must be such that "--t-fast" * 10 > "--t-ext".\n' + 
                '       => An instant folding time must be at least 10x shorter than ' +
                'the time of nucleotide extension. You may use --force to ignore this setting.')

    if version.parse(RNA.__version__) < version.parse(_MIN_VRNA_VERSION):
        raise VersionError('ViennaRNA', RNA.__version__, _MIN_VRNA_VERSION)

    ############################
    # ~~~~~~~~~~~~~~~~~~~~~~~~ #
    # Start with DrTransformer #
    # ~~~~~~~~~~~~~~~~~~~~~~~~ #
    ############################

    if args.paramFile:
        RNA.read_parameter_file(args.paramFile)

    # Set model details.
    vrna_md = RNA.md()
    vrna_md.noLP = 1
    vrna_md.logML = 0
    vrna_md.temperature = args.temp
    vrna_md.dangles = args.dangles
    vrna_md.special_hp = not args.noTetra
    vrna_md.noGU = args.noGU
    vrna_md.noGUclosure = args.noClosingGU

    # Write logging output
    if args.stdout == 'log' or lfh:
        fdata  = "# File generated using DrTransformer v{}\n".format(__version__)
        fdata += "#\n"
        fdata += "# >{}\n# {} \n".format(name, fullseq)
        fdata += "#\n"
        fdata += "# Co-transcriptional folding parameters:\n"
        fdata += "# --t-ext: {} sec\n".format(args.t_ext)
        fdata += "# --t-end: {} sec\n".format(args.t_end)
        fdata += "# --start: {}\n".format(args.start)
        fdata += "# --stop: {}\n".format(args.stop)
        fdata += "#\n"
        fdata += "# Algorithm parameters:\n"
        fdata += "# --o-prune: {}\n".format(args.o_prune)
        fdata += "# --t-fast: {} sec\n".format(args.t_fast)
        #fdata += "# --findpath-search-width: {}\n".format(args.findpath_search_width)
        fdata += "# --min-fraying: {} nuc\n".format(args.min_fraying)
        fdata += "# --k0: {}\n".format(args.k0)
        fdata += "#\n"
        fdata += "# ViennaRNA model details:\n"
        fdata += "# --temp: {} C\n".format(args.temp)
        fdata += "# --dangles: {}\n".format(args.dangles)
        fdata += "# --paramFile: {}\n".format(args.paramFile)
        fdata += "# --noTetra: {}\n".format(args.noTetra)
        fdata += "# --noGU: {}\n".format(args.noGU)
        fdata += "# --noClosingGU: {}\n".format(args.noClosingGU)
        fdata += "#\n"
        fdata += "#\n"
        fdata += "# Results:\n"
        fdata += "# Tanscription Step | Energy-sorted structure count | Structure | Energy "
        fdata += "| [Occupancy-t0 Occupancy-t8] | Structure ID (-> Plotting ID)\n"
        write_output(fdata, stdout=(args.stdout == 'log'), fh = lfh)

    # Write DrForna output
    if args.stdout == 'drf' or dfh:
        # Dictionary to store time course data.
        all_courses = dict()
        fdata = "id time conc struct energy\n"
        write_output(fdata, stdout = (args.stdout == 'drf'), fh = dfh)

    # initialize a directed conformation graph
    TL = TrafoLandscape(fullseq, vrna_md)
    TL.k0 = args.k0
    TL.minh = int(round(args.minh*100)) if args.minh else 0
    TL._transcript_length = args.start - 1

    psites = np.full(args.stop - args.start + 1, args.t_ext, dtype = float)
    if args.pause_sites:
        for term in args.pause_sites:
            site, pause = term.split('=')
            psites[int(site)] = float(pause)
    tr_end = sum(psites)
    #import statprof
    #statprof.start()

    #############
    # ~~~~~~~~~ #
    # Main loop #
    # ~~~~~~~~~ #
    #############
    for tlen in range(args.start, args.stop):
        time = TL.total_time

        logger.info(f'** Transcription step {tlen} **')
        logger.info(f'Before expansion:      {len(list(TL.active_local_mins)):3d} active lmins, {len(list(TL.active_nodes)):3d} active structures, {len(list(TL.inactive_nodes)):3d} inactive structures.')
        itime = datetime.now()

        # Get new nodes and connect them.
        nn, on = TL.expand(mfree = args.min_fraying)
        logger.info(f'After expansion:                         {len(list(TL.active_nodes)):3d} active structures, {len(list(TL.inactive_nodes)):3d} inactive structures.' +
                    f' (Found {len(nn)} new nodes and revisited {len(on)} pruned nodes.)')
        etime = datetime.now()

        cn, ce = TL.get_coarse_network()
        logger.info(f'After coarse graining: {len(list(TL.active_local_mins)):3d} active lmins, {len(list(TL.active_nodes)):3d} active structures, {len(list(TL.inactive_nodes)):3d} inactive structures.' +
                    f' (Simulation network size: nodes = {cn}, edges = {ce}.)')
        ctime = datetime.now()

        # Adjust the lenght of the lin-time simulation:
        t0, t1 = 0, psites[tlen]
        # Adjust the lenght of the log-time simulation:
        t8 = t1 + args.t_end if tlen == args.stop - 1 else t1 + sum(psites[tlen+1:])

        if np.isclose(t0, t1) or np.isclose(t1, t8): # only lin or log-part!
            if np.isclose(t1, t8):
                times = np.array(np.linspace(t0, t1, args.t_lin))
            else:
                times = np.array(np.logspace(np.log10(t1), np.log10(t8), num=args.t_log))
        else:
            lin_times = np.array(np.linspace(t0, t1, args.t_lin))
            log_times = np.array(np.logspace(np.log10(t1), np.log10(t8), num=args.t_log))
            log_times = np.delete(log_times, 0)
            times = np.concatenate([lin_times, log_times])
        if tlen != args.start:
            times = np.delete(times, 0)

        snodes, p0 = TL.get_occupancies()
        assert np.isclose(sum(p0), 1)
        if args.plot_minh:
            assert args.plot_minh > TL.minh, "Plot-minh must be greater than minh."
            # Provide coarse network to get even more coarse network
            ndata = {n: d for n, d in TL.nodes.items() if n in snodes} # only active.
            edata = TL.cg_edges
            assert (n in snodes for n in ndata)
            assert (n in ndata for n in snodes)
            _, _, plot_cgm = top_down_coarse_graining(ndata, edata, args.plot_minh)

            mapping = {hn: set() for hn in snodes}
            for lmin, hidden in plot_cgm.items():
                assert lmin in snodes
                for hn in hidden:
                    assert hn in snodes
                    mapping[hn].add(lmin)
        else:
            plot_cgm, mapping = None, None

        if args.tmpdir:
            _fname = _tmpdir + '/' + name + '-' + str(tlen)
            # Produce input for treekin simulation
            nlist, bbfile, brfile, bofile, callfile = get_tkn_simulation_files(TL, _fname) 

        ti, p8, pf = 0, np.zeros(len(p0)), np.zeros(len(p0))
        for (t, pt) in TL.simulate(snodes, p0, times):
            if tlen < args.stop - 1 and t > t1:
                # NOTE: this is extra time to determine whether structures will
                # become important in the timeframe of transcription.
                pf = np.maximum(pf, pt)
                continue
            if args.stdout == 'drf' or dfh:
                tt = time + t
                for e, occu in enumerate(pt):
                    node = snodes[e]
                    if plot_cgm and node not in plot_cgm:
                        continue
                    elif plot_cgm:
                        for n in plot_cgm[node]:
                            ci = snodes.index(n)
                            occu += pt[ci]/len(mapping[n])
                    ni = TL.nodes[node]['identity']
                    if occu < 0.001 and ni not in all_courses:
                        continue
                    ne = TL.nodes[node]['energy']/100
                    all_courses[ni] = [(tt, occu)]
                    fdata = f"{ni} {tt:03.9f} {occu:03.4f} {node[:tlen]} {ne:6.2f}\n"
                    write_output(fdata, stdout = (args.stdout == 'drf'), fh = dfh)
            ti, p8 = t, pt
        keep = [n for e, n in enumerate(snodes) if pf[e] > args.o_prune] if args.o_prune else []
        TL.set_occupancies(snodes, p8)
        TL.total_time += ti

        # ~~~~~~~~~~~~~~~~~~~~~~~~~ #
        # Output simulation results #
        # ~~~~~~~~~~~~~~~~~~~~~~~~~ #

        # Print the current state *after* the simulation.
        if args.stdout == 'log' or lfh:
            for e, node in enumerate(snodes):
                ni = TL.nodes[node]['identity']
                ne = TL.nodes[node]['energy']/100
                po = p0[e]
                no = p8[e]
                if args.plot_minh:
                    lmins = [TL.nodes[n]['identity'] for n in sorted(mapping[node], 
                                                                     key = lambda x: TL.nodes[x]['identity'])]
                    if lmins:
                        ni +=  f" -> {', '.join(lmins)}"
                fdata = f"{tlen:4d} {e+1:4d} {node[:tlen]} {ne:6.2f} [{po:6.4f} -> {no:6.4f}] ID = {ni}\n"
                write_output(fdata, stdout = (args.stdout == 'log'), fh = lfh)
        stime = datetime.now()

        # ~~~~~ #
        # Prune #
        # ~~~~~ #
        if args.o_prune > 0:
            delth = 10 # delete structures if inactive for more than 10 rounds.
            pn, dn = TL.prune(args.o_prune, delth, keep) 
            logger.info(f'After pruning:         {len(list(TL.active_local_mins)):3d} active lmins, {len(list(TL.active_nodes)):3d} active structures, {len(list(TL.inactive_nodes)):3d} inactive structures.' +
                    f' (Pruned {len(pn)} nodes with combined occupancy below {args.o_prune:.2f}, kept {len(keep)} nodes due to lookahead, deleted {len(dn)} nodes inactive for {delth} transcription steps.)')

        ptime = datetime.now()
        if args.verbose:
            exptime = (etime - itime).total_seconds() 
            cgntime = (ctime - etime).total_seconds()
            simtime = (stime - ctime).total_seconds()
            prntime = (ptime - stime).total_seconds()
            tottime = (ptime - itime).total_seconds()
            logger.info(f'{tlen=}, {tottime=}, {exptime=}, {cgntime=}, {simtime=}, {prntime=}.')
            #print(tlen, tottime, exptime, cgntime, simtime, prntime)
            #sys.stdout.flush()

        if args.draw_graphs:
            raise NotImplementedError('Cannot draw graphs right now.')
            TL.graph_to_json(_fname)

    # Write the last results
    if args.stdout == 'log' or lfh:
        fdata  = "# Distribution of structures at the end:\n"
        fdata += "#         {}\n".format(TL.transcript)
        if args.plot_minh:
            for e, node in enumerate(snodes):
                if node not in plot_cgm:
                    continue
                ne = TL.nodes[node]['energy']/100
                no = p8[e] + sum(p8[snodes.index(n)]/len(mapping[n]) for n in plot_cgm[node])
                ni = TL.nodes[node]['identity']
                nids = [TL.nodes[n]['identity'] for n in sorted(plot_cgm[node], key = lambda x: TL.nodes[x]['identity'])]
                if nids:
                    ni +=  f" + {' + '.join(nids)}"
                fdata += f"{tlen:4d} {e+1:4d} {node[:tlen]} {ne:6.2f} {no:6.4f} ID = {ni}\n"
            write_output(fdata, stdout = (args.stdout == 'log'), fh = lfh)
        else:
            for e, node in enumerate(sorted(TL.active_local_mins, key = lambda x: TL.nodes[x]['energy'])):
                ne = TL.nodes[node]['energy']/100
                no = p8[e]
                ni = TL.nodes[node]['identity']
                fdata += f"{tlen:4d} {e:4d} {node[:tlen]} {ne:6.2f} {no:6.4f} ID = {ni}\n"
            write_output(fdata, stdout=(args.stdout == 'log'), fh = lfh)

    #statprof.stop()
    #statprof.display()

    # CLEANUP file handles
    if lfh: lfh.close()
    if dfh: dfh.close()
    return

if __name__ == '__main__':
    main()
